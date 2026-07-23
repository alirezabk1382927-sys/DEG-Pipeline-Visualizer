"""
DEG Pipeline & Visualizer — Multi-Cancer Analysis
==================================================
Differential Expression Analysis Tool with publication-quality plots.

Steps:
  Step 1: Read raw count files + sample sheet → PyDESeq2 DEG analysis
  Step 2: Generate Volcano Plot, MA Plot, Bar Chart, Heatmap (with optional gene labeling)

Author:  Alireza Balaei
GitHub:  https://github.com/alirezabk1382927-sys
LinkedIn: https://ir.linkedin.com/in/alireza-balaei-kahnamoei-aa8216344

--------------------------------------------------------------------------
IMPORTANT — PyInstaller / frozen build notice
--------------------------------------------------------------------------
PyDESeq2 uses joblib internally. When the app is frozen into a single .exe,
joblib workers can accidentally re-launch the full GUI. We prevent this by
setting LOKY_MAX_CPU_COUNT=1 before any imports, and calling
multiprocessing.freeze_support() at the entry point.
"""

import os
import sys
import re
import gzip
import queue
import threading
import traceback
import multiprocessing
import subprocess
import webbrowser
from pathlib import Path
from datetime import datetime

# ==================================================================
# 0) Frozen-executable detection + joblib/loky fix
# ==================================================================
IS_FROZEN = getattr(sys, "frozen", False)

if IS_FROZEN:
    SAFE_N_CPUS = 1
else:
    SAFE_N_CPUS = max(1, min(4, (os.cpu_count() or 1)))

os.environ["LOKY_MAX_CPU_COUNT"] = str(SAFE_N_CPUS)
os.environ.setdefault("OMP_NUM_THREADS", str(SAFE_N_CPUS))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(SAFE_N_CPUS))
os.environ.setdefault("NUMBA_NUM_THREADS", str(SAFE_N_CPUS))

# ==================================================================
# 1) Dependency check
# ==================================================================
REQUIRED_PACKAGES = [
    ("pandas", "pandas"),
    ("numpy", "numpy"),
    ("matplotlib", "matplotlib"),
    ("seaborn", "seaborn"),
    ("adjustText", "adjustText"),
    ("pydeseq2", "pydeseq2"),
]


def _pip_install(pip_spec: str) -> bool:
    cmd = [sys.executable, "-m", "pip", "install", "--quiet", pip_spec]
    kwargs = {}
    if os.name == "nt":
        kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
    try:
        subprocess.check_call(cmd, timeout=300, **kwargs)
        return True
    except Exception:
        return False


def _ensure_package(import_name: str, pip_spec: str) -> bool:
    import importlib
    try:
        importlib.import_module(import_name)
        return True
    except ImportError:
        pass

    if IS_FROZEN:
        return False

    if _pip_install(pip_spec):
        importlib.invalidate_caches()
        try:
            importlib.import_module(import_name)
            return True
        except ImportError:
            return False
    return False


_missing_packages = []
for _imp, _spec in REQUIRED_PACKAGES:
    if not _ensure_package(_imp, _spec):
        _missing_packages.append(_imp)

if _missing_packages:
    _msg = (
        "The following required packages are missing: "
        f"{', '.join(_missing_packages)}.\n\n"
        "If you are running the packaged .exe, this build is missing "
        "dependencies and needs to be rebuilt (see README).\n"
        "If you are running from source, install them with:\n"
        f"    pip install {' '.join(_missing_packages)}"
    )
    if IS_FROZEN:
        try:
            import tkinter as tk
            from tkinter import messagebox
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror("Missing dependencies", _msg)
        except Exception:
            pass
        sys.exit(1)
    else:
        print(_msg)
        sys.exit(1)

# ==================================================================
# Now safe to import the rest
# ==================================================================
import numpy as np
import pandas as pd
import matplotlib

try:
    matplotlib.use("TkAgg")
except ImportError:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import seaborn as sns

try:
    from adjustText import adjust_text  # noqa: F401
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

APP_AUTHOR = "Alireza Balaei"
APP_GITHUB = "https://github.com/alirezabk1382927-sys"
APP_LINKEDIN = "https://ir.linkedin.com/in/alireza-balaei-kahnamoei-aa8216344"
APP_VERSION = "2.0.2"


# ==================================================================
# Helper: icon path
# ==================================================================
def get_icon_path():
    if getattr(sys, "frozen", False):
        base_path = Path(sys._MEIPASS)
        for name in ("icon.ico", "icon.png"):
            test_path = base_path / name
            if test_path.exists():
                return test_path
        exe_path = Path(sys.executable).parent / "icon.ico"
        if exe_path.exists():
            return exe_path
    else:
        base_path = Path(__file__).parent

    for path in (
        base_path / "icon.ico",
        base_path / "icon.png",
        Path(sys.executable).parent / "icon.ico",
        Path.cwd() / "icon.ico",
        Path.cwd() / "icon.png",
    ):
        if path.exists():
            return path
    return None


# ==================================================================
# 2) Publication style
# ==================================================================
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 12,
    "axes.linewidth": 0.9,
    "axes.edgecolor": "#4d4d4d",
    "axes.titleweight": "bold",
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "savefig.dpi": 600,
    "figure.dpi": 105,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

CUTOFF_COLOR = "#2E7D32"
GRID_COLOR = "#e2e2e2"
PALETTE = {"Up": "#FF0000", "Down": "#0000FF", "NS": "#B3B3B3"}
PALETTE_ALT = {"Up": "#d62828", "Down": "#1d3557", "NS": "#c9c9c9"}

STD_NAMES = {
    "all": "All_DEGs_results.csv",
    "cache": "analysis_cache.csv",
    "sig": "Significant_DEGs.csv",
    "up": "Upregulated_in_Tumor.csv",
    "down": "Downregulated_in_Tumor.csv",
    "norm": "normalized_counts.csv",
    "groups": "sample_grouping.csv",
    "dupes": "duplicate_samples_dropped.csv",
    "log": "analysis_log.txt",
    "gene_list": "all_genes_complete_list.csv",
}


def apply_journal_style(ax):
    ax.set_facecolor("white")
    ax.grid(True, color=GRID_COLOR, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#4d4d4d")
        spine.set_linewidth(0.9)


def get_output_dirs(base_dir: Path):
    if base_dir is None:
        raise ValueError("Output directory not set.")
    base_dir = Path(base_dir)
    dirs = {
        "root": base_dir,
        "data": base_dir / "Data",
        "pictures": base_dir / "Pictures",
        "logs": base_dir / "Logs",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


# ==================================================================
# 3) STEP 1 — raw count files -> PyDESeq2 DEG analysis
# ==================================================================
def read_one_count_file(fpath: Path):
    """Returns (count_series, gene_name_series) or None on failure."""
    if not fpath.exists():
        return None

    open_func = gzip.open if fpath.suffix == ".gz" else open
    try:
        with open_func(fpath, "rt", encoding="utf-8", errors="replace") as f:
            raw_lines = []
            for _ in range(50):
                line = f.readline()
                if not line:
                    break
                raw_lines.append(line)
    except Exception:
        return None

    header_row_idx = None
    for idx, line in enumerate(raw_lines):
        if "gene_id" in line.lower():
            header_row_idx = idx
            break
    if header_row_idx is None:
        return None

    try:
        df = pd.read_csv(fpath, sep="\t", dtype=str, skiprows=header_row_idx)
    except Exception:
        return None

    df.columns = [c.strip() for c in df.columns]
    if "gene_id" not in df.columns or "unstranded" not in df.columns:
        return None
    if "gene_name" not in df.columns:
        df["gene_name"] = df["gene_id"]

    df = df[~df["gene_id"].astype(str).str.startswith(("N_", "__"))]
    dups = df.duplicated(subset="gene_id", keep="first")
    if dups.any():
        df = df[~dups]

    counts = pd.to_numeric(df["unstranded"], errors="coerce")
    nan_mask = counts.isna()
    if nan_mask.any():
        df = df[~nan_mask]
        counts = counts[~nan_mask]

    counts.index = df["gene_id"].values
    return counts, df.set_index("gene_id")["gene_name"]


def normalize_group(raw_series: pd.Series, log):
    def _map(v):
        if pd.isna(v):
            return v
        s = str(v).strip().lower()
        if "primary tumor" in s or "tumor" in s or "tumour" in s:
            return "Tumor"
        if "solid tissue normal" in s or "normal" in s:
            return "Normal"
        return "Other"

    mapped = raw_series.map(_map)
    log("Sample groups: Tumor={}, Normal={}, Other={}".format(
        (mapped == "Tumor").sum(),
        (mapped == "Normal").sum(),
        (mapped == "Other").sum(),
    ))
    return mapped


def find_case_id_column(columns):
    cols_lower = {c.lower().strip(): c for c in columns}
    for cand in ["case id", "case ids", "caseid", "case_id"]:
        if cand in cols_lower:
            return cols_lower[cand]
    for c, orig in cols_lower.items():
        if "case" in c:
            return orig
    return None


def _build_inference(n_cpus: int, log):
    try:
        from pydeseq2.default_inference import DefaultInference
        return DefaultInference(n_cpus=n_cpus), {}
    except Exception as e:
        log(f"DefaultInference not available ({e}); using n_cpus fallback.")
        return None, {"n_cpus": n_cpus}


def run_deg_pipeline(counts_root, sample_sheet_path, output_dir, label,
                      min_count=10, log2fc_thr=1.0, padj_thr=0.05,
                      log=print, progress_callback=None):
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as e:
        raise RuntimeError(f"PyDESeq2 is not installed. Error: {e}")

    counts_root = Path(counts_root)
    sample_sheet_path = Path(sample_sheet_path)
    dirs = get_output_dirs(output_dir)
    data_dir = dirs["data"]

    def update_progress(pct, msg=""):
        if progress_callback:
            progress_callback(pct, msg)

    update_progress(0, "Reading sample sheet...")
    log(f"Reading sample sheet: {sample_sheet_path}")
    sample_sheet = pd.read_csv(sample_sheet_path, sep="\t", dtype=str)
    sample_sheet.columns = [c.strip() for c in sample_sheet.columns]

    need_cols_map = {"File ID": "File ID", "File Name": "File Name",
                     "Tissue Type": "Tissue Type", "Sample ID": "Sample ID"}
    sheet_cols_lower = {c.lower(): c for c in sample_sheet.columns}
    missing = []
    for k, v in need_cols_map.items():
        if v.lower() not in sheet_cols_lower:
            missing.append(v)
        else:
            sample_sheet.rename(columns={sheet_cols_lower[v.lower()]: v}, inplace=True)
    if missing:
        raise ValueError(f"Missing columns from sample sheet: {missing}")

    case_id_col = find_case_id_column(sample_sheet.columns)
    if case_id_col is None:
        log("[WARN] No 'Case ID' column found. Using 'Sample ID' as Case ID.")
        sample_sheet["Case ID"] = sample_sheet["Sample ID"]
    else:
        if case_id_col != "Case ID":
            sample_sheet.rename(columns={case_id_col: "Case ID"}, inplace=True)

    sample_sheet["Sample ID"] = sample_sheet["Sample ID"].astype(str).str.strip()
    log(f"Sample sheet rows: {len(sample_sheet)}")

    dup_mask = sample_sheet.duplicated(subset="Sample ID", keep="first")
    if dup_mask.any():
        dropped = sample_sheet[dup_mask].copy()
        dropped.to_csv(data_dir / f"{label}_{STD_NAMES['dupes']}", index=False)
        log(f"[WARN] {dup_mask.sum()} duplicate Sample IDs found. Keeping first occurrence only.")
        sample_sheet = sample_sheet[~dup_mask].reset_index(drop=True)

    log(f"Reading {len(sample_sheet)} count files...")

    # Build count matrix efficiently using dictionary of arrays
    count_dict = {}
    gene_name_map = {}
    valid_samples = []
    failed_files = []
    all_gene_ids = set()

    total = len(sample_sheet)
    for i, row in sample_sheet.iterrows():
        fpath = counts_root / row["File ID"] / row["File Name"]
        result = read_one_count_file(fpath)
        if result is None:
            failed_files.append(row["Sample ID"])
            continue
        counts, gnames = result
        # Store as plain dict for fast alignment later
        count_dict[row["Sample ID"]] = dict(zip(counts.index, counts.values))
        all_gene_ids.update(counts.index)
        for gid, gname in gnames.items():
            if gid not in gene_name_map:
                gene_name_map[gid] = gname
        valid_samples.append(row["Sample ID"])
        pct = 5 + int((i + 1) / total * 15)  # 5% to 20%
        if (i + 1) % 50 == 0 or (i + 1) == total:
            update_progress(pct, f"Reading files: {i + 1}/{total}")
            log(f"  -> {i + 1}/{total} files read")

    if failed_files:
        log(f"[WARN] {len(failed_files)} files failed and were skipped.")
    if not valid_samples:
        raise RuntimeError("No count files were read successfully.")

    all_gene_ids = sorted(all_gene_ids)
    log(f"Total unique gene IDs: {len(all_gene_ids)}")

    # Build count matrix efficiently using DataFrame from dict
    update_progress(20, "Building count matrix...")
    log("Building count matrix...")

    # Create DataFrame directly from the dictionaries
    count_matrix = pd.DataFrame.from_dict(count_dict, orient='index').T
    count_matrix = count_matrix.reindex(all_gene_ids).fillna(0).astype(int)
    count_matrix = count_matrix[valid_samples]

    # Build metadata
    sample_sheet_valid = sample_sheet[sample_sheet["Sample ID"].isin(valid_samples)].copy()
    metadata = sample_sheet_valid[["Sample ID", "Case ID", "Tissue Type"]].copy()
    metadata = metadata.set_index("Sample ID")
    metadata.rename(columns={"Tissue Type": "Group_raw"}, inplace=True)
    metadata["Group"] = normalize_group(metadata["Group_raw"], log)

    # Keep only Tumor/Normal
    metadata = metadata[metadata["Group"].isin(["Tumor", "Normal"])]
    common_samples = list(set(count_matrix.columns) & set(metadata.index))
    count_matrix = count_matrix[common_samples]
    metadata = metadata.loc[common_samples]

    log(f"Final matrix: {count_matrix.shape[1]} samples x {count_matrix.shape[0]} genes")

    group_counts = metadata["Group"].value_counts()
    if len(group_counts) < 2 or group_counts.min() < 2:
        raise ValueError(f"Need at least 2 samples in BOTH groups. Current: {group_counts.to_dict()}")

    metadata.to_csv(data_dir / f"{label}_{STD_NAMES['groups']}")

    # Filter low-count genes
    update_progress(25, "Filtering low-count genes...")
    min_group_size = metadata["Group"].value_counts().min()
    keep_genes = (count_matrix >= min_count).sum(axis=1) >= min_group_size
    log(f"Genes before filter: {count_matrix.shape[0]} | after: {keep_genes.sum()}")
    count_matrix_filtered = count_matrix.loc[keep_genes]

    # Prepare for DESeq2
    counts_for_deseq = count_matrix_filtered.T
    counts_for_deseq.index = counts_for_deseq.index.astype(str)
    metadata.index = metadata.index.astype(str)
    metadata_deseq = metadata[["Group"]].copy()

    if not counts_for_deseq.index.equals(metadata_deseq.index):
        common_idx = counts_for_deseq.index.intersection(metadata_deseq.index)
        counts_for_deseq = counts_for_deseq.loc[common_idx]
        metadata_deseq = metadata_deseq.loc[common_idx]

    update_progress(30, "Running PyDESeq2...")
    log(f"Running PyDESeq2 (capped at {SAFE_N_CPUS} CPU core(s))...")

    inference, extra_kwargs = _build_inference(SAFE_N_CPUS, log)
    dds_kwargs = dict(
        counts=counts_for_deseq,
        metadata=metadata_deseq,
        design_factors="Group",
        refit_cooks=True,
    )
    if inference is not None:
        dds_kwargs["inference"] = inference
    else:
        dds_kwargs.update(extra_kwargs)

    dds = DeseqDataSet(**dds_kwargs)
    update_progress(40, "DESeq2 fitting...")
    dds.deseq2()

    update_progress(70, "Computing statistics...")
    stats_kwargs = dict(contrast=["Group", "Tumor", "Normal"])
    if inference is not None:
        stats_kwargs["inference"] = inference
    else:
        stats_kwargs.update(extra_kwargs)
    stat_res = DeseqStats(dds, **stats_kwargs)
    stat_res.summary()

    update_progress(85, "Saving results...")
    res = stat_res.results_df.copy().reset_index().rename(columns={"index": "gene_id"})
    if gene_name_map:
        res["gene_name"] = res["gene_id"].map(gene_name_map)
    else:
        res["gene_name"] = res["gene_id"]
    res = res.sort_values("padj")

    res_clean = res.dropna(subset=["padj"])
    degs = res_clean[(res_clean["log2FoldChange"].abs() >= log2fc_thr) &
                     (res_clean["padj"] < padj_thr)]
    degs_up = degs[degs["log2FoldChange"] > 0]
    degs_down = degs[degs["log2FoldChange"] < 0]

    log("=========== Summary ===========")
    log(f"Total genes tested: {len(res_clean)}")
    log(f"DEGs (|log2FC|>={log2fc_thr}, padj<{padj_thr}): {len(degs)}")
    log(f"  Upregulated: {len(degs_up)}")
    log(f"  Downregulated: {len(degs_down)}")

    paths = {
        "all": data_dir / f"{label}_{STD_NAMES['all']}",
        "cache": data_dir / f"{label}_{STD_NAMES['cache']}",
        "sig": data_dir / f"{label}_{STD_NAMES['sig']}",
        "up": data_dir / f"{label}_{STD_NAMES['up']}",
        "down": data_dir / f"{label}_{STD_NAMES['down']}",
        "norm": data_dir / f"{label}_{STD_NAMES['norm']}",
        "groups": data_dir / f"{label}_{STD_NAMES['groups']}",
    }
    res.to_csv(paths["all"], index=False)
    res.to_csv(paths["cache"], index=False)
    degs.to_csv(paths["sig"], index=False)
    degs_up.to_csv(paths["up"], index=False)
    degs_down.to_csv(paths["down"], index=False)

    normalized_counts = pd.DataFrame(
        dds.layers["normed_counts"],
        index=counts_for_deseq.index,
        columns=counts_for_deseq.columns,
    )
    normalized_counts.to_csv(paths["norm"])

    gene_list = res_clean[["gene_id", "gene_name", "log2FoldChange",
                           "padj", "pvalue", "baseMean"]].copy()
    gene_list["regulation"] = "NS"
    gene_list.loc[degs_up.index, "regulation"] = "Up"
    gene_list.loc[degs_down.index, "regulation"] = "Down"
    gene_list = gene_list.sort_values("padj")
    gene_list_path = data_dir / f"{label}_{STD_NAMES['gene_list']}"
    gene_list.to_csv(gene_list_path, index=False)
    paths["gene_list"] = gene_list_path

    update_progress(100, "Complete!")
    log(f"\nAll result files saved to: {data_dir}")
    return paths


# ==================================================================
# 4) DEG table prep + gene matching + plotting
# ==================================================================
def load_csv(path):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def prep_deg_table(df, log2fc_thr, padj_thr):
    if {"gene", "log2FC", "sig"}.issubset(df.columns):
        out = df[["gene", "log2FC", "sig"]].copy()
    else:
        lfc_col = pick_col(df, ["log2FoldChange", "log2FC", "logFC"])
        padj_col = pick_col(df, ["padj", "FDR", "adj.P.Val", "qvalue"])
        name_col = pick_col(df, ["gene_name", "gene_symbol", "hgnc_symbol", "gene_id"])
        if lfc_col is None:
            raise ValueError("Could not find a log2FoldChange column.")
        if padj_col is None:
            raise ValueError("Could not find a padj/pvalue column.")
        if name_col is None:
            raise ValueError("Could not find a gene name column.")
        out = pd.DataFrame({
            "gene": df[name_col].astype(str),
            "log2FC": pd.to_numeric(df[lfc_col], errors="coerce"),
            "sig": pd.to_numeric(df[padj_col], errors="coerce"),
        })

    out["log2FC"] = pd.to_numeric(out["log2FC"], errors="coerce")
    out["sig"] = pd.to_numeric(out["sig"], errors="coerce")
    out = out.dropna(subset=["log2FC", "sig"])
    out["neglog10"] = -np.log10(out["sig"].clip(lower=1e-300))
    out["class"] = "NS"
    out.loc[(out["log2FC"] >= log2fc_thr) & (out["sig"] < padj_thr), "class"] = "Up"
    out.loc[(out["log2FC"] <= -log2fc_thr) & (out["sig"] < padj_thr), "class"] = "Down"
    return out


def parse_gene_list(raw_text: str):
    parts = re.split(r"[,\n;]+", raw_text)
    seen, out = set(), []
    for p in parts:
        g = p.strip()
        if g and g not in seen:
            seen.add(g)
            out.append(g)
    return out


def find_matching_genes(deg_df, requested):
    if "gene" not in deg_df.columns:
        raise ValueError("DataFrame must have a 'gene' column for labeling.")
    gene_series = deg_df["gene"].astype(str)
    exact_set = set(gene_series)
    lower_map = {}
    for g in gene_series:
        if g.lower() not in lower_map:
            lower_map[g.lower()] = g

    found, not_found = [], []
    for req in requested:
        if req in exact_set:
            found.append(req)
        elif req.lower() in lower_map:
            found.append(lower_map[req.lower()])
        else:
            not_found.append(req)
    return found, not_found


def save_processed_data_and_summary(deg_df, data_dir, base_name, log2fc_thr, padj_thr):
    data_path = Path(data_dir) / f"{base_name}_processed_data.csv"
    deg_df.to_csv(data_path, index=False)
    up = deg_df[deg_df["class"] == "Up"].sort_values("sig")
    down = deg_df[deg_df["class"] == "Down"].sort_values("sig")

    lines = [
        f"Analysis summary - {base_name}",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Thresholds: |log2FC| >= {log2fc_thr}, padj < {padj_thr}",
        "",
        f"Total genes analyzed: {len(deg_df)}",
        f"Upregulated: {len(up)}",
        f"Downregulated: {len(down)}",
        f"Not significant: {(deg_df['class']=='NS').sum()}",
        "",
        "Top 10 upregulated (by significance):",
    ]
    for _, r in up.head(10).iterrows():
        lines.append(f"  {r['gene']:<15} log2FC={r['log2FC']:.3f}  sig={r['sig']:.3e}")
    lines.append("")
    lines.append("Top 10 downregulated (by significance):")
    for _, r in down.head(10).iterrows():
        lines.append(f"  {r['gene']:<15} log2FC={r['log2FC']:.3f}  sig={r['sig']:.3e}")

    summary_path = Path(data_dir) / f"{base_name}_summary.txt"
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    return data_path, summary_path


# ==================================================================
# 5) Plotting functions
# ==================================================================
def add_gene_labels_stacked(ax, points_df, x_col, y_col, gene_col, log,
                            side="right", fontsize=6):
    if points_df.empty:
        return

    points = points_df.sort_values(by=y_col, ascending=False).reset_index(drop=True)
    n = len(points)

    for _, row in points.iterrows():
        ax.scatter(row[x_col], row[y_col], s=60, facecolors="none",
                   edgecolors="black", linewidths=1.2, zorder=10)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    if side == "right":
        label_x = xmax - 0.05 * (xmax - xmin)
    else:
        label_x = xmin + 0.05 * (xmax - xmin)

    margin = 0.03 * (ymax - ymin)
    usable_height = (ymax - ymin) - 2 * margin
    if n > 1:
        dy = min(0.04 * (ymax - ymin), usable_height / (n - 1))
    else:
        dy = 0
    total_height = (n - 1) * dy
    y_start = ymax - margin - (usable_height - total_height) / 2
    if n == 1:
        y_start = (ymax + ymin) / 2

    for i, (_, row) in enumerate(points.iterrows()):
        px, py = row[x_col], row[y_col]
        ly = y_start - i * dy

        ax.plot([px, label_x], [py, ly], color="black", linewidth=0.6, zorder=9,
                clip_on=False)

        ax.text(label_x, ly, str(row[gene_col]),
                fontsize=fontsize, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor="black", alpha=0.9),
                zorder=11, ha="left" if side == "right" else "right",
                va="center", clip_on=False)

    log(f"Labeled {n} gene(s).")


def plot_volcano(ax, deg_df, log2fc_thr, padj_thr, palette, title,
                 label_genes=None, log_func=print):
    label_genes = label_genes or []
    order = ["Down", "NS", "Up"]
    display_names = {"Down": "Down-regulated", "NS": "Not Significant", "Up": "Up-regulated"}
    for cls in order:
        sub = deg_df[deg_df["class"] == cls]
        ax.scatter(sub["log2FC"], sub["neglog10"], s=16, c=palette[cls], alpha=0.75,
                    linewidths=0, label=display_names[cls], zorder=3)
    ax.axvline(log2fc_thr, ls="--", lw=1.2, color=CUTOFF_COLOR, zorder=2)
    ax.axvline(-log2fc_thr, ls="--", lw=1.2, color=CUTOFF_COLOR, zorder=2)
    ax.axhline(-np.log10(padj_thr), ls="--", lw=1.2, color=CUTOFF_COLOR, zorder=2)
    ax.set_xlabel("Log2 Fold Change", fontsize=13)
    ax.set_ylabel("-Log10(FDR)", fontsize=13)
    ax.set_title(title, loc="left", fontsize=16, fontweight="bold", pad=12)
    apply_journal_style(ax)
    ax.legend(title="Status", loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
              fontsize=11, title_fontsize=12, markerscale=1.8, frameon=False)

    if label_genes:
        targets = deg_df[deg_df["gene"].astype(str).isin(label_genes)]
        if not targets.empty:
            add_gene_labels_stacked(ax, targets, "log2FC", "neglog10", "gene",
                                    log_func, side="right", fontsize=6)


def plot_ma(ax, df_raw, log2fc_thr, padj_thr, palette, title,
            label_genes=None, log_func=print):
    lfc_col = pick_col(df_raw, ["log2FoldChange", "log2FC", "logFC"])
    padj_col = pick_col(df_raw, ["padj", "FDR", "adj.P.Val", "qvalue"])
    mean_col = pick_col(df_raw, ["baseMean", "AveExpr", "meanExpr"])
    if lfc_col is None or padj_col is None or mean_col is None:
        raise ValueError("Required columns missing for MA plot.")

    d = pd.DataFrame({
        "mean": pd.to_numeric(df_raw[mean_col], errors="coerce"),
        "log2FC": pd.to_numeric(df_raw[lfc_col], errors="coerce"),
        "padj": pd.to_numeric(df_raw[padj_col], errors="coerce"),
    })
    d = d[d["mean"] > 0].dropna()
    d["class"] = "NS"
    d.loc[(d["log2FC"] >= log2fc_thr) & (d["padj"] < padj_thr), "class"] = "Up"
    d.loc[(d["log2FC"] <= -log2fc_thr) & (d["padj"] < padj_thr), "class"] = "Down"
    order = ["Down", "NS", "Up"]
    display_names = {"Down": "Down-regulated", "NS": "Not Significant", "Up": "Up-regulated"}
    for cls in order:
        sub = d[d["class"] == cls]
        ax.scatter(sub["mean"], sub["log2FC"], s=14, c=palette[cls], alpha=0.75,
                    linewidths=0, label=display_names[cls], zorder=3)
    ax.set_xscale("log")
    ax.axhline(0, color=CUTOFF_COLOR, lw=1.2, ls="--", zorder=2)
    ax.set_xlabel("Mean of normalized counts", fontsize=13)
    ax.set_ylabel("Log2 Fold Change", fontsize=13)
    ax.set_title(title, loc="left", fontsize=16, fontweight="bold", pad=12)
    apply_journal_style(ax)
    ax.legend(title="Status", loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
              fontsize=11, title_fontsize=12, markerscale=1.8, frameon=False)

    if label_genes:
        gene_col = pick_col(df_raw, ["gene_name", "gene_id"])
        if gene_col is None:
            return
        df_pos = df_raw.loc[d.index]
        gene_names = df_pos[gene_col].astype(str)
        mask = gene_names.isin(label_genes)
        targets = df_pos.loc[mask]
        if not targets.empty:
            targets_df = pd.DataFrame({
                "gene": targets[gene_col],
                "mean": targets[mean_col],
                "log2FC": targets[lfc_col],
            })
            add_gene_labels_stacked(ax, targets_df, "mean", "log2FC", "gene",
                                    log_func, side="right", fontsize=6)


def plot_summary_bar(ax, up_n, down_n, title):
    cats, vals = ["Upregulated", "Downregulated"], [up_n, down_n]
    colors = [PALETTE["Up"], PALETTE["Down"]]
    bars = ax.bar(cats, vals, color=colors, width=0.5, edgecolor="#333333", linewidth=1, zorder=3)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center", va="bottom",
                fontsize=12, fontweight="bold")
    ax.set_ylabel("Number of genes", fontsize=13)
    ax.set_title(title, loc="left", fontsize=16, fontweight="bold", pad=12)
    apply_journal_style(ax)


def plot_heatmap_expression(ax, norm_counts, sample_groups, gene_ids, gene_labels, title):
    present_ids = [g for g in gene_ids if g in norm_counts.index]
    if not present_ids:
        ax.text(0.5, 0.5, "No matching genes", transform=ax.transAxes, ha="center")
        return
    mat = norm_counts.loc[present_ids].copy()
    mat = np.log2(mat.clip(lower=0) + 1)
    row_means = mat.mean(axis=1)
    row_stds = mat.std(axis=1)
    row_stds = row_stds.replace(0, 1)
    z = mat.sub(row_means, axis=0).div(row_stds, axis=0)
    order = [s for s in sample_groups.sort_values("Group").index if s in z.columns]
    z = z[order]
    row_labels = [gene_labels.get(g, g) for g in z.index]
    sns.heatmap(z, ax=ax, cmap="RdBu_r", center=0, cbar_kws={"label": "Z-score"},
                yticklabels=row_labels, xticklabels=False)
    ax.set_title(title, loc="left", fontsize=15, fontweight="bold", pad=10)
    ax.set_xlabel(f"Samples (n={len(order)})", fontsize=12)


def plot_heatmap_rank(ax, deg_df, top_n, title):
    d = deg_df[deg_df["class"] != "NS"].copy()
    if d.empty:
        ax.text(0.5, 0.5, "No significant genes", transform=ax.transAxes, ha="center")
        return
    d = d.reindex(d["log2FC"].abs().sort_values(ascending=False).index).head(top_n).sort_values("log2FC")
    colors = [PALETTE["Up"] if v > 0 else PALETTE["Down"] for v in d["log2FC"]]
    ax.barh(d["gene"], d["log2FC"], color=colors, edgecolor="#333333", linewidth=0.6, zorder=3)
    ax.axvline(0, color="#333333", lw=1)
    ax.set_xlabel("Log2 Fold Change", fontsize=13)
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=10)
    apply_journal_style(ax)


PLOT_DESCRIPTIONS = {
    "Volcano Plot": (
        "Every gene's fold-change (x) vs. significance (y).\n"
        "Load the cache file (analysis_cache.csv or All_DEGs_results.csv).\n"
        "List genes in the box below to auto-save a labeled version."
    ),
    "MA Plot": (
        "Fold-change vs. average expression level.\n"
        "Requires a file with baseMean column - use cache file, not filtered files."
    ),
    "Summary Bar Chart": (
        "Bar chart of upregulated vs. downregulated gene counts."
    ),
    "Heatmap (Top DEGs)": (
        "Expression heatmap of top DEGs across samples.\n"
        "Best with: Significant_DEGs.csv + normalized_counts.csv + sample_grouping.csv.\n"
        "Without those, a fold-change rank heatmap is drawn instead."
    ),
}


# ==================================================================
# 6) Scrollable frame helper
# ==================================================================
class ScrollableFrame(ttk.Frame):
    def __init__(self, parent, *a, **kw):
        super().__init__(parent, *a, **kw)
        self.canvas = tk.Canvas(self, highlightthickness=0)
        vscroll = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)

        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(self._win, width=e.width))
        self.canvas.configure(yscrollcommand=vscroll.set)

        self.canvas.pack(side="left", fill="both", expand=True)
        vscroll.pack(side="right", fill="y")

        self.canvas.bind("<Enter>", self._bind_wheel)
        self.canvas.bind("<Leave>", self._unbind_wheel)

    def _on_wheel(self, event):
        delta = -1 * (event.delta // 120) if event.delta else (1 if event.num == 5 else -1)
        self.canvas.yview_scroll(int(delta), "units")

    def _bind_wheel(self, _event):
        self.canvas.bind_all("<MouseWheel>", self._on_wheel)
        self.canvas.bind_all("<Button-4>", self._on_wheel)
        self.canvas.bind_all("<Button-5>", self._on_wheel)

    def _unbind_wheel(self, _event):
        self.canvas.unbind_all("<MouseWheel>")
        self.canvas.unbind_all("<Button-4>")
        self.canvas.unbind_all("<Button-5>")


# ==================================================================
# 7) Welcome dialog with GitHub & LinkedIn buttons
# ==================================================================
class WelcomeDialog(tk.Toplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title("Welcome")
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()

        self.update_idletasks()
        pw, ph = parent.winfo_width(), parent.winfo_height()
        px, py = parent.winfo_x(), parent.winfo_y()
        w, h = 460, 240
        x = px + (pw - w) // 2
        y = py + (ph - h) // 2
        self.geometry(f"{w}x{h}+{x}+{y}")

        frame = ttk.Frame(self, padding=20)
        frame.pack(fill="both", expand=True)

        ttk.Label(frame, text="DEG Pipeline & Visualizer",
                  font=("", 14, "bold")).pack(pady=(0, 4))
        ttk.Label(frame, text=f"Version {APP_VERSION}",
                  font=("", 9)).pack()
        ttk.Label(frame, text="",
                  font=("", 6)).pack()

        ttk.Label(frame, text="If you use this tool in your research,\nplease cite it and visit:",
                  justify="center", font=("", 10)).pack(pady=6)

        btn_frame = ttk.Frame(frame)
        btn_frame.pack(pady=8)

        github_btn = ttk.Button(btn_frame, text="🌐 GitHub Repository",
                                command=lambda: webbrowser.open(APP_GITHUB))
        github_btn.pack(side="left", padx=10)

        linkedin_btn = ttk.Button(btn_frame, text="💼 LinkedIn Profile",
                                  command=lambda: webbrowser.open(APP_LINKEDIN))
        linkedin_btn.pack(side="left", padx=10)

        ttk.Label(frame, text="",
                  font=("", 4)).pack()
        ttk.Button(frame, text="Close", command=self.destroy).pack(pady=6)

        self.protocol("WM_DELETE_WINDOW", self.destroy)


# ==================================================================
# 8) Main application
# ==================================================================
class DEGPipelineApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(f"DEG Pipeline & Visualizer — Multi-Cancer Analysis  (v{APP_VERSION})")

        icon_path = get_icon_path()
        if icon_path:
            try:
                if icon_path.suffix.lower() == ".ico":
                    self.iconbitmap(str(icon_path))
                if icon_path.suffix.lower() == ".png":
                    img = tk.PhotoImage(file=str(icon_path))
                    self.iconphoto(True, img)
                    self._icon_image = img
            except Exception as e:
                print(f"Could not set icon: {e}")

        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        w, h = min(1360, sw - 80), min(900, sh - 80)
        self.geometry(f"{w}x{h}+{(sw - w)//2}+{(sh - h)//2}")
        self.minsize(900, 560)

        self.log_queue = queue.Queue()
        self.after(150, self._poll_log_queue)

        self._build_menu()

        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True)

        tab1_wrap = ScrollableFrame(nb)
        tab2_wrap = ScrollableFrame(nb)
        nb.add(tab1_wrap, text="  Step 1 — Raw Data -> DEG Analysis  ")
        nb.add(tab2_wrap, text="  Step 2 — Q1 Figures  ")

        self._build_step1(tab1_wrap.inner)
        self._build_step2(tab2_wrap.inner)

        self.after(400, self._show_welcome)

    def _show_welcome(self):
        WelcomeDialog(self)

    def _build_menu(self):
        menubar = tk.Menu(self)
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About / Citation", command=self._show_about)
        help_menu.add_command(label="GitHub Repository",
                              command=lambda: webbrowser.open(APP_GITHUB))
        help_menu.add_command(label="LinkedIn Profile",
                              command=lambda: webbrowser.open(APP_LINKEDIN))
        menubar.add_cascade(label="Help", menu=help_menu)
        self.config(menu=menubar)

    def _show_about(self):
        messagebox.showinfo(
            "About — DEG Pipeline & Visualizer",
            f"DEG Pipeline & Visualizer  v{APP_VERSION}\n"
            f"Author: {APP_AUTHOR}\n\n"
            "If you use this tool in your research, please cite it and\n"
            "link back to the repository below.\n\n"
            f"GitHub:   {APP_GITHUB}\n"
            f"LinkedIn: {APP_LINKEDIN}"
        )

    # ============================================================
    # STEP 1 UI
    # ============================================================
    def _build_step1(self, parent):
        parent.configure(padding=16)

        intro = ttk.Label(parent, wraplength=980, justify="left", font=("", 10),
                           text="Point this at your raw RNA-seq download (Gene_Expression_Quantification) "
                                "folder with per-sample subfolders, each containing one *.tsv file, and "
                                "your GDC sample sheet. It reads every file, runs PyDESeq2 (Tumor vs Normal), "
                                "and saves result files into Data / Pictures / Logs subfolders.\n"
                                "This step does NOT create plots. Go to Step 2, load the "
                                "analysis_cache.csv, and request figures.")
        intro.pack(fill="x", pady=(0, 12))

        box = ttk.LabelFrame(parent, text="Inputs", padding=12)
        box.pack(fill="x", pady=(0, 12))
        box.columnconfigure(1, weight=1)

        self.s1_counts_dir = tk.StringVar()
        self.s1_sample_sheet = tk.StringVar()
        self.s1_output_dir = tk.StringVar()
        self.s1_label = tk.StringVar(value="TCGA-Project")
        self.s1_min_count = tk.IntVar(value=10)
        self.s1_log2fc = tk.DoubleVar(value=1.0)
        self.s1_padj = tk.DoubleVar(value=0.05)

        r = 0
        ttk.Label(box, text="Gene_Expression_Quantification folder\n(contains per-sample subfolders)",
                  justify="left").grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(box, textvariable=self.s1_counts_dir).grid(row=r, column=1, sticky="ew", padx=6)
        ttk.Button(box, text="Browse...",
                    command=lambda: self._browse_dir(self.s1_counts_dir)).grid(row=r, column=2)
        r += 1

        ttk.Label(box, text="GDC sample sheet (.tsv)").grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(box, textvariable=self.s1_sample_sheet).grid(row=r, column=1, sticky="ew", padx=6)
        ttk.Button(box, text="Browse...",
                    command=lambda: self._browse_file(self.s1_sample_sheet,
                                                       [("TSV files", "*.tsv"), ("All files", "*.*")])
                    ).grid(row=r, column=2)
        r += 1

        ttk.Label(box, text="Output folder\n(Data / Pictures / Logs will be created inside it)",
                  justify="left").grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(box, textvariable=self.s1_output_dir).grid(row=r, column=1, sticky="ew", padx=6)
        ttk.Button(box, text="Browse...",
                    command=lambda: self._browse_dir(self.s1_output_dir)).grid(row=r, column=2)
        r += 1

        ttk.Label(box, text="Project label (e.g. TCGA-COAD, TCGA-READ, or any custom name)"
                  ).grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(box, textvariable=self.s1_label).grid(row=r, column=1, sticky="ew", padx=6)
        r += 1

        ttk.Label(box, text="Minimum read count per gene"
                  ).grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(box, textvariable=self.s1_min_count, width=10).grid(row=r, column=1, sticky="w", padx=6)
        r += 1

        ttk.Label(box, text="|log2FC| threshold"
                  ).grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(box, textvariable=self.s1_log2fc, width=10).grid(row=r, column=1, sticky="w", padx=6)
        r += 1

        ttk.Label(box, text="Adjusted P-value (FDR) threshold"
                  ).grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(box, textvariable=self.s1_padj, width=10).grid(row=r, column=1, sticky="w", padx=6)
        r += 1

        out_names = ttk.LabelFrame(parent, text="Output files (Data/)", padding=10)
        out_names.pack(fill="x", pady=(0, 12))
        for key in ["all", "cache", "sig", "up", "down", "norm", "groups", "gene_list"]:
            note = "  <- load THIS into Step 2" if key == "cache" else ""
            note = "  <- publication gene list" if key == "gene_list" else note
            ttk.Label(out_names, text=f"Data/<label>_{STD_NAMES[key]}{note}",
                      font=("Consolas", 9)).pack(anchor="w")
        ttk.Label(out_names, text=f"Data/<label>_{STD_NAMES['log']}",
                  font=("Consolas", 9), foreground="#555555").pack(anchor="w")

        btn_row = ttk.Frame(parent)
        btn_row.pack(fill="x", pady=(0, 10))
        self.s1_run_btn = ttk.Button(btn_row, text="Run Analysis", command=self._start_step1)
        self.s1_run_btn.pack(side="left")
        ttk.Button(btn_row, text="📋 Copy Log", command=self._copy_s1_log).pack(side="left", padx=10)

        # Progress bar
        progress_frame = ttk.LabelFrame(parent, text="Progress", padding=6)
        progress_frame.pack(fill="x", pady=(0, 6))
        self.progress_var = tk.DoubleVar(value=0)
        self.progress_bar = ttk.Progressbar(progress_frame, variable=self.progress_var,
                                             maximum=100, mode="determinate")
        self.progress_bar.pack(fill="x", pady=(2, 0))
        self.progress_label = ttk.Label(progress_frame, text="Ready", font=("", 9))
        self.progress_label.pack(anchor="w", pady=(2, 0))

        # Style the progress bar to be green
        style = ttk.Style()
        style.theme_use("clam")
        style.configure("green.Horizontal.TProgressbar",
                        troughcolor="#e0e0e0",
                        background="#4CAF50",
                        lightcolor="#4CAF50",
                        darkcolor="#388E3C")
        self.progress_bar.configure(style="green.Horizontal.TProgressbar")

        log_box = ttk.LabelFrame(parent, text="Progress Log", padding=6)
        log_box.pack(fill="both", expand=True)
        self.s1_log = tk.Text(log_box, height=22, wrap="word", font=("Consolas", 9))
        self.s1_log.pack(fill="both", expand=True)

    def _update_progress(self, pct, msg=""):
        self.log_queue.put(("progress", (pct, msg)))

    def _browse_dir(self, var):
        p = filedialog.askdirectory()
        if p:
            var.set(p)

    def _browse_file(self, var, filetypes):
        p = filedialog.askopenfilename(filetypes=filetypes)
        if p:
            var.set(p)

    def _copy_s1_log(self):
        self.clipboard_clear()
        self.clipboard_append(self.s1_log.get("1.0", "end-1c"))
        messagebox.showinfo("Copied", "Log content copied to clipboard.")

    def _start_step1(self):
        if not self.s1_counts_dir.get() or not self.s1_sample_sheet.get() or not self.s1_output_dir.get():
            messagebox.showwarning("Missing input", "Please fill in all fields.")
            return
        self.s1_run_btn.configure(state="disabled")
        self.s1_log.delete("1.0", "end")
        self.progress_var.set(0)
        self.progress_label.configure(text="Starting...")
        t = threading.Thread(target=self._run_step1_worker, daemon=True)
        t.start()

    def _run_step1_worker(self):
        log_lines = []

        def qlog(msg):
            log_lines.append(str(msg))
            self.log_queue.put(("s1", msg))

        def progress_cb(pct, msg=""):
            self.log_queue.put(("progress", (pct, msg)))

        try:
            paths = run_deg_pipeline(
                counts_root=self.s1_counts_dir.get(),
                sample_sheet_path=self.s1_sample_sheet.get(),
                output_dir=self.s1_output_dir.get(),
                label=self.s1_label.get().strip().replace(" ", "_"),
                min_count=self.s1_min_count.get(),
                log2fc_thr=self.s1_log2fc.get(),
                padj_thr=self.s1_padj.get(),
                log=qlog,
                progress_callback=progress_cb,
            )
            try:
                dirs = get_output_dirs(Path(self.s1_output_dir.get()))
                label = self.s1_label.get().strip().replace(" ", "_")
                (dirs["logs"] / f"{label}_{STD_NAMES['log']}").write_text(
                    "\n".join(log_lines), encoding="utf-8")
            except Exception:
                pass
            self.log_queue.put(("s1_done", paths))
        except Exception as e:
            qlog("[FATAL] " + str(e))
            qlog(traceback.format_exc())
            self.log_queue.put(("s1_failed", None))

    def _poll_log_queue(self):
        try:
            for _ in range(50):
                kind, payload = self.log_queue.get_nowait()
                if kind == "s1":
                    self.s1_log.insert("end", str(payload) + "\n")
                    self.s1_log.see("end")
                elif kind == "progress":
                    pct, msg = payload
                    self.progress_var.set(pct)
                    self.progress_label.configure(text=f"{msg} ({pct:.0f}%)")
                elif kind == "s1_done":
                    self.s1_run_btn.configure(state="normal")
                    self.progress_var.set(100)
                    self.progress_label.configure(text="Complete! (100%)")
                    self._on_step1_done(payload)
                elif kind == "s1_failed":
                    self.s1_run_btn.configure(state="normal")
                    self.progress_label.configure(text="Failed!")
                    messagebox.showerror("Analysis failed", "See the progress log for details.")
        except queue.Empty:
            pass
        self.after(150, self._poll_log_queue)

    def _on_step1_done(self, paths):
        self.s1_log.insert("end", "\nDONE. Result files saved. No plots generated.\n")
        if paths:
            self.s2_files["all"].set(str(paths.get("cache", "")))
            self.s2_files["sig"].set(str(paths.get("sig", "")))
            self.s2_files["norm"].set(str(paths.get("norm", "")))
            self.s2_files["groups"].set(str(paths.get("groups", "")))
            self.s2_output_dir.set(str(Path(paths.get("all", "")).parent.parent))
            self.s2_label.set(self.s1_label.get())
            messagebox.showinfo(
                "Analysis complete",
                "DEG analysis finished.\n\n"
                "Go to 'Step 2 — Q1 Figures' tab. Cache file is loaded. "
                "List genes to label, then click 'Generate All Standard Plots'.")

    # ============================================================
    # STEP 2 UI
    # ============================================================
    def _build_step2(self, parent):
        parent.configure(padding=16)
        outer = ttk.Panedwindow(parent, orient="horizontal")
        outer.pack(fill="both", expand=True)
        left = ttk.Frame(outer, padding=6)
        right = ttk.Frame(outer, padding=6)
        outer.add(left, weight=0)
        outer.add(right, weight=1)

        self.s2_files = {
            "all": tk.StringVar(), "sig": tk.StringVar(),
            "norm": tk.StringVar(), "groups": tk.StringVar(),
        }
        self.s2_output_dir = tk.StringVar()
        self.s2_label = tk.StringVar(value="TCGA-Project")
        self.plot_type = tk.StringVar(value="Volcano Plot")
        self.log2fc_thr = tk.DoubleVar(value=1.0)
        self.padj_thr = tk.DoubleVar(value=0.05)
        self.top_n = tk.IntVar(value=20)
        self.palette_choice = tk.StringVar(value="Reference (red/blue/grey)")
        self.fig_w = tk.DoubleVar(value=9.0)
        self.fig_h = tk.DoubleVar(value=7.0)
        self.dpi = tk.IntVar(value=600)
        self.fmt_png = tk.BooleanVar(value=True)
        self.fmt_pdf = tk.BooleanVar(value=True)
        self.fmt_tiff = tk.BooleanVar(value=False)
        self.fmt_svg = tk.BooleanVar(value=False)
        self.export_processed = tk.BooleanVar(value=True)
        self.current_fig = None
        self.current_deg_df = None

        self._s2_file_section(left)
        self._s2_plot_selector(left)
        self._s2_params_section(left)
        self._s2_output_section(left)
        self._s2_buttons(left)
        self._s2_log(left)
        self._s2_preview(right)
        self._on_plot_type_change()

        outer.sashpos(0, 550)

    def _s2_file_row(self, parent, label, var, row, fname_hint, required=False):
        req = " (required)" if required else " (optional)"
        ttk.Label(parent, text=f"{label}{req}\n({fname_hint})", justify="left",
                  font=("", 9)).grid(row=row, column=0, sticky="w", pady=3)
        ttk.Entry(parent, textvariable=var, width=24).grid(row=row, column=1, sticky="ew", padx=4)
        ttk.Button(parent, text="Browse...",
                    command=lambda: self._browse_file(var, [("CSV", "*.csv"), ("All", "*.*")])
                    ).grid(row=row, column=2)

    def _s2_file_section(self, parent):
        box = ttk.LabelFrame(parent, text="1) Load data — no plots until you ask", padding=10)
        box.pack(fill="x", pady=(0, 10))
        box.columnconfigure(1, weight=1)

        ttk.Label(box, text="Analysis results / cache file", font=("", 10, "bold")
                  ).grid(row=0, column=0, sticky="w", pady=(0, 2))
        self._s2_file_row(box, "File", self.s2_files["all"], 1,
                           f"<label>_{STD_NAMES['cache']} or <label>_{STD_NAMES['all']}", required=True)
        ttk.Label(box, text="This file drives Volcano, MA, and Bar Chart.",
                  font=("", 8, "italic"), foreground="#555555", wraplength=380
                  ).grid(row=2, column=0, columnspan=3, sticky="w", pady=(0, 8))

        ttk.Separator(box).grid(row=3, column=0, columnspan=3, sticky="ew", pady=6)
        ttk.Label(box, text="Optional — for expression Heatmap:", font=("", 9, "italic")
                  ).grid(row=4, column=0, columnspan=3, sticky="w")
        self._s2_file_row(box, "Significant DEGs", self.s2_files["sig"], 5, STD_NAMES["sig"])
        self._s2_file_row(box, "Normalized counts", self.s2_files["norm"], 6, STD_NAMES["norm"])
        self._s2_file_row(box, "Sample grouping", self.s2_files["groups"], 7, STD_NAMES["groups"])

        btn_frame = ttk.Frame(box)
        btn_frame.grid(row=8, column=0, columnspan=3, pady=10)
        ttk.Button(btn_frame, text="📂 Load from analysis_log.txt",
                   command=self._load_from_log).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="📁 Load from Data Folder",
                   command=self._load_from_data_folder).pack(side="left", padx=5)

    def _load_from_log(self):
        log_path = filedialog.askopenfilename(
            title="Select analysis_log.txt",
            filetypes=[("Log files", "*.txt"), ("All files", "*.*")]
        )
        if not log_path:
            return
        try:
            with open(log_path, "r", encoding="utf-8") as f:
                content = f.read()
            match_data = re.search(r"All result files saved to:\s*(.+)", content)
            if not match_data:
                raise ValueError("Could not find 'All result files saved to:' line in log.")
            data_dir = Path(match_data.group(1).strip())
            output_dir = data_dir.parent
            self.s2_output_dir.set(str(output_dir))
            cache_match = re.search(r"\[cache\]\s+(\S+)_analysis_cache\.csv", content)
            if cache_match:
                label = cache_match.group(1)
                self.s2_label.set(label)
                self.s2_files["all"].set(str(data_dir / f"{label}_{STD_NAMES['cache']}"))
                self.s2_files["sig"].set(str(data_dir / f"{label}_{STD_NAMES['sig']}"))
                self.s2_files["norm"].set(str(data_dir / f"{label}_{STD_NAMES['norm']}"))
                self.s2_files["groups"].set(str(data_dir / f"{label}_{STD_NAMES['groups']}"))
            else:
                self.s2_files["all"].set(str(data_dir / STD_NAMES["cache"]))
            self.log2("Loaded analysis from log file.")
        except Exception as e:
            messagebox.showerror("Load failed", str(e))

    def _load_from_data_folder(self):
        folder = filedialog.askdirectory(title="Select Data folder")
        if not folder:
            return
        data_path = Path(folder)
        if not data_path.exists():
            messagebox.showerror("Error", "Folder does not exist.")
            return
        cache_candidates = list(data_path.glob(f"*{STD_NAMES['cache']}"))
        if not cache_candidates:
            cache_candidates = list(data_path.glob(f"*{STD_NAMES['all']}"))
        if not cache_candidates:
            messagebox.showerror("Error", f"No file matching *{STD_NAMES['cache']} or *{STD_NAMES['all']} found.")
            return
        if len(cache_candidates) > 1:
            chosen = filedialog.askopenfilename(initialdir=folder, title="Select cache file",
                                                filetypes=[("CSV", "*.csv")])
            if not chosen:
                return
            cache_file = Path(chosen)
        else:
            cache_file = cache_candidates[0]
        stem = cache_file.stem
        if stem.endswith("_analysis_cache"):
            label = stem[:-len("_analysis_cache")]
        elif stem.endswith("_All_DEGs_results"):
            label = stem[:-len("_All_DEGs_results")]
        else:
            label = "TCGA-Project"
        self.s2_label.set(label)
        self.s2_files["all"].set(str(cache_file))

        for key, fname in [("sig", STD_NAMES["sig"]), ("norm", STD_NAMES["norm"]), ("groups", STD_NAMES["groups"])]:
            candidates = list(data_path.glob(f"*{fname}"))
            if candidates:
                self.s2_files[key].set(str(max(candidates, key=lambda p: p.stat().st_mtime)))
            else:
                self.s2_files[key].set("")
        self.s2_output_dir.set(str(data_path.parent))
        self.log2(f"Loaded data from folder: {folder}")

    def _s2_plot_selector(self, parent):
        box = ttk.LabelFrame(parent, text="2) Plot type", padding=10)
        box.pack(fill="x", pady=(0, 10))
        combo = ttk.Combobox(box, textvariable=self.plot_type, state="readonly",
                              values=list(PLOT_DESCRIPTIONS.keys()))
        combo.pack(fill="x")
        combo.bind("<<ComboboxSelected>>", lambda e: self._on_plot_type_change())
        self.desc_label = tk.Text(box, height=5, wrap="word", bg=self.cget("bg"),
                                   relief="flat", font=("", 9))
        self.desc_label.pack(fill="x", pady=(6, 0))
        self.desc_label.configure(state="disabled")

    def _on_plot_type_change(self):
        self.desc_label.configure(state="normal")
        self.desc_label.delete("1.0", "end")
        self.desc_label.insert("1.0", PLOT_DESCRIPTIONS[self.plot_type.get()])
        self.desc_label.configure(state="disabled")

    def _s2_params_section(self, parent):
        box = ttk.LabelFrame(parent, text="3) Parameters", padding=10)
        box.pack(fill="x", pady=(0, 10))
        box.columnconfigure(1, weight=1)
        r = 0
        ttk.Label(box, text="Project label").grid(row=r, column=0, sticky="w")
        ttk.Entry(box, textvariable=self.s2_label).grid(row=r, column=1, sticky="ew", padx=4)
        r += 1
        ttk.Label(box, text="|log2FC| threshold").grid(row=r, column=0, sticky="w")
        ttk.Entry(box, textvariable=self.log2fc_thr, width=10).grid(row=r, column=1, sticky="w", padx=4)
        r += 1
        ttk.Label(box, text="Adjusted P-value (FDR) threshold").grid(row=r, column=0, sticky="w")
        ttk.Entry(box, textvariable=self.padj_thr, width=10).grid(row=r, column=1, sticky="w", padx=4)
        r += 1
        ttk.Label(box, text="Top N genes (heatmap only)").grid(row=r, column=0, sticky="w")
        ttk.Entry(box, textvariable=self.top_n, width=10).grid(row=r, column=1, sticky="w", padx=4)
        r += 1

        ttk.Label(box, text="Genes to label on plots", font=("", 9, "bold")
                  ).grid(row=r, column=0, columnspan=2, sticky="w", pady=(8, 0))
        r += 1
        ttk.Label(box, text="One per line, or comma/semicolon separated. Leave empty for no labels.",
                  font=("", 8, "italic"), foreground="#555555"
                  ).grid(row=r, column=0, columnspan=2, sticky="w")
        r += 1
        self.gene_text = tk.Text(box, height=5, width=30, font=("Consolas", 9))
        self.gene_text.grid(row=r, column=0, columnspan=2, sticky="ew", pady=(2, 6))
        r += 1

        ttk.Label(box, text="Color palette").grid(row=r, column=0, sticky="w")
        ttk.Combobox(box, textvariable=self.palette_choice, state="readonly",
                      values=["Reference (red/blue/grey)", "Deep (dark red/navy/grey)"]
                      ).grid(row=r, column=1, sticky="ew", padx=4)
        r += 1
        ttk.Label(box, text="Figure size (W x H inches)").grid(row=r, column=0, sticky="w")
        sizef = ttk.Frame(box)
        sizef.grid(row=r, column=1, sticky="w", padx=4)
        ttk.Entry(sizef, textvariable=self.fig_w, width=6).pack(side="left")
        ttk.Label(sizef, text=" x ").pack(side="left")
        ttk.Entry(sizef, textvariable=self.fig_h, width=6).pack(side="left")
        r += 1
        ttk.Label(box, text="Export DPI").grid(row=r, column=0, sticky="w")
        ttk.Entry(box, textvariable=self.dpi, width=10).grid(row=r, column=1, sticky="w", padx=4)
        r += 1
        ttk.Checkbutton(box, text="Also export processed data CSV + summary TXT",
                         variable=self.export_processed).grid(row=r, column=0, columnspan=2,
                                                              sticky="w", pady=(4, 0))

    def _s2_output_section(self, parent):
        box = ttk.LabelFrame(parent, text="4) Output", padding=10)
        box.pack(fill="x", pady=(0, 10))
        box.columnconfigure(1, weight=1)
        ttk.Label(box, text="Output folder\n(Pictures/ and Data/ created inside)",
                  justify="left").grid(row=0, column=0, sticky="w")
        ttk.Entry(box, textvariable=self.s2_output_dir).grid(row=0, column=1, sticky="ew", padx=4)
        ttk.Button(box, text="Browse...",
                    command=lambda: self._browse_dir(self.s2_output_dir)).grid(row=0, column=2)
        fmt_row = ttk.Frame(box)
        fmt_row.grid(row=1, column=0, columnspan=3, sticky="w", pady=(6, 0))
        ttk.Checkbutton(fmt_row, text="PNG", variable=self.fmt_png).pack(side="left", padx=4)
        ttk.Checkbutton(fmt_row, text="PDF (vector)", variable=self.fmt_pdf).pack(side="left", padx=4)
        ttk.Checkbutton(fmt_row, text="TIFF", variable=self.fmt_tiff).pack(side="left", padx=4)
        ttk.Checkbutton(fmt_row, text="SVG", variable=self.fmt_svg).pack(side="left", padx=4)

    def _s2_buttons(self, parent):
        top_row = ttk.Frame(parent)
        top_row.pack(fill="x", pady=(0, 4))
        ttk.Button(top_row, text="Generate All Standard Plots (Volcano + MA + Bar Chart)",
                    command=self.generate_all_standard).pack(fill="x")

        row = ttk.Frame(parent)
        row.pack(fill="x", pady=(4, 10))
        ttk.Button(row, text="Preview selected plot", command=lambda: self.generate_preview()
                    ).pack(side="left", expand=True, fill="x", padx=2)
        ttk.Button(row, text="Save selected plot", command=self.save_figure
                    ).pack(side="left", expand=True, fill="x", padx=2)

    def _s2_log(self, parent):
        box = ttk.LabelFrame(parent, text="Log", padding=6)
        box.pack(fill="both", expand=True)
        self.s2_log_text = tk.Text(box, height=10, wrap="word", font=("Consolas", 9))
        self.s2_log_text.pack(fill="both", expand=True)

    def _s2_preview(self, parent):
        self.fig = plt.Figure(figsize=(9, 7), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def log2(self, msg):
        self.s2_log_text.insert("end", msg + "\n")
        self.s2_log_text.see("end")
        self.update_idletasks()

    def _get_palette(self):
        return PALETTE if self.palette_choice.get().startswith("Reference") else PALETTE_ALT

    def _get_label_gene_list(self):
        raw = self.gene_text.get("1.0", "end")
        return parse_gene_list(raw)

    def _validate_numeric(self, value, min_val=None, max_val=None, name="value"):
        v = value.get()
        try:
            v = float(v)
        except (ValueError, tk.TclError):
            raise ValueError(f"Invalid {name}: must be numeric.")
        if min_val is not None and v < min_val:
            raise ValueError(f"{name} must be >= {min_val}")
        if max_val is not None and v > max_val:
            raise ValueError(f"{name} must be <= {max_val}")
        return v

    def generate_preview(self, labeled=False):
        try:
            fig_w = self._validate_numeric(self.fig_w, 1, 50, "Figure width")
            fig_h = self._validate_numeric(self.fig_h, 1, 50, "Figure height")
            log2fc_thr = self._validate_numeric(self.log2fc_thr, 0, 10, "log2FC threshold")
            padj_thr = self._validate_numeric(self.padj_thr, 1e-10, 1, "padj threshold")
            top_n = int(self._validate_numeric(self.top_n, 1, 500, "Top N"))

            self.fig.clf()
            ax = self.fig.add_subplot(111)
            plot_type = self.plot_type.get()
            palette = self._get_palette()
            label = self.s2_label.get().strip() or "Sample"
            main_path = self.s2_files["all"].get().strip()

            if plot_type == "Volcano Plot":
                if not main_path:
                    raise ValueError("Select the Analysis results / cache file.")
                df = load_csv(main_path)
                deg = prep_deg_table(df, log2fc_thr, padj_thr)
                requested = self._get_label_gene_list() if labeled else []
                genes = []
                if labeled and requested:
                    genes, not_found = find_matching_genes(deg, requested)
                    if not_found:
                        self.log2(f"[WARN] Genes not found: {', '.join(not_found)}")
                plot_volcano(ax, deg, log2fc_thr, padj_thr, palette,
                             f"Volcano Plot - {label} (Tumor vs Normal)",
                             label_genes=genes, log_func=self.log2)
                self.current_deg_df = deg

            elif plot_type == "MA Plot":
                if not main_path:
                    raise ValueError("Select the Analysis results / cache file.")
                df = load_csv(main_path)
                requested = self._get_label_gene_list() if labeled else []
                genes = []
                if labeled and requested:
                    gene_col = pick_col(df, ["gene_name", "gene_id"])
                    if gene_col:
                        gene_series = df[gene_col].astype(str)
                        found, not_found = find_matching_genes(pd.DataFrame({"gene": gene_series}), requested)
                        if not_found:
                            self.log2(f"[WARN] Genes not found in MA: {', '.join(not_found)}")
                        if found:
                            genes = found
                plot_ma(ax, df, log2fc_thr, padj_thr, palette,
                        f"MA Plot - {label}", label_genes=genes, log_func=self.log2)
                self.current_deg_df = prep_deg_table(df, log2fc_thr, padj_thr)

            elif plot_type == "Summary Bar Chart":
                if not main_path:
                    raise ValueError("Select the Analysis results / cache file.")
                df = load_csv(main_path)
                deg = prep_deg_table(df, log2fc_thr, padj_thr)
                up_n = int((deg["class"] == "Up").sum())
                down_n = int((deg["class"] == "Down").sum())
                plot_summary_bar(ax, up_n, down_n, f"Differentially Expressed Genes - {label}")
                self.current_deg_df = deg

            elif plot_type == "Heatmap (Top DEGs)":
                sig_path = self.s2_files["sig"].get().strip() or main_path
                if not sig_path:
                    raise ValueError("Select a Significant DEGs file.")
                sig_df = load_csv(sig_path)
                deg = prep_deg_table(sig_df, log2fc_thr, padj_thr)
                self.current_deg_df = deg
                norm_path = self.s2_files["norm"].get().strip()
                group_path = self.s2_files["groups"].get().strip()
                if norm_path and group_path:
                    norm_counts = pd.read_csv(norm_path, index_col=0)
                    sample_groups = pd.read_csv(group_path, index_col=0)
                    id_col = pick_col(sig_df, ["gene_id"])
                    name_col = pick_col(sig_df, ["gene_name", "gene_id"])
                    top = sig_df.copy()
                    padj_col = pick_col(top, ["padj", "FDR", "pvalue"])
                    top = top.sort_values(padj_col).head(top_n)
                    gene_ids = top[id_col].values if id_col else top[name_col].values
                    gene_labels = dict(zip(top[id_col], top[name_col])) if id_col and name_col else {}
                    plot_heatmap_expression(ax, norm_counts, sample_groups, gene_ids,
                                            gene_labels, f"Top {top_n} DEGs - {label}")
                else:
                    plot_heatmap_rank(ax, deg, top_n, f"Top {top_n} DEGs - {label}")

            self.fig.set_size_inches(fig_w, fig_h)
            self.fig.tight_layout()
            self.canvas.draw_idle()
            self.current_fig = self.fig

        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.log2(f"[ERROR] {e}")

    def _export_current(self, pictures_dir, base_name):
        formats = []
        if self.fmt_png.get():
            formats.append("png")
        if self.fmt_pdf.get():
            formats.append("pdf")
        if self.fmt_tiff.get():
            formats.append("tiff")
        if self.fmt_svg.get():
            formats.append("svg")
        if not formats:
            raise ValueError("Check at least one export format.")
        saved = []
        for fmt in formats:
            out_path = Path(pictures_dir) / f"{base_name}.{fmt}"
            self.current_fig.savefig(out_path, dpi=self.dpi.get(), bbox_inches="tight",
                                      format=fmt if fmt != "tiff" else "tiff")
            saved.append(str(out_path))
        return saved

    def _do_save_current(self):
        outdir = self.s2_output_dir.get().strip()
        if not outdir:
            raise ValueError("Choose an output folder first.")
        dirs = get_output_dirs(Path(outdir))
        label = self.s2_label.get().strip().replace(" ", "_") or "sample"
        plot_key = self.plot_type.get().lower().replace(" ", "_").replace("(", "").replace(")", "")

        self.generate_preview(labeled=False)
        if self.current_fig is None:
            raise RuntimeError("Could not generate the plot.")

        base_name = f"{label}_{plot_key}"
        all_saved = self._export_current(dirs["pictures"], base_name)

        if self.export_processed.get() and self.current_deg_df is not None:
            dp, sp = save_processed_data_and_summary(self.current_deg_df, dirs["data"], base_name,
                                                       self.log2fc_thr.get(), self.padj_thr.get())
            all_saved += [str(dp), str(sp)]

        gene_list = self._get_label_gene_list()
        if gene_list and self.plot_type.get() in ("Volcano Plot", "MA Plot"):
            self.generate_preview(labeled=True)
            all_saved += self._export_current(dirs["pictures"], f"{label}_{plot_key}_labeled")

        return all_saved

    def save_figure(self):
        try:
            saved = self._do_save_current()
            self.log2("Saved:\n  " + "\n  ".join(saved))
            messagebox.showinfo("Saved", f"Figure(s) saved to:\n{self.s2_output_dir.get()}")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.log2(f"[ERROR] {e}")

    def generate_all_standard(self):
        if not self.s2_files["all"].get().strip():
            messagebox.showwarning("Missing file", "Select the cache file first.")
            return
        if not self.s2_output_dir.get().strip():
            messagebox.showwarning("Missing output", "Choose an output folder first.")
            return

        main_path = self.s2_files["all"].get().strip()
        label = self.s2_label.get().strip() or "Sample"
        palette = self._get_palette()
        log2fc_thr = self.log2fc_thr.get()
        padj_thr = self.padj_thr.get()

        try:
            df = load_csv(main_path)
            deg = prep_deg_table(df, log2fc_thr, padj_thr)
        except Exception as e:
            messagebox.showerror("Data error", str(e))
            return

        outdir = self.s2_output_dir.get().strip()
        dirs = get_output_dirs(Path(outdir))
        pictures_dir = dirs["pictures"]

        summary_lines = []
        for pt in ["Volcano Plot", "MA Plot", "Summary Bar Chart"]:
            self.plot_type.set(pt)
            try:
                saved = self._do_save_current()
                line = f"{pt}: {len(saved)} file(s) saved"
            except Exception as e:
                line = f"{pt}: skipped — {e}"
            summary_lines.append(line)
            self.log2(line)

        plt.close(self.fig)
        self.fig = plt.Figure(figsize=(18, 6))
        axes = self.fig.subplots(1, 3)
        plot_volcano(axes[0], deg, log2fc_thr, padj_thr, palette, f"Volcano Plot - {label}")
        plot_ma(axes[1], df, log2fc_thr, padj_thr, palette, f"MA Plot - {label}")
        up_n = int((deg["class"] == "Up").sum())
        down_n = int((deg["class"] == "Down").sum())
        plot_summary_bar(axes[2], up_n, down_n, f"DEGs - {label}")
        self.fig.tight_layout()
        self.canvas.figure = self.fig
        self.canvas.draw()
        self.current_fig = self.fig

        base_name = f"{label}_combined_standard"
        try:
            combined_saved = self._export_current(pictures_dir, base_name)
            summary_lines.append(f"Combined figure: {len(combined_saved)} file(s) saved")
        except Exception as e:
            summary_lines.append(f"Combined figure: skipped — {e}")

        requested = self._get_label_gene_list()
        if requested:
            volcano_genes, v_not = find_matching_genes(deg, requested)
            if v_not:
                self.log2(f"[WARN] Genes not found (Volcano): {', '.join(v_not)}")
            gene_col = pick_col(df, ["gene_name", "gene_id"])
            ma_genes = []
            if gene_col:
                gene_series = df[gene_col].astype(str)
                ma_genes, m_not = find_matching_genes(pd.DataFrame({"gene": gene_series}), requested)
                if m_not:
                    self.log2(f"[WARN] Genes not found (MA): {', '.join(m_not)}")

            if volcano_genes or ma_genes:
                fig_labeled, axes_lbl = plt.subplots(1, 3, figsize=(18, 6))
                plot_volcano(axes_lbl[0], deg, log2fc_thr, padj_thr, palette,
                             f"Volcano Plot - {label}", label_genes=volcano_genes, log_func=self.log2)
                plot_ma(axes_lbl[1], df, log2fc_thr, padj_thr, palette,
                        f"MA Plot - {label}", label_genes=ma_genes, log_func=self.log2)
                plot_summary_bar(axes_lbl[2], up_n, down_n, f"DEGs - {label}")
                fig_labeled.tight_layout()
                old_fig = self.current_fig
                self.current_fig = fig_labeled
                try:
                    labeled_saved = self._export_current(pictures_dir, f"{label}_combined_standard_labeled")
                    summary_lines.append(f"Combined labeled figure: {len(labeled_saved)} file(s) saved")
                except Exception as e:
                    summary_lines.append(f"Combined labeled figure: skipped — {e}")
                self.current_fig = old_fig
                plt.close(fig_labeled)

        messagebox.showinfo("Batch generation finished", "\n".join(summary_lines))


# ==================================================================
# 9) Entry point
# ==================================================================
if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = DEGPipelineApp()
    app.mainloop()