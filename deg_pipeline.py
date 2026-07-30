"""
================================================================================
  DEG Pipeline & Visualizer  -  Multi-Cancer Analysis Suite
================================================================================
  Author  : Alireza Balaei
  GitHub  : https://github.com/alirezabk1382927-sys
  LinkedIn: https://ir.linkedin.com/in/alireza-balaei-kahnamoei-aa8216344
  Version : 3.1.0
================================================================================
"""

# ---------------------------------------------------------------------------
# 0)  Frozen-executable detection  &  BLAS / loky thread caps
# ---------------------------------------------------------------------------
import os
import sys
import multiprocessing

IS_FROZEN = getattr(sys, "frozen", False)
SAFE_N_CPUS = 1 if IS_FROZEN else max(1, min(4, (os.cpu_count() or 1)))

os.environ["LOKY_MAX_CPU_COUNT"] = str(SAFE_N_CPUS)
os.environ.setdefault("OMP_NUM_THREADS", str(SAFE_N_CPUS))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(SAFE_N_CPUS))
os.environ.setdefault("MKL_NUM_THREADS", str(SAFE_N_CPUS))
os.environ.setdefault("NUMBA_NUM_THREADS", str(SAFE_N_CPUS))

# ---------------------------------------------------------------------------
# 1)  Dependency check
# ---------------------------------------------------------------------------
import subprocess
import importlib

REQUIRED_PACKAGES = [
    ("pandas",        "pandas"),
    ("numpy",         "numpy"),
    ("matplotlib",    "matplotlib"),
    ("seaborn",       "seaborn"),
    ("adjustText",    "adjustText"),
    ("pydeseq2",      "pydeseq2"),
    ("customtkinter", "customtkinter>=5.2.0"),
    ("PIL",           "Pillow"),
]


def _pip_install(pip_spec):
    cmd = [sys.executable, "-m", "pip", "install", "--quiet", pip_spec]
    kwargs = {}
    if os.name == "nt":
        kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
    try:
        subprocess.check_call(cmd, timeout=600, **kwargs)
        return True
    except Exception:
        return False


def _ensure_package(import_name, pip_spec):
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


_missing = [imp for imp, spec in REQUIRED_PACKAGES if not _ensure_package(imp, spec)]
if _missing:
    print("Missing packages:", ", ".join(_missing))
    print("Install with: pip install " + " ".join(_missing))
    sys.exit(1)

# ---------------------------------------------------------------------------
# 2)  Imports
# ---------------------------------------------------------------------------
import re
import gzip
import queue
import threading
import traceback
import webbrowser
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass
from typing import Callable, Optional, List, Dict

import numpy as np
import pandas as pd
import matplotlib

try:
    matplotlib.use("TkAgg")
except ImportError:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch
import seaborn as sns

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False

import tkinter as tk
from tkinter import filedialog, messagebox
import customtkinter as ctk
from PIL import Image

# ---------------------------------------------------------------------------
# 3)  Constants
# ---------------------------------------------------------------------------
APP_AUTHOR   = "Alireza Balaei"
APP_GITHUB   = "https://github.com/alirezabk1382927-sys"
APP_LINKEDIN = "https://www.linkedin.com/in/alireza-balaei-kahnamoei-aa8216344/"
APP_ORCID    = "https://orcid.org/0009-0009-9746-6571"
APP_VERSION  = "3.4.0"
APP_NAME     = "DEG Pipeline & Visualizer"


# ---- Scientific Dark theme colors -----------------------------------------
class Theme:
    """Scientific Dark theme  -  high contrast, publication-friendly."""

    # ---- Base ----
    BG_DARK        = "#12151C"      # main app background (very dark navy)
    BG_MID         = "#1A1F2A"      # cards / panels
    BG_LIGHT       = "#232A38"      # elevated surfaces (inputs, textboxes)
    BG_LIGHTER     = "#2C3444"      # hover states, active items

    # ---- Text ----
    TEXT_PRIMARY   = "#E5EBF5"
    TEXT_SECONDARY = "#9CA6B8"
    TEXT_MUTED     = "#6B7488"

    # ---- Borders / Dividers ----
    BORDER         = "#2A3140"
    BORDER_STRONG  = "#3A4356"

    # ---- Accent (Cyan  -  scientific tools feel) ----
    ACCENT         = "#22D3EE"      # bright cyan
    ACCENT_HOVER   = "#0EA5C9"
    ACCENT_SOFT    = "#164E5B"      # dark cyan for backgrounds
    ACCENT_TEXT    = "#67E8F9"

    # ---- Semantic ----
    SUCCESS        = "#22C55E"
    SUCCESS_HOVER  = "#16A34A"
    WARNING        = "#F59E0B"
    DANGER         = "#EF4444"
    DANGER_HOVER   = "#DC2626"
    INFO           = "#3B82F6"

    # ---- Light theme (fallback  -  fresh clean look) ----
    L_BG_DARK      = "#F5F7FA"
    L_BG_MID       = "#FFFFFF"
    L_BG_LIGHT     = "#FFFFFF"
    L_BG_LIGHTER   = "#E4E8EF"
    L_TEXT_PRIMARY = "#0F172A"
    L_TEXT_SECONDARY = "#475569"
    L_TEXT_MUTED   = "#94A3B8"
    L_BORDER       = "#E2E8F0"
    L_BORDER_STRONG= "#CBD5E1"
    L_ACCENT       = "#0891B2"
    L_ACCENT_HOVER = "#0E7490"
    L_ACCENT_SOFT  = "#CFFAFE"

    # ---- Plots ----
    PLOT_UP        = "#E63946"
    PLOT_DOWN      = "#1D4ED8"
    PLOT_NS        = "#B0B4BA"
    PLOT_UP_DEEP   = "#B02A37"
    PLOT_DOWN_DEEP = "#0F2E85"
    PLOT_CUTOFF    = "#16A34A"
    PLOT_GRID      = "#e4e6ea"


# Tuples for dual-theme colors: (light, dark)
def dt(light, dark):
    """Helper to build CTk dual-color tuple."""
    return (light, dark)


COL_BG        = dt(Theme.L_BG_DARK,     Theme.BG_DARK)
COL_CARD      = dt(Theme.L_BG_MID,      Theme.BG_MID)
COL_INPUT     = dt(Theme.L_BG_MID,      Theme.BG_LIGHT)
COL_HOVER     = dt(Theme.L_BG_LIGHTER,  Theme.BG_LIGHTER)
COL_TEXT      = dt(Theme.L_TEXT_PRIMARY,   Theme.TEXT_PRIMARY)
COL_TEXT_SEC  = dt(Theme.L_TEXT_SECONDARY, Theme.TEXT_SECONDARY)
COL_TEXT_MUT  = dt(Theme.L_TEXT_MUTED,     Theme.TEXT_MUTED)
COL_BORDER    = dt(Theme.L_BORDER,         Theme.BORDER)
COL_BORDER_STRONG = dt(Theme.L_BORDER_STRONG, Theme.BORDER_STRONG)
COL_ACCENT    = dt(Theme.L_ACCENT,         Theme.ACCENT)
COL_ACCENT_H  = dt(Theme.L_ACCENT_HOVER,   Theme.ACCENT_HOVER)
COL_ACCENT_S  = dt(Theme.L_ACCENT_SOFT,    Theme.ACCENT_SOFT)


PALETTE     = {"Up": Theme.PLOT_UP,      "Down": Theme.PLOT_DOWN,      "NS": Theme.PLOT_NS}
PALETTE_ALT = {"Up": Theme.PLOT_UP_DEEP, "Down": Theme.PLOT_DOWN_DEEP, "NS": Theme.PLOT_NS}

STD_NAMES = {
    "all":       "All_DEGs_results.csv",
    "cache":     "analysis_cache.csv",
    "sig":       "Significant_DEGs.csv",
    "up":        "Upregulated_in_Tumor.csv",
    "down":      "Downregulated_in_Tumor.csv",
    "norm":      "normalized_counts.csv",
    "groups":    "sample_grouping.csv",
    "dupes":     "duplicate_samples_dropped.csv",
    "log":       "analysis_log.txt",
    "gene_list": "all_genes_complete_list.csv",
}

# Plot type registry: name -> (description, needs_baseMean, needs_heatmap_files)
PLOT_TYPES = {
    "All Plots":         "Generate every plot type at once (recommended). "
                         "Use the arrows in the preview to switch between them.",
    "Volcano Plot":      "Every gene's fold-change (x) vs. significance (y). "
                         "Uses the analysis cache file. Add gene symbols in the "
                         "'Genes to label' box to auto-save a labelled copy.",
    "MA Plot":           "Fold-change vs. average expression level. Requires a "
                         "file with a baseMean column (use the analysis cache, "
                         "not the filtered file).",
    "Summary Bar Chart": "Bar chart summarising the number of significantly "
                         "up- and down-regulated genes at the current thresholds.",
    "Heatmap (Top DEGs)":"Expression heatmap of the top DEGs across samples. "
                         "Best with Significant_DEGs + normalized_counts + "
                         "sample_grouping. Otherwise a fold-change rank heatmap "
                         "is drawn.",
}


# ---------------------------------------------------------------------------
# 4)  Matplotlib style
# ---------------------------------------------------------------------------
def apply_publication_style():
    plt.rcParams.update({
        "font.family":         "sans-serif",
        "font.sans-serif":     ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":           12,
        "axes.linewidth":      0.9,
        "axes.edgecolor":      "#4d4d4d",
        "axes.titleweight":    "bold",
        "xtick.direction":     "out",
        "ytick.direction":     "out",
        "legend.frameon":      False,
        "savefig.dpi":         600,
        "figure.dpi":          105,
        "figure.facecolor":    "white",
        "axes.facecolor":      "white",
    })


apply_publication_style()


def apply_journal_style(ax):
    ax.set_facecolor("white")
    ax.grid(True, color=Theme.PLOT_GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#4d4d4d")
        spine.set_linewidth(0.9)


# ---------------------------------------------------------------------------
# 5)  Icon / folder helpers
# ---------------------------------------------------------------------------
def get_icon_path() -> Optional[Path]:
    """Search for the window/taskbar icon (icon.png / icon.ico)."""
    return _find_asset(("icon.png", "icon.ico"))


def get_main_icon_path() -> Optional[Path]:
    """Search for the big header/welcome logo (main_icon.png).

    Falls back to icon.png if a dedicated main_icon.png is not found."""
    p = _find_asset(("main_icon.png", "main icon.png", "main-icon.png"))
    if p is not None:
        return p
    return get_icon_path()


def _find_asset(names) -> Optional[Path]:
    """Look for any of the given asset filenames in common locations."""
    candidates = []
    if getattr(sys, "frozen", False):
        base = Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
        for n in names:
            candidates.append(base / n)
            candidates.append(Path(sys.executable).parent / n)
    try:
        script_dir = Path(__file__).parent
        for n in names:
            candidates.append(script_dir / n)
    except NameError:
        pass
    for n in names:
        candidates.append(Path.cwd() / n)

    for p in candidates:
        try:
            if p.exists():
                return p
        except Exception:
            continue
    return None


def get_output_dirs(base_dir) -> dict:
    if base_dir is None:
        raise ValueError("Output directory not set.")
    base_dir = Path(base_dir)
    dirs = {
        "root":     base_dir,
        "data":     base_dir / "Data",
        "pictures": base_dir / "Pictures",
        "logs":     base_dir / "Logs",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


# ===========================================================================
# 6)  STEP 1  -  DESeq2 pipeline
# ===========================================================================
def read_one_count_file(fpath: Path):
    if not fpath.exists():
        return None
    open_func = gzip.open if fpath.suffix == ".gz" else open
    try:
        with open_func(fpath, "rt", encoding="utf-8", errors="replace") as f:
            raw_lines = [f.readline() for _ in range(50)]
    except Exception:
        return None

    header_row_idx = None
    for idx, line in enumerate(raw_lines):
        if line and "gene_id" in line.lower():
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
    df = df[~df.duplicated(subset="gene_id", keep="first")]

    counts = pd.to_numeric(df["unstranded"], errors="coerce")
    mask   = counts.notna()
    df, counts = df[mask], counts[mask]
    counts.index = df["gene_id"].values
    return counts, df.set_index("gene_id")["gene_name"]


def normalize_group(raw_series, log):
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
    log(f"Sample groups: Tumor={ (mapped=='Tumor').sum() }, "
        f"Normal={ (mapped=='Normal').sum() }, Other={ (mapped=='Other').sum() }")
    return mapped


def find_case_id_column(columns):
    cols_lower = {c.lower().strip(): c for c in columns}
    for cand in ("case id", "case ids", "caseid", "case_id"):
        if cand in cols_lower:
            return cols_lower[cand]
    for c, orig in cols_lower.items():
        if "case" in c:
            return orig
    return None


def _build_inference(n_cpus, log):
    try:
        from pydeseq2.default_inference import DefaultInference
        return DefaultInference(n_cpus=n_cpus), {}
    except Exception as e:
        log(f"DefaultInference not available ({e}); falling back to n_cpus kwarg.")
        return None, {"n_cpus": n_cpus}


@dataclass
class PipelineParams:
    counts_root:       Path
    sample_sheet_path: Path
    output_dir:        Path
    label:             str = "TCGA-Project"
    min_count:         int = 10
    log2fc_thr:        float = 1.0
    padj_thr:          float = 0.05


def run_deg_pipeline(params, log=print, progress_callback=None):
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds  import DeseqStats
    except ImportError as e:
        raise RuntimeError(f"PyDESeq2 is not installed. Error: {e}")

    def upd(pct, msg=""):
        if progress_callback:
            progress_callback(pct, msg)

    counts_root       = Path(params.counts_root)
    sample_sheet_path = Path(params.sample_sheet_path)
    dirs              = get_output_dirs(params.output_dir)
    data_dir          = dirs["data"]
    label             = params.label

    upd(2, "Reading sample sheet")
    log(f"Reading sample sheet: {sample_sheet_path}")
    sample_sheet = pd.read_csv(sample_sheet_path, sep="\t", dtype=str)
    sample_sheet.columns = [c.strip() for c in sample_sheet.columns]

    need_cols = {"File ID", "File Name", "Tissue Type", "Sample ID"}
    lower_map = {c.lower(): c for c in sample_sheet.columns}
    missing = []
    for v in need_cols:
        if v.lower() not in lower_map:
            missing.append(v)
        else:
            sample_sheet.rename(columns={lower_map[v.lower()]: v}, inplace=True)
    if missing:
        raise ValueError(f"Missing columns from sample sheet: {missing}")

    case_id_col = find_case_id_column(sample_sheet.columns)
    if case_id_col is None:
        log("[WARN] No 'Case ID' column found - using 'Sample ID' as Case ID.")
        sample_sheet["Case ID"] = sample_sheet["Sample ID"]
    elif case_id_col != "Case ID":
        sample_sheet.rename(columns={case_id_col: "Case ID"}, inplace=True)

    sample_sheet["Sample ID"] = sample_sheet["Sample ID"].astype(str).str.strip()
    log(f"Sample sheet rows: {len(sample_sheet)}")

    dup_mask = sample_sheet.duplicated(subset="Sample ID", keep="first")
    if dup_mask.any():
        sample_sheet[dup_mask].to_csv(
            data_dir / f"{label}_{STD_NAMES['dupes']}", index=False
        )
        log(f"[WARN] {dup_mask.sum()} duplicate Sample IDs - keeping first only.")
        sample_sheet = sample_sheet[~dup_mask].reset_index(drop=True)

    log(f"Reading {len(sample_sheet)} count files...")
    count_dict, gene_name_map = {}, {}
    valid_samples, failed_files, all_gene_ids = [], [], set()
    total = len(sample_sheet)

    for i, row in sample_sheet.iterrows():
        fpath  = counts_root / row["File ID"] / row["File Name"]
        result = read_one_count_file(fpath)
        if result is None:
            failed_files.append(row["Sample ID"])
            continue
        counts, gnames = result
        count_dict[row["Sample ID"]] = dict(zip(counts.index, counts.values))
        all_gene_ids.update(counts.index)
        for gid, gname in gnames.items():
            gene_name_map.setdefault(gid, gname)
        valid_samples.append(row["Sample ID"])

        pct = 5 + int((i + 1) / total * 15)
        if (i + 1) % 25 == 0 or (i + 1) == total:
            upd(pct, f"Reading files: {i + 1}/{total}")
            log(f"  {i + 1}/{total} files read")

    if failed_files:
        log(f"[WARN] {len(failed_files)} files failed and were skipped.")
    if not valid_samples:
        raise RuntimeError("No count files were read successfully.")

    all_gene_ids = sorted(all_gene_ids)
    log(f"Total unique gene IDs: {len(all_gene_ids)}")

    upd(22, "Building count matrix")
    log("Building count matrix...")
    count_matrix = pd.DataFrame.from_dict(count_dict, orient="index").T
    count_matrix = count_matrix.reindex(all_gene_ids).fillna(0).astype(int)
    count_matrix = count_matrix[valid_samples]

    sheet_valid = sample_sheet[sample_sheet["Sample ID"].isin(valid_samples)].copy()
    metadata = sheet_valid[["Sample ID", "Case ID", "Tissue Type"]].copy()
    metadata = metadata.set_index("Sample ID")
    metadata.rename(columns={"Tissue Type": "Group_raw"}, inplace=True)
    metadata["Group"] = normalize_group(metadata["Group_raw"], log)
    metadata = metadata[metadata["Group"].isin(["Tumor", "Normal"])]

    common = list(set(count_matrix.columns) & set(metadata.index))
    count_matrix = count_matrix[common]
    metadata     = metadata.loc[common]

    log(f"Final matrix: {count_matrix.shape[1]} samples x {count_matrix.shape[0]} genes")

    group_counts = metadata["Group"].value_counts()
    if len(group_counts) < 2 or group_counts.min() < 2:
        raise ValueError(
            f"Need at least 2 samples in BOTH groups. Current: {group_counts.to_dict()}"
        )
    metadata.to_csv(data_dir / f"{label}_{STD_NAMES['groups']}")

    upd(26, "Filtering low-count genes")
    min_group_size = int(metadata["Group"].value_counts().min())
    keep_genes = (count_matrix >= params.min_count).sum(axis=1) >= min_group_size
    log(f"Genes before filter: {count_matrix.shape[0]} | after: {int(keep_genes.sum())}")
    count_matrix_f = count_matrix.loc[keep_genes]

    counts_for_deseq = count_matrix_f.T
    counts_for_deseq.index = counts_for_deseq.index.astype(str)
    metadata.index = metadata.index.astype(str)
    meta_deseq = metadata[["Group"]].copy()

    if not counts_for_deseq.index.equals(meta_deseq.index):
        common_idx = counts_for_deseq.index.intersection(meta_deseq.index)
        counts_for_deseq = counts_for_deseq.loc[common_idx]
        meta_deseq       = meta_deseq.loc[common_idx]

    upd(30, "Running PyDESeq2")
    log(f"Running PyDESeq2 (capped at {SAFE_N_CPUS} CPU core(s))...")
    inference, extra_kwargs = _build_inference(SAFE_N_CPUS, log)

    dds_kwargs = dict(
        counts=counts_for_deseq,
        metadata=meta_deseq,
        design_factors="Group",
        refit_cooks=True,
    )
    if inference is not None:
        dds_kwargs["inference"] = inference
    else:
        dds_kwargs.update(extra_kwargs)

    dds = DeseqDataSet(**dds_kwargs)
    upd(40, "DESeq2 fitting")
    dds.deseq2()

    upd(72, "Computing statistics")
    stats_kwargs = dict(contrast=["Group", "Tumor", "Normal"])
    if inference is not None:
        stats_kwargs["inference"] = inference
    else:
        stats_kwargs.update(extra_kwargs)
    stat_res = DeseqStats(dds, **stats_kwargs)
    stat_res.summary()

    upd(86, "Saving results")
    res = stat_res.results_df.copy().reset_index().rename(columns={"index": "gene_id"})
    res["gene_name"] = res["gene_id"].map(gene_name_map) if gene_name_map else res["gene_id"]
    res = res.sort_values("padj")

    res_clean = res.dropna(subset=["padj"])
    degs = res_clean[
        (res_clean["log2FoldChange"].abs() >= params.log2fc_thr)
        & (res_clean["padj"] < params.padj_thr)
    ]
    degs_up   = degs[degs["log2FoldChange"] > 0]
    degs_down = degs[degs["log2FoldChange"] < 0]

    log("-" * 60)
    log("SUMMARY")
    log(f"  Total genes tested : {len(res_clean)}")
    log(f"  DEGs (|log2FC|>={params.log2fc_thr}, padj<{params.padj_thr}): {len(degs)}")
    log(f"     Upregulated    : {len(degs_up)}")
    log(f"     Downregulated  : {len(degs_down)}")
    log("-" * 60)

    paths = {
        "all":    data_dir / f"{label}_{STD_NAMES['all']}",
        "cache":  data_dir / f"{label}_{STD_NAMES['cache']}",
        "sig":    data_dir / f"{label}_{STD_NAMES['sig']}",
        "up":     data_dir / f"{label}_{STD_NAMES['up']}",
        "down":   data_dir / f"{label}_{STD_NAMES['down']}",
        "norm":   data_dir / f"{label}_{STD_NAMES['norm']}",
        "groups": data_dir / f"{label}_{STD_NAMES['groups']}",
    }

    res.to_csv(paths["all"],   index=False)
    res.to_csv(paths["cache"], index=False)
    degs.to_csv(paths["sig"],  index=False)
    degs_up.to_csv(paths["up"],     index=False)
    degs_down.to_csv(paths["down"], index=False)

    normalized_counts = pd.DataFrame(
        dds.layers["normed_counts"],
        index=counts_for_deseq.index,
        columns=counts_for_deseq.columns,
    )
    normalized_counts.to_csv(paths["norm"])

    gene_list = res_clean[
        ["gene_id", "gene_name", "log2FoldChange", "padj", "pvalue", "baseMean"]
    ].copy()
    gene_list["regulation"] = "NS"
    gene_list.loc[degs_up.index,   "regulation"] = "Up"
    gene_list.loc[degs_down.index, "regulation"] = "Down"
    gene_list = gene_list.sort_values("padj")
    gl_path = data_dir / f"{label}_{STD_NAMES['gene_list']}"
    gene_list.to_csv(gl_path, index=False)
    paths["gene_list"] = gl_path

    upd(100, "Complete")
    log(f"\nAll result files saved to: {data_dir}")
    return paths


# ===========================================================================
# 7)  DEG helpers
# ===========================================================================
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
        lfc_col  = pick_col(df, ["log2FoldChange", "log2FC", "logFC"])
        padj_col = pick_col(df, ["padj", "FDR", "adj.P.Val", "qvalue"])
        name_col = pick_col(df, ["gene_name", "gene_symbol", "hgnc_symbol", "gene_id"])
        if lfc_col is None:
            raise ValueError("Could not find a log2FoldChange column.")
        if padj_col is None:
            raise ValueError("Could not find a padj / FDR column.")
        if name_col is None:
            raise ValueError("Could not find a gene name column.")
        out = pd.DataFrame({
            "gene":   df[name_col].astype(str),
            "log2FC": pd.to_numeric(df[lfc_col],  errors="coerce"),
            "sig":    pd.to_numeric(df[padj_col], errors="coerce"),
        })

    out["log2FC"] = pd.to_numeric(out["log2FC"], errors="coerce")
    out["sig"]    = pd.to_numeric(out["sig"],    errors="coerce")
    out = out.dropna(subset=["log2FC", "sig"])
    out["neglog10"] = -np.log10(out["sig"].clip(lower=1e-300))
    out["class"] = "NS"
    out.loc[(out["log2FC"] >=  log2fc_thr) & (out["sig"] < padj_thr), "class"] = "Up"
    out.loc[(out["log2FC"] <= -log2fc_thr) & (out["sig"] < padj_thr), "class"] = "Down"
    return out


def parse_gene_list(raw_text):
    parts = re.split(r"[,\n;\t\s]+", raw_text or "")
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
    exact_set   = set(gene_series)
    lower_map   = {}
    for g in gene_series:
        lower_map.setdefault(g.lower(), g)

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

    up   = deg_df[deg_df["class"] == "Up"].sort_values("sig")
    down = deg_df[deg_df["class"] == "Down"].sort_values("sig")

    lines = [
        f"Analysis summary - {base_name}",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Thresholds: |log2FC| >= {log2fc_thr}, padj < {padj_thr}",
        "",
        f"Total genes analysed : {len(deg_df)}",
        f"Upregulated          : {len(up)}",
        f"Downregulated        : {len(down)}",
        f"Not significant      : {(deg_df['class']=='NS').sum()}",
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


# ===========================================================================
# 8)  Gene labelling  -  clean margin labels with connecting lines
# ===========================================================================
def _add_margin_labels(ax, points_df, x_col, y_col, gene_col, log=print,
                       fontsize=9):
    """Draw labels stacked in a clean column on the RIGHT margin, each connected
    to its data point by a thin line. Text z-order is very high so it always
    sits on top of the scatter points."""
    if points_df.empty:
        return

    pts = points_df.copy()
    pts["_x"] = pd.to_numeric(pts[x_col], errors="coerce")
    pts["_y"] = pd.to_numeric(pts[y_col], errors="coerce")
    pts = pts.dropna(subset=["_x", "_y"]).sort_values("_y", ascending=False).reset_index(drop=True)
    n = len(pts)
    if n == 0:
        return

    # Highlight rings around chosen points
    ax.scatter(pts["_x"], pts["_y"], s=90, facecolors="none",
               edgecolors="black", linewidths=1.4, zorder=15)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    x_range = xmax - xmin
    y_range = ymax - ymin
    if x_range == 0 or y_range == 0:
        return

    # Reserve label column on the right, INSIDE the axes so it doesn't clip
    label_x = xmax - 0.02 * x_range      # right-anchored inside axes
    top     = ymax - 0.03 * y_range
    bottom  = ymin + 0.05 * y_range
    if n == 1:
        ys = np.array([(top + bottom) / 2])
    else:
        ys = np.linspace(top, bottom, n)

    for i, row in pts.iterrows():
        gx, gy = row["_x"], row["_y"]
        ly = ys[i]

        # Connector line (behind label, in front of scatter data)
        ax.plot([gx, label_x], [gy, ly],
                color="#111111", linewidth=0.7,
                alpha=0.85, zorder=16, solid_capstyle="round")

        # Label with rounded white bbox
        ax.text(label_x, ly, str(row[gene_col]),
                fontsize=fontsize, fontweight="bold",
                ha="right", va="center", zorder=20,
                bbox=dict(boxstyle="round,pad=0.28,rounding_size=0.35",
                          facecolor="white",
                          edgecolor="#111111",
                          linewidth=0.8, alpha=0.97))

    log(f"Labelled {n} gene(s).")


# ===========================================================================
# 9)  Plotting functions
# ===========================================================================
def plot_volcano(ax, deg_df, log2fc_thr, padj_thr, palette, title,
                 label_genes=None, log_func=print):
    label_genes = label_genes or []
    display = {"Down": "Down-regulated", "NS": "Not Significant", "Up": "Up-regulated"}

    for cls in ("Down", "NS", "Up"):
        sub = deg_df[deg_df["class"] == cls]
        ax.scatter(sub["log2FC"], sub["neglog10"], s=18, c=palette[cls],
                   alpha=0.78, linewidths=0, label=display[cls], zorder=3)

    ax.axvline( log2fc_thr, ls="--", lw=1.2, color=Theme.PLOT_CUTOFF, zorder=2)
    ax.axvline(-log2fc_thr, ls="--", lw=1.2, color=Theme.PLOT_CUTOFF, zorder=2)
    ax.axhline(-np.log10(padj_thr), ls="--", lw=1.2, color=Theme.PLOT_CUTOFF, zorder=2)

    ax.set_xlabel("Log$_{2}$ Fold Change", fontsize=13)
    ax.set_ylabel("-Log$_{10}$(FDR)",       fontsize=13)
    ax.set_title(title, loc="left", fontsize=16, fontweight="bold", pad=12)
    apply_journal_style(ax)
    ax.legend(title="Status", loc="upper center",
              bbox_to_anchor=(0.5, -0.13), ncol=3,
              fontsize=11, title_fontsize=12, markerscale=1.8, frameon=False)

    if label_genes:
        targets = deg_df[deg_df["gene"].astype(str).isin(label_genes)]
        if not targets.empty:
            _add_margin_labels(ax, targets, "log2FC", "neglog10", "gene",
                               log=log_func, fontsize=9)


def plot_ma(ax, df_raw, log2fc_thr, padj_thr, palette, title,
            label_genes=None, log_func=print):
    lfc_col  = pick_col(df_raw, ["log2FoldChange", "log2FC", "logFC"])
    padj_col = pick_col(df_raw, ["padj", "FDR", "adj.P.Val", "qvalue"])
    mean_col = pick_col(df_raw, ["baseMean", "AveExpr", "meanExpr"])
    if lfc_col is None or padj_col is None or mean_col is None:
        raise ValueError(
            "MA plot needs baseMean / log2FC / padj columns. "
            "Load the analysis cache file (not the filtered file)."
        )

    d = pd.DataFrame({
        "mean":   pd.to_numeric(df_raw[mean_col], errors="coerce"),
        "log2FC": pd.to_numeric(df_raw[lfc_col],  errors="coerce"),
        "padj":   pd.to_numeric(df_raw[padj_col], errors="coerce"),
    })
    d = d[d["mean"] > 0].dropna()
    d["class"] = "NS"
    d.loc[(d["log2FC"] >=  log2fc_thr) & (d["padj"] < padj_thr), "class"] = "Up"
    d.loc[(d["log2FC"] <= -log2fc_thr) & (d["padj"] < padj_thr), "class"] = "Down"

    display = {"Down": "Down-regulated", "NS": "Not Significant", "Up": "Up-regulated"}
    for cls in ("Down", "NS", "Up"):
        sub = d[d["class"] == cls]
        ax.scatter(sub["mean"], sub["log2FC"], s=15, c=palette[cls],
                   alpha=0.78, linewidths=0, label=display[cls], zorder=3)

    ax.set_xscale("log")
    ax.axhline(0, color=Theme.PLOT_CUTOFF, lw=1.2, ls="--", zorder=2)
    ax.set_xlabel("Mean of normalised counts", fontsize=13)
    ax.set_ylabel("Log$_{2}$ Fold Change",     fontsize=13)
    ax.set_title(title, loc="left", fontsize=16, fontweight="bold", pad=12)
    apply_journal_style(ax)
    ax.legend(title="Status", loc="upper center",
              bbox_to_anchor=(0.5, -0.13), ncol=3,
              fontsize=11, title_fontsize=12, markerscale=1.8, frameon=False)

    if label_genes:
        gene_col = pick_col(df_raw, ["gene_name", "gene_id"])
        if gene_col is None:
            return
        df_pos    = df_raw.loc[d.index]
        gene_names = df_pos[gene_col].astype(str)
        targets   = df_pos.loc[gene_names.isin(label_genes)]
        if not targets.empty:
            tdf = pd.DataFrame({
                "gene":   targets[gene_col],
                "mean":   targets[mean_col],
                "log2FC": targets[lfc_col],
            })
            _add_margin_labels(ax, tdf, "mean", "log2FC", "gene",
                               log=log_func, fontsize=9)


def plot_summary_bar(ax, up_n, down_n, title):
    cats   = ["Upregulated", "Downregulated"]
    vals   = [up_n, down_n]
    colors = [Theme.PLOT_UP, Theme.PLOT_DOWN]
    bars = ax.bar(cats, vals, color=colors, width=0.55,
                  edgecolor="#333333", linewidth=1, zorder=3)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}",
                ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.set_ylabel("Number of genes", fontsize=13)
    ax.set_title(title, loc="left", fontsize=16, fontweight="bold", pad=12)
    apply_journal_style(ax)
    ymax = max(vals) if vals else 1
    ax.set_ylim(0, ymax * 1.18 + 1)


def plot_heatmap_expression(ax, norm_counts, sample_groups, gene_ids,
                            gene_labels, title):
    """Expression heatmap. Robust to string/int index mismatch."""
    # Coerce to strings so ID matching is consistent
    norm_counts = norm_counts.copy()
    norm_counts.index = norm_counts.index.astype(str)
    gene_ids_str = [str(g) for g in gene_ids]

    present = [g for g in gene_ids_str if g in norm_counts.index]
    if not present:
        ax.text(0.5, 0.5, "No matching genes for heatmap.\n"
                          "Check that gene_ids in Significant DEGs\n"
                          "match the index of normalized_counts.csv.",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=11, color="#666666")
        ax.set_title(title, loc="left", fontsize=15, fontweight="bold", pad=10)
        ax.set_xticks([]); ax.set_yticks([])
        return

    mat = np.log2(norm_counts.loc[present].clip(lower=0) + 1)
    row_std = mat.std(axis=1).replace(0, 1)
    z = mat.sub(mat.mean(axis=1), axis=0).div(row_std, axis=0)

    # Ensure sample_groups index is string
    sample_groups = sample_groups.copy()
    sample_groups.index = sample_groups.index.astype(str)
    z.columns = z.columns.astype(str)

    if "Group" in sample_groups.columns:
        order = [s for s in sample_groups.sort_values("Group").index if s in z.columns]
    else:
        order = list(z.columns)
    if not order:
        order = list(z.columns)
    z = z[order]

    row_labels = [str(gene_labels.get(g, g)) for g in z.index]

    sns.heatmap(
        z, ax=ax, cmap="RdBu_r", center=0,
        cbar_kws={"label": "Z-score", "shrink": 0.7},
        yticklabels=row_labels, xticklabels=False,
        linewidths=0, rasterized=True,
    )
    ax.set_title(title, loc="left", fontsize=15, fontweight="bold", pad=10)
    ax.set_xlabel(f"Samples (n={len(order)})", fontsize=12)
    ax.set_ylabel("")


def plot_heatmap_rank(ax, deg_df, top_n, title):
    d = deg_df[deg_df["class"] != "NS"].copy()
    if d.empty:
        ax.text(0.5, 0.5, "No significant genes for heatmap.",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=11, color="#666666")
        ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=10)
        ax.set_xticks([]); ax.set_yticks([])
        return
    d = (d.reindex(d["log2FC"].abs().sort_values(ascending=False).index)
           .head(top_n).sort_values("log2FC"))
    colors = [Theme.PLOT_UP if v > 0 else Theme.PLOT_DOWN for v in d["log2FC"]]
    ax.barh(d["gene"], d["log2FC"], color=colors,
            edgecolor="#333333", linewidth=0.6, zorder=3)
    ax.axvline(0, color="#333333", lw=1)
    ax.set_xlabel("Log$_{2}$ Fold Change", fontsize=13)
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=10)
    apply_journal_style(ax)


# ===========================================================================
# 10)  Plot dispatcher  -  makes ONE figure per plot type
# ===========================================================================
def make_figure_for_plot_type(plot_type, ctx):
    """Return a matplotlib Figure for the given plot_type using values from
    the context dict `ctx`. `ctx` contains: main_df, deg_df, sig_df,
    norm_counts, sample_groups, log2fc, padj, palette, label, top_n,
    label_genes (list), log_func."""
    fig = Figure(figsize=(ctx["fig_w"], ctx["fig_h"]), dpi=100, facecolor="white")
    ax  = fig.add_subplot(111)
    label   = ctx["label"]
    palette = ctx["palette"]
    log2fc, padj = ctx["log2fc"], ctx["padj"]

    if plot_type == "Volcano Plot":
        deg = ctx["deg_df"]
        genes = []
        if ctx["label_genes"]:
            found, notf = find_matching_genes(deg, ctx["label_genes"])
            if notf:
                ctx["log_func"](f"[WARN] Not found (Volcano): {', '.join(notf)}")
            genes = found
        plot_volcano(ax, deg, log2fc, padj, palette,
                     f"Volcano Plot  -  {label} (Tumor vs Normal)",
                     label_genes=genes, log_func=ctx["log_func"])

    elif plot_type == "MA Plot":
        df = ctx["main_df"]
        genes = []
        if ctx["label_genes"]:
            gene_col = pick_col(df, ["gene_name", "gene_id"])
            if gene_col:
                gs = df[gene_col].astype(str)
                f, notf = find_matching_genes(pd.DataFrame({"gene": gs}), ctx["label_genes"])
                if notf:
                    ctx["log_func"](f"[WARN] Not found (MA): {', '.join(notf)}")
                genes = f
        plot_ma(ax, df, log2fc, padj, palette,
                f"MA Plot  -  {label}",
                label_genes=genes, log_func=ctx["log_func"])

    elif plot_type == "Summary Bar Chart":
        deg = ctx["deg_df"]
        up_n   = int((deg["class"] == "Up").sum())
        down_n = int((deg["class"] == "Down").sum())
        plot_summary_bar(ax, up_n, down_n,
                         f"Differentially Expressed Genes  -  {label}")

    elif plot_type == "Heatmap (Top DEGs)":
        top_n = ctx["top_n"]
        sig_df = ctx.get("sig_df")
        norm_counts = ctx.get("norm_counts")
        sample_groups = ctx.get("sample_groups")
        if sig_df is not None and norm_counts is not None and sample_groups is not None:
            id_col   = pick_col(sig_df, ["gene_id"])
            name_col = pick_col(sig_df, ["gene_name", "gene_id"])
            top = sig_df.copy()
            padj_col = pick_col(top, ["padj", "FDR", "pvalue"])
            if padj_col:
                top = top.sort_values(padj_col)
            top = top.head(top_n)
            if id_col:
                gene_ids = top[id_col].astype(str).values
                gene_labels = (dict(zip(top[id_col].astype(str),
                                        top[name_col].astype(str)))
                               if name_col else {})
            else:
                gene_ids = top[name_col].astype(str).values
                gene_labels = {}
            plot_heatmap_expression(ax, norm_counts, sample_groups, gene_ids,
                                    gene_labels, f"Top {top_n} DEGs  -  {label}")
        else:
            plot_heatmap_rank(ax, ctx["deg_df"], top_n, f"Top {top_n} DEGs  -  {label}")

    else:
        ax.text(0.5, 0.5, f"Unknown plot type: {plot_type}",
                ha="center", va="center", transform=ax.transAxes)

    fig.tight_layout()
    return fig


def make_combined_figure(ctx):
    """Build a wide 1x3 combined figure: Volcano + MA + Bar."""
    label   = ctx["label"]
    palette = ctx["palette"]
    log2fc, padj = ctx["log2fc"], ctx["padj"]
    label_genes = ctx.get("label_genes", [])
    deg = ctx["deg_df"]
    df  = ctx["main_df"]

    # Wide layout, kept tall enough for legends below
    fig = Figure(figsize=(20, 6.5), dpi=100, facecolor="white")
    axes = fig.subplots(1, 3)

    # 1) Volcano
    volcano_genes = []
    if label_genes:
        found, notf = find_matching_genes(deg, label_genes)
        volcano_genes = found
        if notf:
            ctx["log_func"](f"[WARN] Not found (Volcano): {', '.join(notf)}")
    plot_volcano(axes[0], deg, log2fc, padj, palette,
                 f"Volcano  -  {label}",
                 label_genes=volcano_genes, log_func=ctx["log_func"])

    # 2) MA
    ma_genes = []
    if label_genes:
        gene_col = pick_col(df, ["gene_name", "gene_id"])
        if gene_col:
            gs = df[gene_col].astype(str)
            f, notf = find_matching_genes(pd.DataFrame({"gene": gs}), label_genes)
            ma_genes = f
            if notf:
                ctx["log_func"](f"[WARN] Not found (MA): {', '.join(notf)}")
    plot_ma(axes[1], df, log2fc, padj, palette,
            f"MA Plot  -  {label}",
            label_genes=ma_genes, log_func=ctx["log_func"])

    # 3) Bar
    up_n   = int((deg["class"] == "Up").sum())
    down_n = int((deg["class"] == "Down").sum())
    plot_summary_bar(axes[2], up_n, down_n, f"DEGs  -  {label}")

    fig.tight_layout()
    return fig


# ===========================================================================
# 11)  Modern UI components
# ===========================================================================
class Card(ctk.CTkFrame):
    """Rounded card with optional title header."""
    def __init__(self, master, title="", subtitle="", **kw):
        super().__init__(master, corner_radius=12, fg_color=COL_CARD,
                         border_width=1, border_color=COL_BORDER, **kw)
        if title:
            head = ctk.CTkFrame(self, fg_color="transparent")
            head.pack(fill="x", padx=18, pady=(14, 2))
            ctk.CTkLabel(head, text=title,
                         font=ctk.CTkFont(size=15, weight="bold"),
                         text_color=COL_TEXT, anchor="w").pack(fill="x")
            if subtitle:
                ctk.CTkLabel(head, text=subtitle,
                             font=ctk.CTkFont(size=11),
                             text_color=COL_TEXT_SEC,
                             anchor="w", justify="left",
                             wraplength=580).pack(fill="x", pady=(2, 0))
        self.body = ctk.CTkFrame(self, fg_color="transparent")
        self.body.pack(fill="both", expand=True, padx=18, pady=(8, 14))


# ---------------------------------------------------------------------------
#  Clipboard helper  -  make copy/paste/cut/select-all work on ANY text widget
#  regardless of keyboard layout (fixes Ctrl+C/V not working with non-English
#  layouts like Persian, Russian, Arabic, etc.) and adds a right-click menu.
# ---------------------------------------------------------------------------
def enable_clipboard(widget):
    """Attach robust clipboard shortcuts + right-click context menu to a
    tk.Text, tk.Entry, or their CTk equivalents (CTkTextbox / CTkEntry).

    Works with any keyboard layout because it binds to *both* the layout-
    dependent letter events (<Control-c>) and the layout-independent keycode
    events (<Control-Key>) that fire even when the current layout is Persian
    or another non-Latin script."""

    # CTkTextbox / CTkEntry wrap the underlying tk widget in `_textbox` /
    # `_entry`. Bind on the actual tk widget so events land correctly.
    tk_widget = widget
    if hasattr(widget, "_textbox"):
        tk_widget = widget._textbox
    elif hasattr(widget, "_entry"):
        tk_widget = widget._entry

    is_text_widget = tk_widget.winfo_class() in ("Text",)

    # --- Actions --------------------------------------------------------
    def _copy(*_):
        try:
            tk_widget.event_generate("<<Copy>>")
        except Exception:
            pass
        return "break"

    def _paste(*_):
        try:
            # Delete current selection first (so paste replaces it)
            if is_text_widget:
                try:
                    tk_widget.delete("sel.first", "sel.last")
                except tk.TclError:
                    pass
            else:
                try:
                    tk_widget.delete("sel.first", "sel.last")
                except tk.TclError:
                    pass
            tk_widget.event_generate("<<Paste>>")
        except Exception:
            pass
        return "break"

    def _cut(*_):
        try:
            tk_widget.event_generate("<<Cut>>")
        except Exception:
            pass
        return "break"

    def _select_all(*_):
        try:
            if is_text_widget:
                tk_widget.tag_add("sel", "1.0", "end-1c")
                tk_widget.mark_set("insert", "1.0")
                tk_widget.see("insert")
            else:
                tk_widget.select_range(0, "end")
                tk_widget.icursor("end")
        except Exception:
            pass
        return "break"

    # --- Keyboard bindings (layout-safe) --------------------------------
    # Bind both the letter events (English layout) and the low-level Key
    # events (any layout) so shortcuts work in Persian too.
    tk_widget.bind("<Control-c>", _copy)
    tk_widget.bind("<Control-C>", _copy)
    tk_widget.bind("<Control-v>", _paste)
    tk_widget.bind("<Control-V>", _paste)
    tk_widget.bind("<Control-x>", _cut)
    tk_widget.bind("<Control-X>", _cut)
    tk_widget.bind("<Control-a>", _select_all)
    tk_widget.bind("<Control-A>", _select_all)

    # Fallback: layout-independent keycode handling (fires even if layout
    # doesn't map a key to 'c', 'v', 'x', 'a')
    def _keycode(event):
        if not (event.state & 0x0004):     # Control modifier
            return
        try:
            kc = event.keycode
        except Exception:
            return
        # Windows keycodes: C=67, V=86, X=88, A=65
        if kc == 67:                    return _copy()
        elif kc == 86:                  return _paste()
        elif kc == 88:                  return _cut()
        elif kc == 65:                  return _select_all()

    tk_widget.bind("<Control-Key>", _keycode, add="+")

    # --- Right-click context menu ---------------------------------------
    menu = tk.Menu(tk_widget, tearoff=0)
    menu.add_command(label="Cut",       accelerator="Ctrl+X", command=_cut)
    menu.add_command(label="Copy",      accelerator="Ctrl+C", command=_copy)
    menu.add_command(label="Paste",     accelerator="Ctrl+V", command=_paste)
    menu.add_separator()
    menu.add_command(label="Select All", accelerator="Ctrl+A", command=_select_all)

    def _show_menu(event):
        try:
            # Focus the widget first so keyboard actions apply to it
            tk_widget.focus_set()
            menu.tk_popup(event.x_root, event.y_root)
        finally:
            menu.grab_release()
        return "break"

    tk_widget.bind("<Button-3>", _show_menu)          # Windows / Linux
    tk_widget.bind("<Button-2>", _show_menu)          # macOS trackpad
    tk_widget.bind("<Control-Button-1>", _show_menu)  # macOS ctrl-click


class FileRow(ctk.CTkFrame):
    """Label + entry + Browse button."""
    def __init__(self, master, label, variable, browse_kind="file",
                 filetypes=None, required=False, hint="", **kw):
        super().__init__(master, fg_color="transparent", **kw)
        self.variable = variable
        self.browse_kind = browse_kind
        self.filetypes = filetypes or [("All files", "*.*")]

        head = ctk.CTkFrame(self, fg_color="transparent")
        head.pack(fill="x")
        ctk.CTkLabel(
            head, text=label + ("  *" if required else ""),
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=(Theme.L_ACCENT if required else Theme.L_TEXT_PRIMARY,
                        Theme.ACCENT_TEXT if required else Theme.TEXT_PRIMARY),
            anchor="w",
        ).pack(side="left")
        if hint:
            ctk.CTkLabel(head, text=f"  {hint}",
                         font=ctk.CTkFont(size=10),
                         text_color=COL_TEXT_MUT,
                         anchor="w").pack(side="left")

        row = ctk.CTkFrame(self, fg_color="transparent")
        row.pack(fill="x", pady=(4, 8))
        self.entry = ctk.CTkEntry(
            row, textvariable=variable, height=36, corner_radius=8,
            fg_color=COL_INPUT, border_color=COL_BORDER,
            text_color=COL_TEXT,
        )
        self.entry.pack(side="left", fill="x", expand=True)
        ctk.CTkButton(row, text="Browse", width=90, height=36,
                      corner_radius=8, command=self._browse,
                      fg_color=COL_ACCENT, hover_color=COL_ACCENT_H,
                      text_color=("white", Theme.BG_DARK),
                      font=ctk.CTkFont(size=12, weight="bold"),
                      ).pack(side="left", padx=(8, 0))

    def _browse(self):
        if self.browse_kind == "dir":
            p = filedialog.askdirectory()
        else:
            p = filedialog.askopenfilename(filetypes=self.filetypes)
        if p:
            self.variable.set(p)


# ===========================================================================
# 12)  Welcome dialog
# ===========================================================================
class WelcomeDialog(ctk.CTkToplevel):
    def __init__(self, parent, icon_image=None):
        super().__init__(parent)
        self.title("Welcome")
        self.resizable(False, False)
        self.transient(parent)
        self.configure(fg_color=COL_BG)

        w, h = 560, 500
        self.update_idletasks()
        px = parent.winfo_x() + (parent.winfo_width() - w) // 2
        py = parent.winfo_y() + (parent.winfo_height() - h) // 2
        self.geometry(f"{w}x{h}+{max(px,0)}+{max(py,0)}")

        wrap = ctk.CTkFrame(self, fg_color="transparent")
        wrap.pack(fill="both", expand=True, padx=32, pady=24)

        # Logo (if icon.png exists)
        if icon_image is not None:
            ctk.CTkLabel(wrap, image=icon_image, text="").pack(pady=(0, 12))

        ctk.CTkLabel(
            wrap, text=APP_NAME,
            font=ctk.CTkFont(size=22, weight="bold"),
            text_color=COL_TEXT,
        ).pack(pady=(0, 4))
        ctk.CTkLabel(
            wrap, text=f"Version {APP_VERSION}   .   by {APP_AUTHOR}",
            font=ctk.CTkFont(size=12),
            text_color=COL_TEXT_SEC,
        ).pack()

        ctk.CTkLabel(
            wrap,
            text=("A modern differential expression analysis toolkit\n"
                  "for TCGA-style RNA-seq data, powered by PyDESeq2."),
            font=ctk.CTkFont(size=12),
            text_color=COL_TEXT,
            justify="center",
        ).pack(pady=(16, 8))

        ctk.CTkLabel(
            wrap,
            text=("If you use this tool in your research, please cite\n"
                  "and star the GitHub repository."),
            font=ctk.CTkFont(size=11),
            text_color=COL_TEXT_SEC,
            justify="center",
        ).pack(pady=(0, 14))

        btnrow = ctk.CTkFrame(wrap, fg_color="transparent")
        btnrow.pack()
        ctk.CTkButton(btnrow, text="GitHub", width=140, height=36,
                      corner_radius=10,
                      fg_color=COL_ACCENT, hover_color=COL_ACCENT_H,
                      text_color=("white", Theme.BG_DARK),
                      font=ctk.CTkFont(size=12, weight="bold"),
                      command=lambda: webbrowser.open(APP_GITHUB)
                      ).pack(side="left", padx=8)
        ctk.CTkButton(btnrow, text="LinkedIn", width=120, height=36,
                      corner_radius=10,
                      fg_color="transparent", border_width=1,
                      border_color=COL_ACCENT, text_color=COL_ACCENT,
                      hover_color=COL_ACCENT_S,
                      font=ctk.CTkFont(size=12, weight="bold"),
                      command=lambda: webbrowser.open(APP_LINKEDIN)
                      ).pack(side="left", padx=6)
        ctk.CTkButton(btnrow, text="ORCID", width=120, height=36,
                      corner_radius=10,
                      fg_color="transparent", border_width=1,
                      border_color=COL_ACCENT, text_color=COL_ACCENT,
                      hover_color=COL_ACCENT_S,
                      font=ctk.CTkFont(size=12, weight="bold"),
                      command=lambda: webbrowser.open(APP_ORCID)
                      ).pack(side="left", padx=6)

        ctk.CTkButton(wrap, text="Get Started", width=200, height=40,
                      corner_radius=10,
                      fg_color=Theme.SUCCESS, hover_color=Theme.SUCCESS_HOVER,
                      text_color="white",
                      font=ctk.CTkFont(size=13, weight="bold"),
                      command=self.destroy
                      ).pack(pady=(18, 0))

        self.after(150, self._make_modal)

    def _make_modal(self):
        try:
            self.grab_set()
            self.focus_set()
        except Exception:
            pass


# ===========================================================================
# 13)  Main Application
# ===========================================================================
class DEGPipelineApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        # Force dark by default (scientific dark palette)
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        self.title(f"{APP_NAME}  -  v{APP_VERSION}")
        self.configure(fg_color=COL_BG)

        # ---- Load icons -------------------------------------------------
        # Two separate images:
        #   * icon.png       -> the app icon (window title bar / taskbar / exe)
        #   * main_icon.png  -> the big banner logo (sidebar header + welcome)
        self._icon_path      = get_icon_path()
        self._main_icon_path = get_main_icon_path()
        self._header_image_ctk = None     # CTkImage for sidebar header
        self._welcome_image_ctk = None    # CTkImage for welcome dialog
        self._load_icons()

        # ---- Geometry ---------------------------------------------------
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        w, h = min(1460, sw - 60), min(940, sh - 80)
        self.geometry(f"{w}x{h}+{(sw - w)//2}+{(sh - h)//2}")
        self.minsize(1180, 740)

        # ---- State ------------------------------------------------------
        self.log_queue = queue.Queue()

        # Preview navigation state
        self.preview_plots: List[Dict] = []   # [{'fig','name','filename_key','savable','skip_reason'}]
        self.preview_figs:  List[Figure] = []
        self.preview_names: List[str]    = []
        self.preview_idx = 0

        # ---- Layout: sidebar + main -------------------------------------
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self._build_sidebar()
        self._build_main_area()

        # ---- Enable clipboard shortcuts on ALL text widgets -------------
        # (Ctrl+C / V / X / A works with any keyboard layout, right-click
        # menu is also added.)
        self.after(200, self._enable_clipboard_on_all)

        # ---- Log queue --------------------------------------------------
        self.after(120, self._poll_log_queue)
        # ---- Welcome dialog --------------------------------------------
        self.after(350, self._show_welcome)

    def _enable_clipboard_on_all(self, root=None):
        """Walk the widget tree and enable clipboard shortcuts on every
        tk.Text / tk.Entry (and their CTk wrappers)."""
        if root is None:
            root = self
        try:
            cls = root.winfo_class()
            if cls in ("Text", "Entry", "TEntry"):
                try:
                    enable_clipboard(root)
                except Exception:
                    pass
        except Exception:
            pass
        try:
            for child in root.winfo_children():
                self._enable_clipboard_on_all(child)
        except Exception:
            pass

    # ================================================================
    # Icon
    # ================================================================
    def _load_icons(self):
        """Load both icon.png (window icon) and main_icon.png (header logo).

        * icon.png       -> set as the window/taskbar icon
        * main_icon.png  -> scaled to CTkImages used in the sidebar header
                            and the welcome dialog
        Both fall back gracefully if a file is missing."""

        # ============================================================
        # 1)  Window / taskbar icon (from icon.png)
        # ============================================================
        if self._icon_path is not None:
            try:
                pil_win = Image.open(str(self._icon_path))
                if pil_win.mode not in ("RGBA", "RGB"):
                    pil_win = pil_win.convert("RGBA")

                if os.name == "nt":
                    ico_path = None
                    if self._icon_path.suffix.lower() == ".ico":
                        ico_path = self._icon_path
                    else:
                        try:
                            ico_path = self._icon_path.with_suffix(".ico")
                            if not ico_path.exists():
                                sq = pil_win.copy()
                                w, h = sq.size
                                side = min(w, h)
                                left = (w - side) // 2
                                top  = (h - side) // 2
                                sq = sq.crop((left, top, left + side, top + side))
                                sq.save(ico_path, format="ICO",
                                        sizes=[(16, 16), (32, 32), (48, 48),
                                               (64, 64), (128, 128), (256, 256)])
                        except Exception as e:
                            print(f"[icon] Could not build .ico: {e}")
                            ico_path = None

                    if ico_path and ico_path.exists():
                        try:
                            self.iconbitmap(str(ico_path))
                        except Exception as e:
                            print(f"[icon] iconbitmap failed: {e}")

                # Cross-platform iconphoto (works on Linux/Mac too)
                try:
                    import tkinter as _tk
                    self._photo_for_wm = _tk.PhotoImage(file=str(self._icon_path))
                    self.iconphoto(True, self._photo_for_wm)
                except Exception as e:
                    print(f"[icon] iconphoto failed: {e}")
            except Exception as e:
                print(f"[icon] Could not load icon.png: {e}")

        # ============================================================
        # 2)  Sidebar-header + welcome banner (from main_icon.png)
        # ============================================================
        if self._main_icon_path is not None:
            try:
                pil_main = Image.open(str(self._main_icon_path))
                if pil_main.mode not in ("RGBA", "RGB"):
                    pil_main = pil_main.convert("RGBA")

                # Sidebar header (wide banner, fits nicely in ~220px sidebar)
                header_pil = pil_main.copy()
                header_pil.thumbnail((210, 140), Image.LANCZOS)
                self._header_image_ctk = ctk.CTkImage(
                    light_image=header_pil, dark_image=header_pil,
                    size=header_pil.size,
                )

                # Welcome dialog (bigger)
                welcome_pil = pil_main.copy()
                welcome_pil.thumbnail((420, 240), Image.LANCZOS)
                self._welcome_image_ctk = ctk.CTkImage(
                    light_image=welcome_pil, dark_image=welcome_pil,
                    size=welcome_pil.size,
                )
            except Exception as e:
                print(f"[icon] Could not load main_icon.png: {e}")

    # ================================================================
    # Welcome
    # ================================================================
    def _show_welcome(self):
        WelcomeDialog(self, icon_image=self._welcome_image_ctk)

    # ================================================================
    # Sidebar
    # ================================================================
    def _build_sidebar(self):
        # Sidebar is a fixed-width column laid out with pack() so items flow
        # naturally from top (brand) to middle (nav) to bottom (theme + links).
        SIDEBAR_W = 260
        side = ctk.CTkFrame(self, width=SIDEBAR_W, corner_radius=0,
                            fg_color=COL_CARD, border_width=0)
        side.grid(row=0, column=0, sticky="nsw")
        side.pack_propagate(False)

        # ---------- TOP: brand banner --------------------------------
        brand = ctk.CTkFrame(side, fg_color="transparent")
        brand.pack(fill="x", padx=16, pady=(22, 6))

        if self._header_image_ctk is not None:
            ctk.CTkLabel(brand, image=self._header_image_ctk, text=""
                         ).pack(pady=(0, 8))
        else:
            ctk.CTkLabel(brand, text="DEG",
                         font=ctk.CTkFont(size=28, weight="bold"),
                         text_color=COL_ACCENT).pack()

        ctk.CTkLabel(brand, text="DEG Pipeline",
                     font=ctk.CTkFont(size=16, weight="bold"),
                     text_color=COL_TEXT).pack(pady=(2, 0))
        ctk.CTkLabel(brand, text=f"v{APP_VERSION}",
                     font=ctk.CTkFont(size=11),
                     text_color=COL_TEXT_MUT).pack()

        # ---------- Divider ------------------------------------------
        ctk.CTkFrame(side, height=1, fg_color=COL_BORDER
                     ).pack(fill="x", padx=18, pady=(16, 12))

        # ---------- Section label ------------------------------------
        ctk.CTkLabel(side, text="  NAVIGATION",
                     font=ctk.CTkFont(size=10, weight="bold"),
                     text_color=COL_TEXT_MUT,
                     anchor="w").pack(fill="x", padx=18, pady=(0, 6))

        # ---------- Nav buttons --------------------------------------
        nav_wrap = ctk.CTkFrame(side, fg_color="transparent")
        nav_wrap.pack(fill="x", padx=12, pady=(0, 6))

        self.nav_buttons = {}
        for key, text in [
            ("step1", "  Step 1  -  DEG Analysis"),
            ("step2", "  Step 2  -  Visualisation"),
            ("about", "  About  &  Tutorial"),
        ]:
            btn = ctk.CTkButton(
                nav_wrap, text=text, anchor="w", height=42, corner_radius=10,
                fg_color="transparent", text_color=COL_TEXT,
                hover_color=COL_HOVER,
                font=ctk.CTkFont(size=13, weight="bold"),
                command=lambda k=key: self._select_page(k),
            )
            btn.pack(fill="x", pady=3)
            self.nav_buttons[key] = btn

        # ---------- Spacer (pushes bottom section down) --------------
        ctk.CTkFrame(side, fg_color="transparent").pack(fill="both", expand=True)

        # ---------- BOTTOM: divider + appearance + links -------------
        ctk.CTkFrame(side, height=1, fg_color=COL_BORDER
                     ).pack(fill="x", padx=18, pady=(6, 12))

        bottom = ctk.CTkFrame(side, fg_color="transparent")
        bottom.pack(fill="x", padx=16, pady=(0, 16))

        ctk.CTkLabel(bottom, text="APPEARANCE",
                     font=ctk.CTkFont(size=10, weight="bold"),
                     text_color=COL_TEXT_MUT,
                     anchor="w").pack(fill="x", pady=(0, 4))

        seg = ctk.CTkSegmentedButton(
            bottom, values=["Light", "Dark", "System"],
            command=self._change_appearance,
            fg_color=COL_INPUT,
            selected_color=COL_ACCENT,
            selected_hover_color=COL_ACCENT_H,
            unselected_color=COL_INPUT,
            unselected_hover_color=COL_HOVER,
            text_color=COL_TEXT,
            font=ctk.CTkFont(size=11, weight="bold"),
            height=30,
        )
        seg.set("Dark")
        seg.pack(fill="x", pady=(0, 12))

        # Quick-links: full-width stacked buttons (clean, no cramping)
        ctk.CTkLabel(bottom, text="QUICK LINKS",
                     font=ctk.CTkFont(size=10, weight="bold"),
                     text_color=COL_TEXT_MUT,
                     anchor="w").pack(fill="x", pady=(0, 4))
        for text, url in [
            ("GitHub Repository", APP_GITHUB),
            ("LinkedIn Profile",  APP_LINKEDIN),
            ("ORCID Profile",     APP_ORCID),
        ]:
            ctk.CTkButton(bottom, text=text, height=30, corner_radius=8,
                          anchor="w",
                          fg_color="transparent", border_width=1,
                          border_color=COL_BORDER,
                          text_color=COL_TEXT,
                          hover_color=COL_HOVER,
                          font=ctk.CTkFont(size=11, weight="bold"),
                          command=lambda u=url: webbrowser.open(u),
                          ).pack(fill="x", pady=2)

    def _change_appearance(self, mode):
        ctk.set_appearance_mode(mode.lower())

    def _select_page(self, key):
        for k, b in self.nav_buttons.items():
            if k == key:
                b.configure(fg_color=COL_ACCENT_S,
                            text_color=COL_ACCENT)
            else:
                b.configure(fg_color="transparent",
                            text_color=COL_TEXT)
        self.pages[key].tkraise()

    def _show_about(self):
        messagebox.showinfo(
            f"About  -  {APP_NAME}",
            f"{APP_NAME}  v{APP_VERSION}\n"
            f"Author: {APP_AUTHOR}\n\n"
            "A modern differential expression analysis toolkit for\n"
            "TCGA-style RNA-seq data, powered by PyDESeq2.\n\n"
            f"GitHub  : {APP_GITHUB}\n"
            f"LinkedIn: {APP_LINKEDIN}\n"
            f"ORCID   : {APP_ORCID}"
        )

    # ================================================================
    # Main area
    # ================================================================
    def _build_main_area(self):
        container = ctk.CTkFrame(self, fg_color=COL_BG, corner_radius=0)
        container.grid(row=0, column=1, sticky="nsew")
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.pages = {}
        for key in ("step1", "step2", "about"):
            page = ctk.CTkFrame(container, fg_color="transparent")
            page.grid(row=0, column=0, sticky="nsew")
            self.pages[key] = page

        self._build_step1_page(self.pages["step1"])
        self._build_step2_page(self.pages["step2"])
        self._build_about_page(self.pages["about"])
        self._select_page("step1")

    # ================================================================
    # STEP 1
    # ================================================================
    def _build_step1_page(self, parent):
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(1, weight=1)

        header = ctk.CTkFrame(parent, fg_color="transparent")
        header.grid(row=0, column=0, sticky="ew", padx=28, pady=(22, 8))
        ctk.CTkLabel(header, text="Step 1  -  DEG Analysis",
                     font=ctk.CTkFont(size=26, weight="bold"),
                     text_color=COL_TEXT, anchor="w").pack(anchor="w")
        ctk.CTkLabel(
            header,
            text=("Point to your Gene Expression Quantification folder "
                  "and GDC sample sheet, then run PyDESeq2. Results are "
                  "saved into Data / Pictures / Logs. Plots are generated "
                  "in Step 2."),
            font=ctk.CTkFont(size=12),
            text_color=COL_TEXT_SEC,
            anchor="w", justify="left", wraplength=900,
        ).pack(anchor="w", pady=(4, 0))

        body = ctk.CTkScrollableFrame(parent, fg_color="transparent")
        body.grid(row=1, column=0, sticky="nsew", padx=20, pady=(8, 20))
        body.grid_columnconfigure(0, weight=1)

        self.s1_counts_dir   = tk.StringVar()
        self.s1_sample_sheet = tk.StringVar()
        self.s1_output_dir   = tk.StringVar()
        self.s1_label        = tk.StringVar(value="TCGA-Project")
        self.s1_min_count    = tk.StringVar(value="10")
        self.s1_log2fc       = tk.StringVar(value="1.0")
        self.s1_padj         = tk.StringVar(value="0.05")

        # Inputs card
        inputs = Card(body, title="Inputs",
                      subtitle="All three paths are required to run the analysis.")
        inputs.pack(fill="x", padx=8, pady=(6, 10))

        FileRow(inputs.body,
                "Gene Expression Quantification folder",
                self.s1_counts_dir, browse_kind="dir", required=True,
                hint="root folder with per-sample subfolders"
                ).pack(fill="x")
        FileRow(inputs.body,
                "GDC Sample Sheet (.tsv)",
                self.s1_sample_sheet, browse_kind="file", required=True,
                filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")],
                hint="the GDC sample sheet you downloaded"
                ).pack(fill="x")
        FileRow(inputs.body,
                "Output folder",
                self.s1_output_dir, browse_kind="dir", required=True,
                hint="Data / Pictures / Logs created inside"
                ).pack(fill="x")

        # Parameters
        params = Card(body, title="Analysis parameters")
        params.pack(fill="x", padx=8, pady=10)

        grid = ctk.CTkFrame(params.body, fg_color="transparent")
        grid.pack(fill="x")
        grid.grid_columnconfigure((0, 1, 2, 3), weight=1)

        def _mini(parent_, label, var, col):
            wrap = ctk.CTkFrame(parent_, fg_color="transparent")
            wrap.grid(row=0, column=col, sticky="ew", padx=6)
            ctk.CTkLabel(wrap, text=label,
                         font=ctk.CTkFont(size=11, weight="bold"),
                         text_color=COL_TEXT, anchor="w").pack(fill="x")
            ctk.CTkEntry(wrap, textvariable=var, height=36,
                         corner_radius=8,
                         fg_color=COL_INPUT, border_color=COL_BORDER,
                         text_color=COL_TEXT,
                         ).pack(fill="x", pady=(4, 0))

        _mini(grid, "Project label",              self.s1_label,     0)
        _mini(grid, "Min. read count per gene",   self.s1_min_count, 1)
        _mini(grid, "|log2FC| threshold",         self.s1_log2fc,    2)
        _mini(grid, "Adj. p-value (FDR)",         self.s1_padj,      3)

        # Output overview
        overview = Card(body, title="Output files",
                        subtitle="These files are produced in <output>/Data/")
        overview.pack(fill="x", padx=8, pady=10)
        for key in ("all", "cache", "sig", "up", "down", "norm", "groups", "gene_list"):
            note = ""
            if key == "cache":     note = "    (load THIS into Step 2)"
            if key == "gene_list": note = "    (publication-ready gene list)"
            row = ctk.CTkFrame(overview.body, fg_color="transparent")
            row.pack(fill="x", pady=1)
            ctk.CTkLabel(row, text=">",
                         font=ctk.CTkFont(size=13, weight="bold"),
                         text_color=COL_ACCENT).pack(side="left", padx=(0, 8))
            ctk.CTkLabel(row,
                         text=f"<label>_{STD_NAMES[key]}{note}",
                         font=ctk.CTkFont(size=11, family="Consolas"),
                         text_color=COL_TEXT_SEC,
                         anchor="w").pack(side="left", fill="x", expand=True)

        # Run + progress
        run_card = Card(body, title="Run")
        run_card.pack(fill="x", padx=8, pady=10)

        runrow = ctk.CTkFrame(run_card.body, fg_color="transparent")
        runrow.pack(fill="x")
        self.s1_run_btn = ctk.CTkButton(
            runrow, text="  Run DEG Analysis",
            height=44, corner_radius=10,
            font=ctk.CTkFont(size=13, weight="bold"),
            fg_color=Theme.SUCCESS, hover_color=Theme.SUCCESS_HOVER,
            text_color="white",
            command=self._start_step1,
        )
        self.s1_run_btn.pack(side="left", padx=(0, 8))
        ctk.CTkButton(runrow, text="Copy log", height=44, width=110,
                      corner_radius=10,
                      fg_color="transparent", border_width=1,
                      border_color=COL_BORDER, text_color=COL_TEXT,
                      hover_color=COL_HOVER,
                      command=self._copy_s1_log).pack(side="left", padx=4)
        ctk.CTkButton(runrow, text="Clear log", height=44, width=110,
                      corner_radius=10,
                      fg_color="transparent", border_width=1,
                      border_color=COL_BORDER, text_color=COL_TEXT,
                      hover_color=COL_HOVER,
                      command=lambda: (self.s1_log.configure(state="normal"),
                                       self.s1_log.delete("1.0", "end"))
                      ).pack(side="left", padx=4)

        self.progress_label = ctk.CTkLabel(
            run_card.body, text="Ready",
            font=ctk.CTkFont(size=11), text_color=COL_TEXT_SEC, anchor="w")
        self.progress_label.pack(fill="x", pady=(14, 4))
        self.progress_bar = ctk.CTkProgressBar(
            run_card.body, height=14, corner_radius=8,
            progress_color=COL_ACCENT, fg_color=COL_INPUT,
        )
        self.progress_bar.pack(fill="x")
        self.progress_bar.set(0)

        # Log
        log_card = Card(body, title="Progress log")
        log_card.pack(fill="both", expand=True, padx=8, pady=(10, 6))
        self.s1_log = ctk.CTkTextbox(
            log_card.body, height=260, corner_radius=8, wrap="word",
            font=ctk.CTkFont(size=11, family="Consolas"),
            fg_color=COL_INPUT, text_color=COL_TEXT,
            border_width=1, border_color=COL_BORDER,
        )
        self.s1_log.pack(fill="both", expand=True)

    # ---- STEP 1 helpers ---------------------------------------------
    def _log1(self, msg):
        self.s1_log.configure(state="normal")
        self.s1_log.insert("end", str(msg) + "\n")
        self.s1_log.see("end")

    def _copy_s1_log(self):
        text = self.s1_log.get("1.0", "end-1c")
        self.clipboard_clear()
        self.clipboard_append(text)
        messagebox.showinfo("Copied", "Progress log copied to clipboard.")

    def _validate_number(self, var, name, cast=float, min_val=None, max_val=None):
        try:
            v = cast(var.get())
        except Exception:
            raise ValueError(f"'{name}' must be a valid {cast.__name__}.")
        if min_val is not None and v < min_val:
            raise ValueError(f"'{name}' must be >= {min_val}.")
        if max_val is not None and v > max_val:
            raise ValueError(f"'{name}' must be <= {max_val}.")
        return v

    def _start_step1(self):
        errors = []
        if not self.s1_counts_dir.get().strip():   errors.append("Counts folder is required.")
        if not self.s1_sample_sheet.get().strip(): errors.append("Sample sheet is required.")
        if not self.s1_output_dir.get().strip():   errors.append("Output folder is required.")
        if not self.s1_label.get().strip():        errors.append("Project label is required.")

        try:
            min_count = int(self._validate_number(self.s1_min_count, "Min read count",
                                                  cast=int, min_val=0))
            log2fc    = self._validate_number(self.s1_log2fc, "|log2FC| threshold",
                                              min_val=0, max_val=20)
            padj      = self._validate_number(self.s1_padj, "padj threshold",
                                              min_val=1e-10, max_val=1)
        except ValueError as e:
            errors.append(str(e))

        counts_root = Path(self.s1_counts_dir.get().strip()) if self.s1_counts_dir.get() else None
        sheet_path  = Path(self.s1_sample_sheet.get().strip()) if self.s1_sample_sheet.get() else None
        if counts_root and not counts_root.exists():
            errors.append(f"Counts folder does not exist:\n  {counts_root}")
        if sheet_path and not sheet_path.exists():
            errors.append(f"Sample sheet does not exist:\n  {sheet_path}")

        if errors:
            messagebox.showerror("Invalid input", "\n\n".join(errors))
            return

        self.s1_run_btn.configure(state="disabled", text="  Running...")
        self.s1_log.configure(state="normal")
        self.s1_log.delete("1.0", "end")
        self.progress_bar.set(0)
        self.progress_bar.configure(progress_color=COL_ACCENT)
        self.progress_label.configure(text="Starting...")

        params = PipelineParams(
            counts_root       = counts_root,
            sample_sheet_path = sheet_path,
            output_dir        = Path(self.s1_output_dir.get().strip()),
            label             = self.s1_label.get().strip().replace(" ", "_"),
            min_count         = min_count,
            log2fc_thr        = log2fc,
            padj_thr          = padj,
        )
        threading.Thread(
            target=self._run_step1_worker, args=(params,), daemon=True
        ).start()

    def _run_step1_worker(self, params):
        log_lines = []

        def qlog(msg):
            log_lines.append(str(msg))
            self.log_queue.put(("s1_log", msg))

        def progress_cb(pct, msg=""):
            self.log_queue.put(("s1_prog", (pct, msg)))

        try:
            paths = run_deg_pipeline(params, log=qlog, progress_callback=progress_cb)
            try:
                dirs = get_output_dirs(params.output_dir)
                (dirs["logs"] / f"{params.label}_{STD_NAMES['log']}"
                 ).write_text("\n".join(log_lines), encoding="utf-8")
            except Exception:
                pass
            self.log_queue.put(("s1_done", paths))
        except Exception as e:
            qlog("[FATAL] " + str(e))
            qlog(traceback.format_exc())
            self.log_queue.put(("s1_failed", str(e)))

    # ================================================================
    # Queue pump
    # ================================================================
    def _poll_log_queue(self):
        try:
            for _ in range(120):
                kind, payload = self.log_queue.get_nowait()
                if kind == "s1_log":
                    self._log1(payload)
                elif kind == "s1_prog":
                    pct, msg = payload
                    self.progress_bar.set(pct / 100.0)
                    self.progress_label.configure(text=f"{msg}   ({pct:.0f}%)")
                elif kind == "s1_done":
                    self._on_step1_done(payload)
                elif kind == "s1_failed":
                    self.s1_run_btn.configure(state="normal", text="  Run DEG Analysis")
                    self.progress_bar.configure(progress_color=Theme.DANGER)
                    self.progress_label.configure(text="Failed - see log")
                    messagebox.showerror("Analysis failed", str(payload))
                elif kind == "s2_log":
                    self._log2(payload)
                elif kind == "s2_prog":
                    pct, msg = payload
                    self.s2_progress.set(pct / 100.0)
                    self.s2_progress_label.configure(text=f"{msg}   ({pct:.0f}%)")
                elif kind == "s2_previews_ready":
                    self._on_previews_ready(payload)
                elif kind == "s2_previews_failed":
                    self._on_previews_failed(payload)
                elif kind == "s2_save_done":
                    self._on_save_done(payload)
                elif kind == "s2_save_failed":
                    self._on_save_failed(payload)
        except queue.Empty:
            pass
        self.after(100, self._poll_log_queue)

    def _on_step1_done(self, paths):
        self.s1_run_btn.configure(state="normal", text="  Run DEG Analysis")
        self.progress_bar.set(1.0)
        self.progress_bar.configure(progress_color=Theme.SUCCESS)
        self.progress_label.configure(text="Complete!   (100%)")
        self._log1("\nDONE. Result files saved. Switch to Step 2 to generate plots.\n")

        if paths:
            self.s2_files["all"].set(str(paths.get("cache", "")))
            self.s2_files["sig"].set(str(paths.get("sig",   "")))
            self.s2_files["norm"].set(str(paths.get("norm",  "")))
            self.s2_files["groups"].set(str(paths.get("groups", "")))
            self.s2_output_dir.set(str(Path(paths.get("all", "")).parent.parent))
            self.s2_label.set(self.s1_label.get())
            self.s2_log2fc.set(self.s1_log2fc.get())
            self.s2_padj.set(self.s1_padj.get())

            if messagebox.askyesno(
                "Analysis complete",
                "DEG analysis finished successfully.\n\n"
                "The cache file has been auto-loaded into Step 2.\n"
                "Do you want to switch to Step 2 now to generate plots?",
            ):
                self._select_page("step2")

    # ================================================================
    # STEP 2  -  build page
    # ================================================================
    def _build_step2_page(self, parent):
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        # State vars
        self.s2_files = {
            "all":    tk.StringVar(),
            "sig":    tk.StringVar(),
            "norm":   tk.StringVar(),
            "groups": tk.StringVar(),
        }
        self.s2_output_dir    = tk.StringVar()
        self.s2_label         = tk.StringVar(value="TCGA-Project")
        self.plot_type        = tk.StringVar(value="All Plots")   # DEFAULT
        self.s2_log2fc        = tk.StringVar(value="1.0")
        self.s2_padj          = tk.StringVar(value="0.05")
        self.top_n            = tk.StringVar(value="20")
        self.palette_choice   = tk.StringVar(value="Reference (red/blue/grey)")
        self.fig_w            = tk.StringVar(value="9.0")
        self.fig_h            = tk.StringVar(value="7.0")
        self.dpi              = tk.StringVar(value="600")
        self.fmt_png          = tk.BooleanVar(value=True)
        self.fmt_pdf          = tk.BooleanVar(value=True)
        self.fmt_tiff         = tk.BooleanVar(value=False)
        self.fmt_svg          = tk.BooleanVar(value=False)
        self.export_processed = tk.BooleanVar(value=True)
        # If any genes are listed to label, also generate the unlabeled copy
        self.save_both_variants = tk.BooleanVar(value=True)
        # Also build a "combined" figure with Volcano + MA + Bar side by side
        self.save_combined      = tk.BooleanVar(value=True)

        # ---- Left panel ------------------------------------------------
        left = ctk.CTkScrollableFrame(parent, width=480, fg_color="transparent")
        left.grid(row=0, column=0, sticky="nsw", padx=(20, 10), pady=20)
        left.grid_columnconfigure(0, weight=1)

        ctk.CTkLabel(left, text="Step 2  -  Visualisation",
                     font=ctk.CTkFont(size=24, weight="bold"),
                     text_color=COL_TEXT,
                     anchor="w").pack(fill="x", padx=4, pady=(0, 4))
        ctk.CTkLabel(left,
                     text="Load a cache file, choose a plot type, then "
                          "preview or export publication-ready figures.",
                     font=ctk.CTkFont(size=11),
                     text_color=COL_TEXT_SEC,
                     anchor="w", justify="left", wraplength=450
                     ).pack(fill="x", padx=4, pady=(0, 12))

        # 1. Data
        data = Card(left, title="1.  Load data")
        data.pack(fill="x", padx=4, pady=6)
        FileRow(data.body, "Analysis cache / results (required)",
                self.s2_files["all"], browse_kind="file", required=True,
                filetypes=[("CSV", "*.csv"), ("All", "*.*")],
                hint="drives Volcano, MA, Bar"
                ).pack(fill="x")
        ctk.CTkLabel(data.body,
                     text="Optional  -  needed for expression Heatmap:",
                     font=ctk.CTkFont(size=11, slant="italic"),
                     text_color=COL_TEXT_MUT, anchor="w"
                     ).pack(fill="x", pady=(6, 2))
        FileRow(data.body, "Significant DEGs", self.s2_files["sig"],
                filetypes=[("CSV", "*.csv")]).pack(fill="x")
        FileRow(data.body, "Normalised counts", self.s2_files["norm"],
                filetypes=[("CSV", "*.csv")]).pack(fill="x")
        FileRow(data.body, "Sample grouping", self.s2_files["groups"],
                filetypes=[("CSV", "*.csv")]).pack(fill="x")

        quick = ctk.CTkFrame(data.body, fg_color="transparent")
        quick.pack(fill="x", pady=(4, 0))
        ctk.CTkButton(quick, text="Load from Data folder",
                      corner_radius=8, height=34,
                      fg_color=COL_ACCENT, hover_color=COL_ACCENT_H,
                      text_color=("white", Theme.BG_DARK),
                      font=ctk.CTkFont(size=12, weight="bold"),
                      command=self._load_from_data_folder
                      ).pack(side="left", padx=(0, 6), expand=True, fill="x")
        ctk.CTkButton(quick, text="Load from log file",
                      corner_radius=8, height=34,
                      fg_color="transparent", border_width=1,
                      border_color=COL_ACCENT, text_color=COL_ACCENT,
                      hover_color=COL_ACCENT_S,
                      font=ctk.CTkFont(size=12, weight="bold"),
                      command=self._load_from_log
                      ).pack(side="left", padx=(6, 0), expand=True, fill="x")

        # 2. Plot type
        pcard = Card(left, title="2.  Plot type")
        pcard.pack(fill="x", padx=4, pady=6)
        ctk.CTkOptionMenu(
            pcard.body, variable=self.plot_type,
            values=list(PLOT_TYPES.keys()),
            command=lambda *_: self._on_plot_type_change(),
            height=38, corner_radius=8,
            fg_color=COL_ACCENT, button_color=COL_ACCENT_H,
            button_hover_color=COL_ACCENT_H,
            text_color=("white", Theme.BG_DARK),
            font=ctk.CTkFont(size=12, weight="bold"),
            dropdown_fg_color=COL_INPUT,
            dropdown_text_color=COL_TEXT,
            dropdown_hover_color=COL_HOVER,
        ).pack(fill="x")

        self.plot_desc = ctk.CTkLabel(
            pcard.body, text="", anchor="w", justify="left",
            wraplength=440, font=ctk.CTkFont(size=11),
            text_color=COL_TEXT_SEC,
        )
        self.plot_desc.pack(fill="x", pady=(6, 0))

        # 3. Parameters
        pp = Card(left, title="3.  Parameters")
        pp.pack(fill="x", padx=4, pady=6)

        pg = ctk.CTkFrame(pp.body, fg_color="transparent")
        pg.pack(fill="x")
        pg.grid_columnconfigure((0, 1), weight=1)

        def _pf(parent_, label, var, row, col):
            wrap = ctk.CTkFrame(parent_, fg_color="transparent")
            wrap.grid(row=row, column=col, sticky="ew", padx=4, pady=4)
            ctk.CTkLabel(wrap, text=label,
                         font=ctk.CTkFont(size=11, weight="bold"),
                         text_color=COL_TEXT, anchor="w").pack(fill="x")
            ctk.CTkEntry(wrap, textvariable=var, height=34, corner_radius=8,
                         fg_color=COL_INPUT, border_color=COL_BORDER,
                         text_color=COL_TEXT
                         ).pack(fill="x", pady=(3, 0))

        _pf(pg, "Project label",         self.s2_label,   0, 0)
        _pf(pg, "|log2FC| threshold",    self.s2_log2fc,  0, 1)
        _pf(pg, "Adj. p-value (FDR)",    self.s2_padj,    1, 0)
        _pf(pg, "Top N (heatmap only)",  self.top_n,      1, 1)

        ctk.CTkLabel(pp.body, text="Genes to label on plots",
                     font=ctk.CTkFont(size=11, weight="bold"),
                     text_color=COL_TEXT, anchor="w"
                     ).pack(fill="x", pady=(10, 2))
        ctk.CTkLabel(pp.body,
                     text="One per line, or comma/semicolon separated. Empty = no labels.",
                     font=ctk.CTkFont(size=10),
                     text_color=COL_TEXT_MUT, anchor="w",
                     justify="left", wraplength=440
                     ).pack(fill="x")
        self.gene_text = ctk.CTkTextbox(
            pp.body, height=90, corner_radius=8,
            font=ctk.CTkFont(size=11, family="Consolas"),
            fg_color=COL_INPUT, text_color=COL_TEXT,
            border_width=1, border_color=COL_BORDER,
        )
        self.gene_text.pack(fill="x", pady=(4, 6))

        pr = ctk.CTkFrame(pp.body, fg_color="transparent")
        pr.pack(fill="x", pady=(4, 0))
        ctk.CTkLabel(pr, text="Colour palette",
                     font=ctk.CTkFont(size=11, weight="bold"),
                     text_color=COL_TEXT, anchor="w").pack(side="left")
        ctk.CTkOptionMenu(
            pr, variable=self.palette_choice,
            values=["Reference (red/blue/grey)", "Deep (dark red/navy/grey)"],
            height=32, corner_radius=8,
            fg_color=COL_INPUT, button_color=COL_ACCENT,
            button_hover_color=COL_ACCENT_H,
            text_color=COL_TEXT,
            dropdown_fg_color=COL_INPUT,
            dropdown_text_color=COL_TEXT,
            dropdown_hover_color=COL_HOVER,
        ).pack(side="right", fill="x", expand=True, padx=(10, 0))

        pr2 = ctk.CTkFrame(pp.body, fg_color="transparent")
        pr2.pack(fill="x", pady=(6, 0))
        pr2.grid_columnconfigure((0, 1, 2), weight=1)
        _pf(pr2, "Fig. width (in)",  self.fig_w, 0, 0)
        _pf(pr2, "Fig. height (in)", self.fig_h, 0, 1)
        _pf(pr2, "Export DPI",       self.dpi,   0, 2)

        # 4. Output
        oc = Card(left, title="4.  Output")
        oc.pack(fill="x", padx=4, pady=6)
        FileRow(oc.body, "Output folder", self.s2_output_dir,
                browse_kind="dir",
                hint="Pictures/ and Data/ created inside"
                ).pack(fill="x")

        fmt = ctk.CTkFrame(oc.body, fg_color="transparent")
        fmt.pack(fill="x", pady=(4, 0))
        ctk.CTkLabel(fmt, text="Export formats:",
                     font=ctk.CTkFont(size=11, weight="bold"),
                     text_color=COL_TEXT, anchor="w").pack(anchor="w")
        chks = ctk.CTkFrame(fmt, fg_color="transparent")
        chks.pack(fill="x", pady=(4, 0))
        for label, var in [
            ("PNG",  self.fmt_png),
            ("PDF",  self.fmt_pdf),
            ("TIFF", self.fmt_tiff),
            ("SVG",  self.fmt_svg),
        ]:
            ctk.CTkCheckBox(chks, text=label, variable=var,
                            corner_radius=6, checkbox_width=18, checkbox_height=18,
                            fg_color=COL_ACCENT, hover_color=COL_ACCENT_H,
                            text_color=COL_TEXT, border_color=COL_BORDER_STRONG,
                            ).pack(side="left", padx=8)
        # Extra options
        ctk.CTkLabel(oc.body, text="Extra options:",
                     font=ctk.CTkFont(size=11, weight="bold"),
                     text_color=COL_TEXT, anchor="w"
                     ).pack(anchor="w", pady=(10, 2))
        ctk.CTkCheckBox(oc.body,
                        text="Also save unlabeled copies when genes are listed",
                        variable=self.save_both_variants,
                        corner_radius=6, checkbox_width=18, checkbox_height=18,
                        fg_color=COL_ACCENT, hover_color=COL_ACCENT_H,
                        text_color=COL_TEXT, border_color=COL_BORDER_STRONG,
                        ).pack(anchor="w", pady=(4, 0))
        ctk.CTkCheckBox(oc.body,
                        text="Also build combined figure (Volcano + MA + Bar)",
                        variable=self.save_combined,
                        corner_radius=6, checkbox_width=18, checkbox_height=18,
                        fg_color=COL_ACCENT, hover_color=COL_ACCENT_H,
                        text_color=COL_TEXT, border_color=COL_BORDER_STRONG,
                        ).pack(anchor="w", pady=(4, 0))
        ctk.CTkCheckBox(oc.body,
                        text="Also export processed data CSV + summary TXT",
                        variable=self.export_processed,
                        corner_radius=6, checkbox_width=18, checkbox_height=18,
                        fg_color=COL_ACCENT, hover_color=COL_ACCENT_H,
                        text_color=COL_TEXT, border_color=COL_BORDER_STRONG,
                        ).pack(anchor="w", pady=(4, 0))

        # 5. Actions
        actions = Card(left, title="5.  Actions")
        actions.pack(fill="x", padx=4, pady=6)

        self.s2_preview_btn = ctk.CTkButton(
            actions.body, text="Generate & Preview",
            height=46, corner_radius=10,
            font=ctk.CTkFont(size=13, weight="bold"),
            fg_color=Theme.SUCCESS, hover_color=Theme.SUCCESS_HOVER,
            text_color="white",
            command=self.generate_previews,
        )
        self.s2_preview_btn.pack(fill="x")

        self.s2_save_btn = ctk.CTkButton(
            actions.body, text="Save all previewed plots",
            height=40, corner_radius=10,
            font=ctk.CTkFont(size=12, weight="bold"),
            fg_color=COL_ACCENT, hover_color=COL_ACCENT_H,
            text_color=("white", Theme.BG_DARK),
            command=self.save_all_previews,
        )
        self.s2_save_btn.pack(fill="x", pady=(6, 0))

        # Progress
        self.s2_progress_label = ctk.CTkLabel(
            actions.body, text="Idle",
            font=ctk.CTkFont(size=11), text_color=COL_TEXT_SEC, anchor="w")
        self.s2_progress_label.pack(fill="x", pady=(12, 3))
        self.s2_progress = ctk.CTkProgressBar(
            actions.body, height=10, corner_radius=6,
            progress_color=COL_ACCENT, fg_color=COL_INPUT,
        )
        self.s2_progress.pack(fill="x")
        self.s2_progress.set(0)

        # Log
        log_card = Card(left, title="Log")
        log_card.pack(fill="both", expand=True, padx=4, pady=(6, 4))
        self.s2_log_text = ctk.CTkTextbox(
            log_card.body, height=120, corner_radius=8, wrap="word",
            font=ctk.CTkFont(size=10, family="Consolas"),
            fg_color=COL_INPUT, text_color=COL_TEXT,
            border_width=1, border_color=COL_BORDER,
        )
        self.s2_log_text.pack(fill="both", expand=True)

        # ---- Right panel (preview) --------------------------------------
        right = ctk.CTkFrame(parent, corner_radius=14, fg_color=COL_CARD,
                             border_width=1, border_color=COL_BORDER)
        right.grid(row=0, column=1, sticky="nsew", padx=(10, 20), pady=20)
        right.grid_rowconfigure(1, weight=1)
        right.grid_columnconfigure(0, weight=1)

        # Header with switch controls
        head = ctk.CTkFrame(right, fg_color="transparent")
        head.grid(row=0, column=0, sticky="ew", padx=16, pady=(14, 4))
        head.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(head, text="Preview",
                     font=ctk.CTkFont(size=16, weight="bold"),
                     text_color=COL_TEXT, anchor="w"
                     ).grid(row=0, column=0, sticky="w")

        # Center: title of current plot
        self.preview_title_label = ctk.CTkLabel(
            head, text="No preview yet - click Generate & Preview",
            font=ctk.CTkFont(size=13, weight="bold"),
            text_color=COL_ACCENT, anchor="center"
        )
        self.preview_title_label.grid(row=0, column=1, sticky="ew", padx=10)

        # Right: navigation
        nav = ctk.CTkFrame(head, fg_color="transparent")
        nav.grid(row=0, column=2, sticky="e")

        self.prev_btn = ctk.CTkButton(
            nav, text="<", width=42, height=34, corner_radius=8,
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color=COL_INPUT, hover_color=COL_HOVER,
            text_color=COL_TEXT, border_width=1, border_color=COL_BORDER,
            command=self._preview_prev, state="disabled",
        )
        self.prev_btn.pack(side="left", padx=(0, 4))

        self.preview_counter = ctk.CTkLabel(
            nav, text="0 / 0",
            font=ctk.CTkFont(size=12, weight="bold"),
            text_color=COL_TEXT_SEC, width=60,
        )
        self.preview_counter.pack(side="left", padx=4)

        self.next_btn = ctk.CTkButton(
            nav, text=">", width=42, height=34, corner_radius=8,
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color=COL_INPUT, hover_color=COL_HOVER,
            text_color=COL_TEXT, border_width=1, border_color=COL_BORDER,
            command=self._preview_next, state="disabled",
        )
        self.next_btn.pack(side="left", padx=(4, 0))

        # Thumbnails / tabs row
        self.thumb_row = ctk.CTkFrame(right, fg_color="transparent")
        self.thumb_row.grid(row=1, column=0, sticky="ew", padx=16, pady=(4, 0))
        self.thumb_buttons: List[ctk.CTkButton] = []

        # Matplotlib canvas  (no toolbar - the arrows above are enough
        # and the default toolbar clashes visually with the dark theme)
        self.preview_frame = ctk.CTkFrame(right, corner_radius=10,
                                          fg_color="white",
                                          border_width=1, border_color=COL_BORDER)
        self.preview_frame.grid(row=2, column=0, sticky="nsew",
                                padx=16, pady=(6, 14))
        right.grid_rowconfigure(2, weight=1)

        self.fig = Figure(figsize=(9, 7), dpi=100, facecolor="white")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.preview_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=6, pady=6)

        self._draw_welcome_plot()
        self._on_plot_type_change()

    # ================================================================
    # Preview navigation
    # ================================================================
    def _draw_welcome_plot(self):
        self.fig.clf()
        ax = self.fig.add_subplot(111)
        ax.text(0.5, 0.58,
                APP_NAME,
                ha="center", va="center", fontsize=22, fontweight="bold",
                color="#333333", transform=ax.transAxes)
        ax.text(0.5, 0.45,
                "Load a cache file and click Generate & Preview.",
                ha="center", va="center", fontsize=12,
                color="#666666", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_facecolor("white")
        self.canvas.draw_idle()

    def _log2(self, msg):
        self.s2_log_text.configure(state="normal")
        self.s2_log_text.insert("end", str(msg) + "\n")
        self.s2_log_text.see("end")

    def _on_plot_type_change(self):
        self.plot_desc.configure(text=PLOT_TYPES.get(self.plot_type.get(), ""))

    def _get_palette(self):
        return PALETTE if self.palette_choice.get().startswith("Reference") else PALETTE_ALT

    def _get_label_gene_list(self):
        return parse_gene_list(self.gene_text.get("1.0", "end"))

    def _clear_thumbs(self):
        for b in self.thumb_buttons:
            b.destroy()
        self.thumb_buttons = []

    def _rebuild_thumbs(self):
        self._clear_thumbs()
        for i, name in enumerate(self.preview_names):
            btn = ctk.CTkButton(
                self.thumb_row, text=name, height=30, corner_radius=8,
                fg_color=COL_INPUT if i != self.preview_idx else COL_ACCENT,
                text_color=COL_TEXT if i != self.preview_idx
                           else ("white", Theme.BG_DARK),
                hover_color=COL_HOVER,
                border_width=1, border_color=COL_BORDER,
                font=ctk.CTkFont(size=11, weight="bold"),
                command=lambda idx=i: self._show_preview(idx),
            )
            btn.pack(side="left", padx=3, pady=4)
            self.thumb_buttons.append(btn)

    def _show_preview(self, idx):
        if not self.preview_figs:
            return
        idx = idx % len(self.preview_figs)
        self.preview_idx = idx
        new_fig = self.preview_figs[idx]

        # Swap the canvas figure and redraw
        self.canvas.figure = new_fig
        self.canvas.draw_idle()

        # Update UI state
        self.preview_title_label.configure(text=self.preview_names[idx])
        self.preview_counter.configure(
            text=f"{idx + 1} / {len(self.preview_figs)}")

        state_prev = "normal" if len(self.preview_figs) > 1 else "disabled"
        self.prev_btn.configure(state=state_prev)
        self.next_btn.configure(state=state_prev)

        # Refresh thumb colors
        for i, b in enumerate(self.thumb_buttons):
            if i == idx:
                b.configure(fg_color=COL_ACCENT,
                            text_color=("white", Theme.BG_DARK))
            else:
                b.configure(fg_color=COL_INPUT,
                            text_color=COL_TEXT)

    def _preview_prev(self):
        if not self.preview_figs:
            return
        self._show_preview((self.preview_idx - 1) % len(self.preview_figs))

    def _preview_next(self):
        if not self.preview_figs:
            return
        self._show_preview((self.preview_idx + 1) % len(self.preview_figs))

    # ================================================================
    # STEP 2 quick-load
    # ================================================================
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
            messagebox.showerror(
                "No cache file",
                f"No file matching *{STD_NAMES['cache']} or *{STD_NAMES['all']} "
                "found in that folder.")
            return

        if len(cache_candidates) > 1:
            chosen = filedialog.askopenfilename(
                initialdir=folder, title="Select cache file",
                filetypes=[("CSV", "*.csv")])
            if not chosen:
                return
            cache_file = Path(chosen)
        else:
            cache_file = cache_candidates[0]

        stem = cache_file.stem
        if   stem.endswith("_analysis_cache"):    label = stem[:-len("_analysis_cache")]
        elif stem.endswith("_All_DEGs_results"):  label = stem[:-len("_All_DEGs_results")]
        else:                                     label = "TCGA-Project"

        self.s2_label.set(label)
        self.s2_files["all"].set(str(cache_file))

        for key, fname in [("sig", STD_NAMES["sig"]),
                           ("norm", STD_NAMES["norm"]),
                           ("groups", STD_NAMES["groups"])]:
            cands = list(data_path.glob(f"*{fname}"))
            self.s2_files[key].set(
                str(max(cands, key=lambda p: p.stat().st_mtime))
                if cands else "")

        self.s2_output_dir.set(str(data_path.parent))
        self._log2(f"Loaded data from folder: {folder}")

    def _load_from_log(self):
        log_path = filedialog.askopenfilename(
            title="Select analysis_log.txt",
            filetypes=[("Log files", "*.txt"), ("All files", "*.*")]
        )
        if not log_path:
            return
        try:
            content = Path(log_path).read_text(encoding="utf-8", errors="replace")
            m = re.search(r"All result files saved to:\s*(.+)", content)
            if not m:
                raise ValueError("Could not find 'All result files saved to:' in log.")
            data_dir = Path(m.group(1).strip())
            self.s2_output_dir.set(str(data_dir.parent))

            cache_files = list(data_dir.glob(f"*{STD_NAMES['cache']}"))
            if cache_files:
                cache_file = cache_files[0]
                label = cache_file.stem[:-len("_analysis_cache")]
                self.s2_label.set(label)
                self.s2_files["all"].set(str(cache_file))
                for key, fname in [("sig", STD_NAMES["sig"]),
                                   ("norm", STD_NAMES["norm"]),
                                   ("groups", STD_NAMES["groups"])]:
                    cands = list(data_dir.glob(f"*{fname}"))
                    self.s2_files[key].set(str(cands[0]) if cands else "")
            self._log2("Loaded analysis from log file.")
        except Exception as e:
            messagebox.showerror("Load failed", str(e))

    # ================================================================
    # STEP 2  -  Generate previews  (async, no GUI freeze)
    # ================================================================
    def _selected_plot_types(self) -> List[str]:
        """Return the list of plot types to generate based on the dropdown."""
        chosen = self.plot_type.get()
        if chosen == "All Plots":
            return ["Volcano Plot", "MA Plot", "Summary Bar Chart", "Heatmap (Top DEGs)"]
        return [chosen]

    def _heatmap_data_available(self, ctx) -> bool:
        """True only if we have enough data for a REAL expression heatmap."""
        return (ctx.get("sig_df") is not None
                and ctx.get("norm_counts") is not None
                and ctx.get("sample_groups") is not None
                and not ctx["sig_df"].empty
                and not ctx["norm_counts"].empty)

    def _build_ctx_for_plots(self):
        """Build a shared ctx dict used by the worker to render plots."""
        fig_w  = self._validate_number(self.fig_w, "Figure width",   min_val=1, max_val=50)
        fig_h  = self._validate_number(self.fig_h, "Figure height",  min_val=1, max_val=50)
        log2fc = self._validate_number(self.s2_log2fc, "log2FC threshold", min_val=0, max_val=20)
        padj   = self._validate_number(self.s2_padj,   "padj threshold",   min_val=1e-10, max_val=1)
        top_n  = int(self._validate_number(self.top_n, "Top N",
                                           cast=int, min_val=1, max_val=500))

        main_path = self.s2_files["all"].get().strip()
        if not main_path:
            raise ValueError("Load the Analysis cache / results file first.")

        main_df = load_csv(main_path)
        deg_df  = prep_deg_table(main_df, log2fc, padj)

        # Optional heatmap inputs
        sig_path    = self.s2_files["sig"].get().strip()
        norm_path   = self.s2_files["norm"].get().strip()
        groups_path = self.s2_files["groups"].get().strip()

        sig_df, norm_counts, sample_groups = None, None, None
        if sig_path:
            try:
                sig_df = load_csv(sig_path)
            except Exception as e:
                self._log2(f"[WARN] Could not read Significant DEGs: {e}")
        if norm_path:
            try:
                norm_counts = pd.read_csv(norm_path, index_col=0)
            except Exception as e:
                self._log2(f"[WARN] Could not read normalized counts: {e}")
        if groups_path:
            try:
                sample_groups = pd.read_csv(groups_path, index_col=0)
            except Exception as e:
                self._log2(f"[WARN] Could not read sample groups: {e}")

        return {
            "main_df":       main_df,
            "deg_df":        deg_df,
            "sig_df":        sig_df,
            "norm_counts":   norm_counts,
            "sample_groups": sample_groups,
            "log2fc":        log2fc,
            "padj":          padj,
            "palette":       self._get_palette(),
            "label":         self.s2_label.get().strip() or "Sample",
            "top_n":         top_n,
            "label_genes":   self._get_label_gene_list(),
            "fig_w":         fig_w,
            "fig_h":         fig_h,
            "log_func":      lambda m: self.log_queue.put(("s2_log", m)),
        }

    def generate_previews(self):
        """Kick off asynchronous generation of all requested plots."""
        try:
            ctx = self._build_ctx_for_plots()
        except Exception as e:
            messagebox.showerror("Cannot generate", str(e))
            return

        types = self._selected_plot_types()
        self.s2_preview_btn.configure(state="disabled", text="Generating...")
        self.s2_save_btn.configure(state="disabled")
        self.s2_progress.set(0)
        self.s2_progress.configure(progress_color=COL_ACCENT)
        self.s2_progress_label.configure(text="Preparing...")

        threading.Thread(
            target=self._preview_worker, args=(ctx, types), daemon=True,
        ).start()

    def _preview_worker(self, ctx, types):
        """Build a list of PlotSpec-like dicts:
             {'fig': Figure, 'name': str, 'filename_key': str}
        Handles unlabeled + labeled + combined + heatmap-availability."""
        try:
            plots = []              # list of dicts
            gene_list = ctx["label_genes"]
            has_labels = bool(gene_list)
            want_both = self.save_both_variants.get()
            want_combined = self.save_combined.get()
            hm_available = self._heatmap_data_available(ctx)

            # ---- Plan the work so the progress bar is meaningful --------
            work = []   # list of ('single'|'combined', plot_type or None, labeled_bool)
            for pt in types:
                if pt == "Heatmap (Top DEGs)":
                    # Real expression heatmap only if data is available.
                    # If not available, we still render ONE placeholder for
                    # the preview but skip it during save.
                    work.append(("single", pt, False))
                    continue

                if has_labels:
                    if want_both:
                        work.append(("single", pt, False))   # plain
                        work.append(("single", pt, True))    # labeled
                    else:
                        work.append(("single", pt, True))    # labeled only
                else:
                    work.append(("single", pt, False))

            # Combined figure (Volcano + MA + Bar side by side)
            combined_needed = (
                want_combined
                and all(t in types for t in ("Volcano Plot", "MA Plot", "Summary Bar Chart"))
            )
            if combined_needed:
                if has_labels and want_both:
                    work.append(("combined", None, False))
                    work.append(("combined", None, True))
                else:
                    work.append(("combined", None, has_labels))

            total = max(len(work), 1)

            for i, (kind, pt, labeled) in enumerate(work):
                pct = int(5 + i / total * 90)
                if kind == "single":
                    tag = pt + (" (labeled)" if labeled else "")
                    self.log_queue.put(("s2_prog", (pct, f"Rendering: {tag}")))
                    self.log_queue.put(("s2_log",  f"Rendering: {tag}"))

                    plot_ctx = dict(ctx)
                    if not labeled:
                        plot_ctx["label_genes"] = []

                    # Skip real heatmap if data missing but keep a
                    # placeholder in the preview
                    is_heatmap = (pt == "Heatmap (Top DEGs)")
                    savable = True
                    if is_heatmap and not hm_available:
                        savable = False   # do NOT save this one

                    try:
                        fig = make_figure_for_plot_type(pt, plot_ctx)
                    except Exception as e:
                        self.log_queue.put(("s2_log", f"[ERROR] {tag}: {e}"))
                        continue

                    key = (pt.lower()
                           .replace(" ", "_")
                           .replace("(", "").replace(")", ""))
                    if labeled:
                        key = key + "_labeled"
                    plots.append({
                        "fig": fig,
                        "name": tag,
                        "filename_key": key,
                        "savable": savable,
                        "skip_reason": None if savable
                        else "No expression data (Significant DEGs + "
                             "normalized_counts + sample_grouping required).",
                    })

                elif kind == "combined":
                    tag = "Combined (Volcano + MA + Bar)" + (
                        " (labeled)" if labeled else "")
                    self.log_queue.put(("s2_prog", (pct, f"Rendering: {tag}")))
                    self.log_queue.put(("s2_log",  f"Rendering: {tag}"))
                    plot_ctx = dict(ctx)
                    if not labeled:
                        plot_ctx["label_genes"] = []
                    try:
                        fig = make_combined_figure(plot_ctx)
                        key = "combined_standard" + ("_labeled" if labeled else "")
                        plots.append({
                            "fig": fig,
                            "name": tag,
                            "filename_key": key,
                            "savable": True,
                            "skip_reason": None,
                        })
                    except Exception as e:
                        self.log_queue.put(("s2_log", f"[ERROR] {tag}: {e}"))
                        continue

            self.log_queue.put(("s2_prog", (100, "Done")))
            self.log_queue.put(("s2_previews_ready", (plots, ctx)))
        except Exception as e:
            self.log_queue.put(("s2_previews_failed",
                                (str(e), traceback.format_exc())))

    def _on_previews_ready(self, payload):
        plots, ctx = payload
        # Close any old figures to free memory
        for spec in self.preview_plots:
            try:
                plt.close(spec["fig"])
            except Exception:
                pass

        self.preview_plots = plots
        self.preview_figs  = [p["fig"]  for p in plots]     # kept for compatibility
        self.preview_names = [p["name"] for p in plots]
        self.preview_idx   = 0
        self.last_ctx      = ctx

        self.s2_preview_btn.configure(state="normal", text="Generate & Preview")
        self.s2_save_btn.configure(state="normal")
        self.s2_progress.configure(progress_color=Theme.SUCCESS)
        self.s2_progress.set(1.0)
        savable = sum(1 for p in plots if p["savable"])
        self.s2_progress_label.configure(
            text=f"Ready - {len(plots)} plot(s), {savable} savable")

        if not plots:
            messagebox.showwarning("No plots",
                                   "No plots were generated. See the log.")
            self._draw_welcome_plot()
            self.preview_title_label.configure(text="No preview")
            self.preview_counter.configure(text="0 / 0")
            self._clear_thumbs()
            return

        self._rebuild_thumbs()
        self._show_preview(0)
        self._log2(f"Generated {len(plots)} plot(s).")

    def _on_previews_failed(self, payload):
        msg, tb = payload
        self.s2_preview_btn.configure(state="normal", text="Generate & Preview")
        self.s2_progress.configure(progress_color=Theme.DANGER)
        self.s2_progress_label.configure(text="Failed - see log")
        self._log2(f"[ERROR] Preview generation failed: {msg}")
        self._log2(tb)
        messagebox.showerror("Preview failed", msg)

    # ================================================================
    # STEP 2  -  Save ALL previews  (async)
    # ================================================================
    def save_all_previews(self):
        if not self.preview_plots:
            messagebox.showwarning("Nothing to save",
                                   "Generate & Preview first.")
            return
        outdir = self.s2_output_dir.get().strip()
        if not outdir:
            messagebox.showwarning("Missing output",
                                   "Choose an output folder first.")
            return

        formats = []
        if self.fmt_png.get():  formats.append("png")
        if self.fmt_pdf.get():  formats.append("pdf")
        if self.fmt_tiff.get(): formats.append("tiff")
        if self.fmt_svg.get():  formats.append("svg")
        if not formats:
            messagebox.showwarning("No format", "Select at least one export format.")
            return

        try:
            dpi = int(self._validate_number(self.dpi, "DPI",
                                            cast=int, min_val=50, max_val=2400))
        except ValueError as e:
            messagebox.showerror("Invalid DPI", str(e))
            return

        self.s2_save_btn.configure(state="disabled", text="Saving...")
        self.s2_preview_btn.configure(state="disabled")
        self.s2_progress.set(0)
        self.s2_progress.configure(progress_color=COL_ACCENT)
        self.s2_progress_label.configure(text="Saving...")

        # Prepare data snapshot
        try:
            dirs = get_output_dirs(Path(outdir))
        except Exception as e:
            messagebox.showerror("Output folder", str(e))
            self.s2_save_btn.configure(state="normal", text="Save all previewed plots")
            self.s2_preview_btn.configure(state="normal")
            return

        # Snapshot plots (deep enough for the worker thread)
        plot_specs = list(self.preview_plots)
        label      = (self.s2_label.get().strip().replace(" ", "_") or "sample")
        ctx        = getattr(self, "last_ctx", None)

        threading.Thread(
            target=self._save_worker,
            args=(plot_specs, dirs, label, formats, dpi, ctx),
            daemon=True,
        ).start()

    def _save_worker(self, plot_specs, dirs, label, formats, dpi, ctx):
        saved_files = []
        skipped     = []
        # Filter to savable plots only (heatmap without data is skipped)
        savable = [p for p in plot_specs if p["savable"]]
        n = max(len(savable), 1)
        try:
            # Log the skipped ones first
            for p in plot_specs:
                if not p["savable"]:
                    reason = p.get("skip_reason", "")
                    skipped.append(p["name"])
                    self.log_queue.put((
                        "s2_log",
                        f"[SKIP] {p['name']} - {reason}"))

            for i, spec in enumerate(savable):
                base_name = f"{label}_{spec['filename_key']}"
                self.log_queue.put(("s2_prog",
                                    (int(5 + i / n * 85),
                                     f"Saving: {spec['name']}")))
                self.log_queue.put(("s2_log", f"Saving {spec['name']}..."))
                for fmt in formats:
                    out_path = Path(dirs["pictures"]) / f"{base_name}.{fmt}"
                    try:
                        spec["fig"].savefig(out_path, dpi=dpi,
                                            bbox_inches="tight", format=fmt,
                                            facecolor="white")
                        saved_files.append(str(out_path))
                    except Exception as e:
                        self.log_queue.put(("s2_log",
                                            f"[ERROR] {out_path.name}: {e}"))

            # Optional: processed data + summary from ctx.deg_df
            if self.export_processed.get() and ctx is not None:
                try:
                    dp, sp = save_processed_data_and_summary(
                        ctx["deg_df"], dirs["data"],
                        f"{label}_processed",
                        ctx["log2fc"], ctx["padj"])
                    saved_files += [str(dp), str(sp)]
                except Exception as e:
                    self.log_queue.put(("s2_log",
                                        f"[WARN] Could not write processed data: {e}"))

            self.log_queue.put(("s2_prog", (100, "Saved")))
            self.log_queue.put(("s2_save_done", (saved_files, skipped)))
        except Exception as e:
            self.log_queue.put(("s2_save_failed",
                                (str(e), traceback.format_exc())))

    def _on_save_done(self, payload):
        saved_files, skipped = payload
        self.s2_save_btn.configure(state="normal", text="Save all previewed plots")
        self.s2_preview_btn.configure(state="normal")
        self.s2_progress.configure(progress_color=Theme.SUCCESS)
        self.s2_progress_label.configure(
            text=f"Saved {len(saved_files)} file(s)"
                 + (f", skipped {len(skipped)}" if skipped else ""))
        self._log2(f"Saved {len(saved_files)} file(s):")
        for p in saved_files:
            self._log2(f"  {p}")
        if skipped:
            self._log2(f"Skipped {len(skipped)} plot(s): {', '.join(skipped)}")
        msg = f"Successfully saved {len(saved_files)} file(s) to:\n{self.s2_output_dir.get()}"
        if skipped:
            msg += (f"\n\nSkipped ({len(skipped)}):\n  "
                    + "\n  ".join(skipped)
                    + "\n\n(These need Significant DEGs + normalized_counts "
                      "+ sample_grouping to be produced.)")
        messagebox.showinfo("Saved", msg)

    def _on_save_failed(self, payload):
        msg, tb = payload
        self.s2_save_btn.configure(state="normal", text="Save all previewed plots")
        self.s2_preview_btn.configure(state="normal")
        self.s2_progress.configure(progress_color=Theme.DANGER)
        self.s2_progress_label.configure(text="Save failed - see log")
        self._log2(f"[ERROR] Save failed: {msg}")
        self._log2(tb)
        messagebox.showerror("Save failed", msg)

    # ================================================================
    # ABOUT / TUTORIAL PAGE
    # ================================================================
    def _build_about_page(self, parent):
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        scroll = ctk.CTkScrollableFrame(parent, fg_color="transparent")
        scroll.grid(row=0, column=0, sticky="nsew", padx=28, pady=22)
        scroll.grid_columnconfigure(0, weight=1)

        # ---- Hero card (logo + title + author) --------------------------
        hero = ctk.CTkFrame(scroll, corner_radius=14, fg_color=COL_CARD,
                            border_width=1, border_color=COL_BORDER)
        hero.pack(fill="x", pady=(0, 14))
        hero_in = ctk.CTkFrame(hero, fg_color="transparent")
        hero_in.pack(fill="x", padx=24, pady=22)

        if self._welcome_image_ctk is not None:
            ctk.CTkLabel(hero_in, image=self._welcome_image_ctk, text=""
                         ).pack(pady=(0, 12))

        ctk.CTkLabel(hero_in, text=APP_NAME,
                     font=ctk.CTkFont(size=26, weight="bold"),
                     text_color=COL_TEXT).pack()
        ctk.CTkLabel(hero_in,
                     text=f"Version {APP_VERSION}   .   by {APP_AUTHOR}",
                     font=ctk.CTkFont(size=13),
                     text_color=COL_ACCENT).pack(pady=(4, 0))

        ctk.CTkLabel(hero_in,
                     text=("A modern, publication-quality Differential Expression "
                           "Analysis toolkit for TCGA-style RNA-seq data, powered "
                           "by PyDESeq2."),
                     font=ctk.CTkFont(size=12), text_color=COL_TEXT_SEC,
                     justify="center", wraplength=760,
                     ).pack(pady=(14, 0))

        # ---- Citation card ---------------------------------------------
        cite = Card(scroll, title="Please cite this work",
                    subtitle="Significant effort went into designing this "
                             "tool. If it helped your research, please cite "
                             "it and star the GitHub repo - it truly matters.")
        cite.pack(fill="x", pady=8)

        citation_text = (
            f"Balaei A. ({datetime.now().year}). "
            f"{APP_NAME} (v{APP_VERSION}) [Computer software].\n"
            f"{APP_GITHUB}"
        )
        cite_box = ctk.CTkTextbox(
            cite.body, height=68, corner_radius=8, wrap="word",
            font=ctk.CTkFont(size=11, family="Consolas"),
            fg_color=COL_INPUT, text_color=COL_TEXT,
            border_width=1, border_color=COL_BORDER,
        )
        cite_box.pack(fill="x")
        cite_box.insert("1.0", citation_text)
        cite_box.configure(state="disabled")

        ctk.CTkButton(
            cite.body, text="Copy citation to clipboard",
            height=36, corner_radius=8,
            fg_color=COL_ACCENT, hover_color=COL_ACCENT_H,
            text_color=("white", Theme.BG_DARK),
            font=ctk.CTkFont(size=12, weight="bold"),
            command=lambda: (self.clipboard_clear(),
                             self.clipboard_append(citation_text),
                             messagebox.showinfo("Copied", "Citation copied to clipboard.")),
        ).pack(fill="x", pady=(8, 0))

        # ---- About the author / social links ---------------------------
        author = Card(scroll, title="About the author",
                      subtitle="Alireza Balaei - developer of this toolkit.")
        author.pack(fill="x", pady=8)

        for label, url in [
            ("GitHub Repository", APP_GITHUB),
            ("LinkedIn Profile",  APP_LINKEDIN),
            ("ORCID Profile",     APP_ORCID),
        ]:
            row = ctk.CTkFrame(author.body, fg_color="transparent")
            row.pack(fill="x", pady=3)
            ctk.CTkLabel(row, text=label + ":",
                         font=ctk.CTkFont(size=12, weight="bold"),
                         text_color=COL_TEXT, anchor="w", width=170,
                         ).pack(side="left")
            link_lbl = ctk.CTkLabel(
                row, text=url,
                font=ctk.CTkFont(size=11, family="Consolas", underline=True),
                text_color=COL_ACCENT, anchor="w", cursor="hand2",
            )
            link_lbl.pack(side="left", fill="x", expand=True)
            link_lbl.bind("<Button-1>", lambda e, u=url: webbrowser.open(u))

            ctk.CTkButton(row, text="Open", width=80, height=28,
                          corner_radius=6,
                          fg_color="transparent", border_width=1,
                          border_color=COL_ACCENT, text_color=COL_ACCENT,
                          hover_color=COL_ACCENT_S,
                          font=ctk.CTkFont(size=11, weight="bold"),
                          command=lambda u=url: webbrowser.open(u),
                          ).pack(side="right")

        # ---- Tutorial ---------------------------------------------------
        tut = Card(scroll, title="Tutorial - how to use this tool",
                   subtitle="Follow these steps for a complete "
                            "differential expression analysis.")
        tut.pack(fill="x", pady=8)

        tutorial_steps = [
            ("STEP 1 - Prepare your data",
             "Download RNA-seq data from GDC (TCGA) as 'Gene Expression "
             "Quantification' STAR-Counts, and its accompanying Sample Sheet "
             "(.tsv). Keep the folder structure the download provides - one "
             "subfolder per sample, each containing a .tsv counts file."),

            ("STEP 2 - Run DEG Analysis",
             "Open the 'Step 1 - DEG Analysis' page. Select:\n"
             "  1. The Gene Expression Quantification folder\n"
             "  2. The GDC Sample Sheet (.tsv)\n"
             "  3. An Output folder\n"
             "Set a project label (e.g. TCGA-COAD), thresholds "
             "(|log2FC| >= 1, padj < 0.05 are typical), then click "
             "'Run DEG Analysis'. Results are saved into Data/, Pictures/, Logs/."),

            ("STEP 3 - Load results into Visualisation",
             "When the analysis finishes, the cache file is auto-loaded "
             "into Step 2. Or open 'Step 2', click 'Load from Data folder' "
             "and pick the Data/ folder produced in Step 2."),

            ("STEP 4 - Choose a plot type",
             "'All Plots' (default) generates Volcano + MA + Bar + Heatmap. "
             "You can also pick a single plot type. The heatmap only "
             "renders as a real expression heatmap if Significant DEGs + "
             "normalized counts + sample grouping files are provided."),

            ("STEP 5 - Label genes (optional)",
             "Paste gene symbols in the 'Genes to label' box (one per line, "
             "or comma/semicolon separated). Both labeled and unlabeled "
             "copies are generated by default; toggle the option in the "
             "Output section if you only want labeled ones."),

            ("STEP 6 - Preview & Save",
             "Click 'Generate & Preview' - use the arrows / thumbnails to "
             "flip through the plots. When you are happy, click 'Save all "
             "previewed plots'. All figures are exported at your chosen "
             "DPI in the selected formats (PNG / PDF / TIFF / SVG). A "
             "combined side-by-side figure (Volcano + MA + Bar) is also "
             "saved when the option is enabled."),

            ("STEP 7 - Cite this work",
             "If you use this tool in a publication or preprint, please "
             "cite it (see the 'Please cite this work' card above) and "
             "star the GitHub repository. Your acknowledgement supports "
             "continued development of free scientific software."),
        ]
        for title, text in tutorial_steps:
            step = ctk.CTkFrame(tut.body, fg_color=COL_INPUT,
                                corner_radius=10,
                                border_width=1, border_color=COL_BORDER)
            step.pack(fill="x", pady=6)
            ctk.CTkLabel(step, text=title,
                         font=ctk.CTkFont(size=13, weight="bold"),
                         text_color=COL_ACCENT, anchor="w"
                         ).pack(fill="x", padx=14, pady=(10, 4))
            ctk.CTkLabel(step, text=text,
                         font=ctk.CTkFont(size=11),
                         text_color=COL_TEXT, anchor="w",
                         justify="left", wraplength=780,
                         ).pack(fill="x", padx=14, pady=(0, 12))

        # ---- Acknowledgement footer -------------------------------------
        foot = ctk.CTkFrame(scroll, fg_color="transparent")
        foot.pack(fill="x", pady=(6, 14))
        ctk.CTkLabel(
            foot,
            text=("This tool was built with significant effort - "
                  "if it helped your research, please cite and share it. Thank you!"),
            font=ctk.CTkFont(size=11, slant="italic"),
            text_color=COL_TEXT_SEC,
            justify="center", wraplength=780,
        ).pack()


# ===========================================================================
# 14)  Entry point
# ===========================================================================
def main():
    multiprocessing.freeze_support()
    app = DEGPipelineApp()
    app.mainloop()


if __name__ == "__main__":
    main()
