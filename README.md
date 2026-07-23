# 🧬 DEG Pipeline & Visualizer: Multi-Cancer Analysis

![Python Version](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Stable-brightgreen)
![Version](https://img.shields.io/badge/Version-2.0.2-orange)

**Author:** Alireza Balaei Kahnamoei  
**GitHub:** [alirezabk1382927-sys](https://github.com/alirezabk1382927-sys)  
**LinkedIn:** [Alireza Balaei Kahnamoei](https://ir.linkedin.com/in/alireza-balaei-kahnamoei-aa8216344)

---

## 📖 The Problem We Are Solving

Let’s be honest—analyzing high-throughput RNA-seq data, especially from large consortia like **TCGA**, is not exactly a walk in the park. 

Traditionally, researchers have to wrestle with **R/Bioconductor**, complex command-line tools, and endless dependency conflicts. If you are a wet-lab biologist, a clinical researcher, or just someone who wants to focus on the *biology* rather than the *syntax*, this process can feel like an insurmountable wall. 

The learning curve is steep, the error messages are cryptic, and one small mistake in the code can cost you days of troubleshooting. 

---

## 🚀 The Solution: A Click-and-Go Pipeline

I built **DEG Pipeline & Visualizer** to bridge that gap. 

This tool takes everything that makes differential expression analysis complicated and wraps it into a **clean, point-and-click Graphical User Interface (GUI)**. You don’t need to write a single line of R or Python code. 

Just point the software to your raw count files, click "Run," and let it handle the heavy lifting using **PyDESeq2**. In return, you get publication-ready statistical results and high-resolution figures (Volcano plots, MA plots, and Heatmaps) in seconds. 

**And the best part?** It is completely **Open Source** and free for everyone to use, modify, and share.

---

## 🖥️ Two Ways to Use This Tool

### Option 1: For Non-Coders (The Standalone .exe)
If you don't want to install Python or deal with terminals, a **pre-compiled Windows executable (.exe)** is available for download. 
> *(Note: Check the "Releases" section of this repository to download the latest `.exe` file. Just download, double-click, and you are ready to go!)*

### Option 2: For Developers / Python Users (Run from Source)
If you prefer to run the raw Python script, you will need to install a few dependencies. Don't worry—I've made it as simple as possible.

1. **Clone the repository** (or download the ZIP):
   ```bash
   git clone https://github.com/alirezabk1382927-sys/DEG-Pipeline-Visualizer.git
   cd DEG-Pipeline-Visualizer
