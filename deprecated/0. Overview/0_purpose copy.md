---
title: 0.1  Purpose
navigation: true      # ← 見出し専用にしてリンクは非表示
---

# 0.1  Purpose



---

## Table of Contents

1. [Overview](#0-overview)
2. [Quick‑Start Tutorial (`interfaces.py`)](#1-interfacespy-quick-start-tutorial)
3. [Calling Functions Programmatically](#2-calling-functions-programmatically)
4. [Function Reference](#3-function-reference)

   * [`download`](#download)
   * [`preprocess`](#preprocess)
5. [Developer Guide](#4-developer-guide)
6. [Appendix](#5-appendix)

---

## 0. Overview

### 0.1  Purpose

Celline streamlines the **entire lifecycle** of public scRNA‑seq/snRNA‑seq datasets: discovery, download, metadata harmonisation, quality‑control, preprocessing and format conversion (AnnData/Seurat) – all from a single Python entry point or CLI.

It is designed for *large‑scale* meta‑analyses where **hundreds of GEO/SRA accessions** must be treated identically and reproducibly.

### 0.2  Core Features

| Feature                   | What it does                                                                             | Why it matters                                     |
| ------------------------- | ---------------------------------------------------------------------------------------- | -------------------------------------------------- |
| Automated data retrieval  | Finds and downloads FASTQ/SRA + rich metadata via GEO, SRA, ENA & AWS mirrors            | Eliminates manual curation, guarantees provenance  |
| Uniform preprocessing     | Wraps **FastQC**, **STAR/Salmon**, **Cell Ranger** & **Seurat** with consistent defaults | Comparable downstream metrics across studies       |
| Batch‑aware integration   | Integrates data with **ContinuousVI**, **Scanorama** & **Harmony**                       | Removes batch effects while preserving biology     |
| Python‑CLI parity         | Every function is callable as `Project.call()` *or* `celline <subcommand>`               | Fits both scripting and pipeline‑automation styles |
| Plugin architecture       | New functions discovered at runtime via entry‑points                                     | Extensible without touching core code              |
| Re‑entrant & checkpointed | Caches each step under versioned directory tree                                          | Crash‑safe, lets you resume mid‑pipeline           |

### 0.3  High‑Level Architecture

```
┌────────────┐     user code/CLI      ┌────────────┐
│ interfaces │   ───────────────►    │ functions/ │
│  Project   │◄───────────────       │  plugins   │
└────────────┘      status           └────────────┘
      ▲                ▲                    │
      │                │                    ▼
      │           Project cache       plugin registry
      │
      └────────► CLI wrapper (celline)
```

---

## 1. `interfaces.py` Quick‑Start Tutorial

### 1.1  Installation

```bash [install.sh]
pip install celline   # latest release
# or bleeding‑edge
pip install git+https://github.com/Kataoka-K-Lab/Celline.git
```

*R dependencies* – Seurat conversion requires `R >= 4.2`, `Seurat`, `SeuratDisk`, `reticulate`.

### 1.2  Initialise a project directory – *3 lines to start*

```py [project.py]
from celline.interfaces import Project
proj = Project("/data/celline_work", proj_name="brain_aging")
proj.info()   # prints location, disk usage, config snapshot
```

Created structure:

```
/data/celline_work/
  raw/        # FASTQ / SRA
  processed/  # .h5ad / .rds
  meta/       # metadata TSV, JSON
  reports/    # HTML QC reports
  logs/       # logfiles per run
  tmp/        # scratch space
```

### 1.3  End‑to‑end pipeline (young mouse dataset)

```py
proj.run_pipeline([
    ("download",   {"accession": "GSE173917"}),
    ("preprocess", {"min_cells": 200, "min_genes": 500}),
    ("integrate",  {"method": "continuousvi"}),
    ("reduce",     {"n_pcs": 50})
])
```

Each tuple is *function name, kwargs*. Progress is logged to **logs/pipeline.log**.

### 1.4  Dive into results

```py
adata = proj.load_count_matrix(sample="GSE173917", format="anndata")
adata.obs.head()
```

### 1.5  CLI parity

```bash
celline download GSE173917 -w /data/celline_work && \
celline preprocess -w /data/celline_work --min-cells 200 --min-genes 500
```

### 1.6  Project object anatomy

| Attribute        | Type             | Purpose                                                 |
| ---------------- | ---------------- | ------------------------------------------------------- |
| `workdir`        | `Path`           | Root of project tree                                    |
| `config`         | `dict`           | Parsed from *settings.toml* – tweak threads, references |
| `logger`         | `logging.Logger` | File‑ & console‑level logs                              |
| `state`          | `dict`           | Last‑known status of each function per sample           |
| `call()`         | method           | Enqueue & run a function; returns `TaskHandle`          |
| `run_pipeline()` | method           | Sugar for sequential calls                              |

`TaskHandle` exposes `.status`, `.returncode`, `.stdout`.

### 1.7  FAQ (selected)

| Symptom                          | Cause                    | Remedy                                                 |
| -------------------------------- | ------------------------ | ------------------------------------------------------ |
| **`SRA accession not found`**    | GEO lacks SRR link       | `--ena` flag forces ENA route                          |
| **R error `cannot find Seurat`** | R libs missing           | Install R deps snippet                                 |
| **Disk full**                    | default tmp on small SSD | set `TMPDIR=/big_disk/tmp` or `proj.config["tmp_dir"]` |

---

## 2. Calling Functions Programmatically

### 2.1  Function categories

| Category                           | Examples                                         | Typical output                       |
| ---------------------------------- | ------------------------------------------------ | ------------------------------------ |
| **I/O**                            | `download`, `sync_DB`, `set_transcriptome`       | FASTQ, reference files               |
| **Pre‑processing**                 | `preprocess`, `count`, `vcount`, `create_seurat` | filtered matrices                    |
| **Integration / Batch correction** | `integrate`, `reduce`, `batch_cor`               | latent space, PCA/Harmony embeddings |
| **Annotation**                     | `predict_celltype`, `interactive`                | cell type labels, HTML dashboards    |
| **Meta utilities**                 | `info`, `bash`, `job`                            | workflow glue                        |

### 2.2  Signatures & handles

```py
handle = proj.call("preprocess", min_cells=500, wait=False)
print(handle.status)  # queued / running / error / done
handle.join(timeout=3600)  # block until finished
```

### 2.3  Dependency graph

Celline does **not** enforce DAGs automatically; you call in the order you need. Future versions will allow defining dependencies so that `integrate` waits for all `preprocess` tasks.

---

## 3. Function Reference

### <a name="download"></a> 3.1  `download`

*(No change – see earlier section)*

### <a name="preprocess"></a> 3.2  `preprocess`

| Field              | Description                                                          |
| ------------------ | -------------------------------------------------------------------- |
| **Purpose**        | QC and barcode filtering, then gene‑by‑cell matrix construction.     |
| **CLI**            | `celline preprocess`                                                 |
| **Inputs**         | Raw FASTQ in `raw/<sample>/`; parameters below                       |
| **Outputs**        | `.h5ad` and `.h5seurat` under `processed/`; QC HTML under `reports/` |
| **Key parameters** |                                                                      |

| Name        | Type               | Default        | Notes                                 |
| ----------- | ------------------ | -------------- | ------------------------------------- |
| `min_cells` | int                | 3              | filter genes expressed in fewer cells |
| `min_genes` | int                | 200            | filter cells with low gene counts     |
| `chemistry` | str\|None          | auto           | 10x chemistry hint (v2, v3)           |
| `aligner`   | {"star", "salmon"} | star           | transcriptomic aligner                |
| `threads`   | int                | config‑default | CPU threads                           |

**Processing flow**

1. `fastqc` each FASTQ – summary to *reports/*.
2. Align reads (`STAR` default) & count with `featureCounts`.
3. Load matrix into **Scanpy**, apply gene/cell filters.
4. Save as `.h5ad`; optional conversion to `.h5seurat` via `SeuratDisk`.

Error handling mirrors the `download` layout.

---

## 4. Developer Guide

> *Goal: enable you to add a new capability in **< 1 hour** and ship it to users.*

### 4.1  Function taxonomy

| Use‑case                                  | Recommended base                     | I/O pattern                         |
| ----------------------------------------- | ------------------------------------ | ----------------------------------- |
| **Run external command** (e.g. STAR)      | `ShellFunction(_base.ShellFunction)` | produce files, parse logs           |
| **Pure Python analysis** (Scanpy metrics) | `PythonFunction`                     | return Python objects + plots       |
| **R‑centric step** (Seurat clustering)    | `RFunction`                          | wrap `subprocess` call to `Rscript` |

### 4.2  File layout & naming

```
celline/
  functions/
    myfunc.py          # new function class
    __init__.py        # add `from .myfunc import myfunc`
  templates/
    r/
      myfunc.R         # if needed
  tests/
    test_myfunc.py
```

Name your class **snake\_case → CamelCase**: *`count_genes` → `CountGenes`*.
CLI subcommand auto‑generates from file name (`celline count-genes`).

### 4.3  Step‑by‑step: writing `qc_plot`

```py
# celline/functions/qc_plot.py
from celline.functions._base import PythonFunction
import scanpy as sc, seaborn as sns

class qc_plot(PythonFunction):
    """Generate PDF of QC violin plots"""
    def __init__(self, sample: str):
        super().__init__(name="qc_plot")
        self.sample = sample
    def run(self, project):
        ad = project.load_count_matrix(self.sample)
        sc.pl.violin(ad, ["n_genes", "percent_mito"], show=False, save="qc.pdf")
        self.out_files.append(project.workdir/"reports"/f"{self.sample}_qc.pdf")
```

Add to `functions/__init__.py`:

```py
from .qc_plot import qc_plot
```

Write a minimal test:

```py
# tests/test_qc_plot.py
from celline.interfaces import Project
proj = Project("tmp")
proj.state["processed"]["SAMPLE"] = "dummy.h5ad"  # fixture
proj.call("qc_plot", sample="SAMPLE")
assert (proj.workdir/"reports"/"SAMPLE_qc.pdf").exists()
```

### 4.4  Advanced patterns

* **Chained functions** – define `requires = ["download"]` so scheduler auto‑orders steps (coming soon).
* **Parallel per‑sample** – override `split_samples()` to return iterable of tasks.
* **Streaming output** – provide generator from `command()` to yield incremental logs.

### 4.5  Suggested roadmap functions

| Name               | Category    | Short description                            |
| ------------------ | ----------- | -------------------------------------------- |
| `doublet_detect`   | QC          | Wrap **Scrublet** to flag doublets           |
| `gene_set_score`   | Analysis    | Compute module scores for hallmark gene sets |
| `cellchat_export`  | Interaction | Export expression matrices for **CellChat**  |
| `differential_age` | Statistics  | Pseudo‑bulk DE test across age groups        |
| `report_html`      | Reporting   | Collate outputs into browsable dashboard     |

### 4.6  Testing & CI tips

* Use **pytest‑fixtures** with tiny 1k‑cell `.h5ad` stub.
* Mock network calls (`pytest‑requests‑mock`) when unit‑testing `download` logic.
* GitHub Actions provides `sra-tools`, `STAR`, `R` via `setup‑conda` matrix.

### 4.7  Contributing checklist

1. `black . && isort .`
2. `pytest -q`
3. Update *Function Reference* in this doc.
4. Add yourself to `AUTHORS.md`.

---

## 5. Appendix

### 5.1  `settings.toml` keys (excerpt)

| Key                | Type | Default    | Explanation                            |
| ------------------ | ---- | ---------- | -------------------------------------- |
| `threads`          | int  | phys‑cores | Global CPU pool size                   |
| `tmp_dir`          | str  | TMPDIR     | scratch space; heavy I/O use           |
| `reference_genome` | str  | "mm10"     | alias in `celline/resources/refs.toml` |
| `rscript_path`     | str  | auto       | custom `Rscript` if not on PATH        |

### 5.2  R Tips

```py
import celline.utils.r as cr
seurat_obj = cr.read_h5seurat("processed/SAMPLE.h5seurat")
cr.run("DimPlot", seurat_obj, group_by="celltype")
```

### 5.3  Known Issues

* `fasterq-dump` RAM spike – set `--mem` or switch to HTTP mirror.
* Windows not officially supported (works via WSL2).

### 5.4  Version Matrix

| Celline | Python   | R    | Seurat |
| ------- | -------- | ---- | ------ |
| v1.4.x  | 3.9–3.12 | ≥4.2 | 5.0    |

---

*End of document.*
