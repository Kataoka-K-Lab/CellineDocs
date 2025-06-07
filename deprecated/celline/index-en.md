# Celline - Single Cell RNA-seq Analysis Pipeline

Celline is a comprehensive, interactive pipeline for single-cell RNA sequencing (scRNA-seq) analysis designed to streamline the workflow from raw FASTQ files to biological insights.

## 🧬 Key Features

### Data Processing
- **Automated Data Processing**: From raw FASTQ files to expression matrices
- **Quality Control**: Built-in QC metrics and filtering
- **Public Database Integration**: Automatic data download from SRA and GEO

### Analysis Capabilities
- **Dimensionality Reduction**: PCA, t-SNE, and UMAP implementations
- **Clustering**: Multiple clustering algorithms
- **Cell Type Prediction**: Automated cell type annotation
- **Batch Effect Correction**: Multiple methods for data integration

### User Interface
- **Interactive Visualization**: Web-based interface for data exploration
- **CLI**: Rich command-line tools
- **Reproducible Workflows**: Containerized environments and version control

## 🚀 Quick Start

### Installation

```bash
pip install celline
```

### Basic Usage

```bash
# Initialize a new project
celline init

# List available functions
celline list

# Run preprocessing
celline run preprocess

# Launch interactive interface
celline interactive
```

### Usage with UV Environment

```bash
# Execution via UV
uv run celline run add GSE123456
uv run celline run download
uv run celline interactive
```

## 📊 Workflow Overview

1. **Project Initialization** - `celline init`
2. **Sample Addition** - `celline run add [accession_id]`
3. **Data Download** - `celline run download`
4. **Count Processing** - `celline run count`
5. **Preprocessing & QC** - `celline run preprocess`
6. **Analysis & Visualization** - `celline interactive`

## 🏗️ Architecture

Celline consists of the following main components:

- **Project System**: Project management and configuration
- **Function Framework**: Extensible analysis function system
- **Database Layer**: Metadata and sample information management
- **Web Interface**: Interactive UI with FastAPI + Vue.js
- **CLI System**: Command-line interface with Rich UI

## 📚 Documentation Structure

### Getting Started
Installation, configuration, and quick start guide

### CLI Reference
Command-line usage and configuration

### API Reference
Detailed reference for functions and classes

### Interactive Mode
Web interface and interactive features

### Architecture
System design and architecture details

### Developer Guide
Creating custom functions and extensions

## 🔧 System Requirements

- **Python**: 3.10+
- **R**: 4.0+ (Seurat package)
- **Cell Ranger**: 7.0+ (optional)
- **Memory**: 8GB+ recommended
- **Storage**: Depends on project size

## 📖 Next Steps

- [Getting Started](/celline/getting-started) - Detailed setup guide
- [CLI Reference](/celline/cli) - Command-line usage
- [Interactive Mode](/celline/interactive) - Web interface
- [Architecture](/celline/architecture) - System design details

---

> **Note**: Celline is an open-source project developed for research purposes.
> For commercial use and support, please check the project repository.