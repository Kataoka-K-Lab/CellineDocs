# Getting Started

This is a complete guide to getting started with single-cell RNA-seq analysis using Celline.
This section provides step-by-step instructions from installation to your first analysis.

## 📋 Prerequisites

Before using Celline, please ensure you meet the following requirements:

### System Requirements

- **Python**: 3.10 or higher
- **R**: 4.0 or higher (Seurat package required)
- **Memory**: 8GB or more recommended
- **Storage**: Sufficient capacity according to project size

### Recommended Tools (Optional)

- **Cell Ranger**: For 10x Genomics data processing (7.0 or higher)
- **Docker**: For reproducible environment setup
- **Git**: For version control

## 🚀 Installation

### Installation via pip

```bash
# Install the latest version
pip install celline

# Install a specific version
pip install celline==0.1.10
```

### Installation in UV Environment

```bash
# Execution in UV environment
uv add celline

# Or execute directly via UV
uv run --with celline celline --help
```

### Development Version Installation

```bash
# Install development version from GitHub
pip install git+https://github.com/YUYA556223/celline.git
```

## ⚙️ Initial Setup

### 1. Project Initialization

```bash
# Start a project in a new directory
mkdir my-scrna-project
cd my-scrna-project

# Initialize Celline project
celline init
```

### 2. Configuration Check

```bash
# Check available functions
celline list

# Check system information
celline info

# Check and modify settings
celline config
```

## 📁 Project Structure

After initialization, the following directory structure is created:

```
my-scrna-project/
├── setting.toml         # Project configuration
├── samples.toml         # Sample information
├── data/                # Analyzed data
├── resources/           # Raw data and metadata
├── results/             # Result files
└── integration/         # Integrated analysis results
```

## 🔧 Configuration Files

### setting.toml

Manages the overall project configuration:

```toml
[project]
name = "my-scrna-project"
version = "1.0.0"
description = "Single cell RNA-seq analysis project"

[execution]
system = "multithreading"  # or "PBS"
nthread = 4
pbs_server = ""

[R]
r_path = "/usr/bin/R"

[fetch]
wait_time = 4
```

### samples.toml

Manages information about samples to be analyzed:

```toml
# Samples are automatically added by the celline run add command
GSM1234567 = "Sample 1 description"
GSM1234568 = "Sample 2 description"
```

## 🎯 Basic Workflow

### Step 1: Add Samples

```bash
# Add samples from GEO
celline run add GSE123456

# Add multiple samples at once
celline run add GSM1234567 GSM1234568

# Add samples from file
celline run add --from-file samples.txt
```

### Step 2: Data Download

```bash
# Download data for added samples
celline run download

# Parallel download (4 threads)
celline run download --nthread 4
```

### Step 3: Count Processing

```bash
# Count processing with Cell Ranger
celline run count
```

### Step 4: Preprocessing

```bash
# Quality control and preprocessing
celline run preprocess

# Target specific cell types only
celline run preprocess --target-celltype Neuron Astrocyte
```

### Step 5: Interactive Analysis

```bash
# Launch web interface
celline interactive
```

## 🌐 Execution Modes

### CLI Execution

```bash
# Standard CLI execution
celline run [function_name] [arguments]

# Display help
celline help [function_name]
```

### UV Execution

```bash
# Execution in UV environment
uv run celline run [function_name] [arguments]

# Execute with project-specific dependencies
uv run --with celline --with scanpy celline run preprocess
```

### Python Execution

```python
from celline import Project
from celline.functions.add import Add
from celline.functions.download import Download

# Create project
project = Project("./my-project")

# Add samples
add_samples = Add([Add.SampleInfo(id="GSM1234567", title="Sample 1")])
project.call(add_samples)

# Download data
download = Download()
project.call(download)
```

## ⚡ Quick Start Example

Here's an example of a complete workflow:

```bash
# 1. Create project
mkdir test-project && cd test-project
celline init

# 2. Add samples (example: GSE115189)
celline run add GSE115189

# 3. Download data
celline run download --nthread 2

# 4. Count processing
celline run count

# 5. Preprocessing
celline run preprocess

# 6. Start interactive analysis
celline interactive
```

## 🔍 Verification and Troubleshooting

### Installation Verification

```bash
# Verify Celline is properly installed
celline --help

# Check version
python -c "import celline; print(celline.__version__)"
```

### Common Issues

1. **R packages not found**
   ```bash
   # Set R path
   celline config --system multithreading
   ```

2. **Parallel processing configuration**
   ```bash
   # Adjust thread count
   celline config --nthread 4
   ```

3. **Memory shortage**
   - Change to fewer threads
   - Run with smaller sample sizes

## 📚 Next Steps

- [Installation](/celline/getting-started/installation) - Detailed installation guide
- [Quick Start](/celline/getting-started/quick-start) - Complete tutorial
- [CLI Reference](/celline/cli) - Command-line details
- [Configuration](/celline/cli/configuration) - Detailed configuration

---

> **Tip**: Once setup is complete, check the available commands in [CLI Reference](/celline/cli).