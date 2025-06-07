# Installation

This document provides detailed instructions for installing Celline.
Please choose the optimal installation method according to your environment.

## 🐍 Python Environment

### System Requirements

| Requirement | Minimum Version | Recommended Version |
|-------------|-----------------|-------------------|
| Python | 3.10 | 3.11+ |
| pip | 21.0 | Latest |
| Memory | 4GB | 8GB+ |
| Storage | 1GB | 10GB+ |

### Creating Python Virtual Environment

```bash
# Using venv
python -m venv celline-env
source celline-env/bin/activate  # Linux/Mac
# or
celline-env\Scripts\activate     # Windows

# Using conda
conda create -n celline python=3.11
conda activate celline
```

## 📦 Installation Methods

### Method 1: Via pip (Recommended)

```bash
# Install latest stable version
pip install celline

# Upgrade
pip install --upgrade celline

# Specify a particular version
pip install celline==0.1.10
```

### Method 2: Via UV

```bash
# If UV is not installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Add Celline to project
uv add celline

# Or use temporarily
uv run --with celline celline --help
```

### Method 3: Development Version Installation

```bash
# Install latest development version from GitHub
pip install git+https://github.com/YUYA556223/celline.git

# Specify a particular branch
pip install git+https://github.com/YUYA556223/celline.git@develop
```

## 🔧 Dependencies

Celline depends on the following major packages:

### Python Dependencies

```text
argparse>=1.4.0
continuousvi>=0.1.5
fastapi>=0.104.0
inquirer>=3.4.0
pandas>=2.2.3
polars>=1.26.0
pyarrow>=19.0.1
pydantic>=2.5.0
requests>=2.31.0
rich>=14.0.0
scanpy>=1.11.1
scrublet>=0.2.3
tqdm>=4.67.1
uvicorn>=0.24.0
```

### R Dependencies

Celline uses R for some functionality. The following R packages are required:

```r
# Required packages
install.packages(c(
  "Seurat",
  "SeuratDisk", 
  "tidyverse",
  "hdf5r"
))

# Optional packages
install.packages(c(
  "scater",
  "scran",
  "SingleCellExperiment",
  "BiocManager"
))
```

## 🚀 Post-Installation Verification

### Basic Operation Check

```bash
# Check if Celline is installed
celline --help

# Check version
python -c "import celline; print('Celline version:', celline.__version__)"

# Display available functions
celline list
```

### System Information Check

```bash
# Display system information
celline info

# Check configuration status
celline config
```

## 🔧 External Tool Installation

### Cell Ranger (Optional)

Required for processing 10x Genomics data:

```bash
# Download Cell Ranger (account registration required)
# https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

# Add to PATH after installation
export PATH=/path/to/cellranger:$PATH
```

### R Environment Setup

```bash
# Check if R is installed
R --version

# Check R path
which R

# Set R path in Celline
celline config
# or
export R_HOME=/usr/lib/R
```

## 🐳 Docker Environment Usage

### Docker Image (Planned for Future Implementation)

```bash
# Docker image usage example (under development)
docker run -it celline/celline:latest

# Mount to use local data
docker run -v $(pwd):/workspace celline/celline:latest
```

## 🌐 Development Environment Setup

### Installation from Source

```bash
# Clone repository
git clone https://github.com/YUYA556223/celline.git
cd celline

# Install including development dependencies
pip install -e ".[dev]"

# Or using uv
uv sync --all-extras
```

### Development Tools Setup

```bash
# Code formatter
pip install ruff black

# Testing tools
pip install pytest pytest-cov

# Documentation generation
pip install sphinx sphinx-rtd-theme
```

## 🔍 Troubleshooting

### Common Issues and Solutions

#### 1. Python Version Related

```bash
# If Python version is too old
python --version  # Confirm it's 3.10+

# Use pyenv to manage Python versions
pyenv install 3.11.0
pyenv global 3.11.0
```

#### 2. Dependency Conflicts

```bash
# Resolve dependency conflicts
pip install --force-reinstall celline

# Or reinstall in clean environment
pip uninstall celline
pip install celline
```

#### 3. R-related Issues

```bash
# If R is not found
sudo apt-get install r-base r-base-dev  # Ubuntu/Debian
brew install r                          # macOS

# If R package installation fails
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

#### 4. Memory-related Issues

```bash
# Reduce memory usage
celline config --nthread 1

# Set up swap file (Linux)
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

#### 5. Permission-related Issues

```bash
# Install at user level
pip install --user celline

# Avoid permission issues
pip install --no-deps celline
```

## ✅ Installation Verification Checklist

- [ ] Python 3.10+ is installed
- [ ] pip is the latest version
- [ ] `celline --help` executes normally
- [ ] `celline list` displays function list
- [ ] R environment is configured (if using R)
- [ ] Required R packages are installed
- [ ] Sufficient disk space is available

## 🔄 Upgrades

### Regular Upgrades

```bash
# Upgrade to latest version
pip install --upgrade celline

# Downgrade to specific version
pip install celline==0.1.9

# Upgrade to development version
pip install --upgrade git+https://github.com/YUYA556223/celline.git
```

### Configuration File Migration

New versions may change the configuration file format:

```bash
# Back up configuration
cp setting.toml setting.toml.backup

# Migrate to new configuration format (if needed)
celline config --migrate
```

---

> **Note**: If you encounter issues during installation, please refer to the [Troubleshooting](/celline/troubleshooting) section.