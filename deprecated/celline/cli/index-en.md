# CLI Reference

This is the complete reference for Celline's command-line interface (CLI).
It provides detailed explanations of all commands, options, and usage examples.

## 🎯 Basic Syntax

```bash
celline [command] [options] [arguments]
```

### Global Options

| Option | Description |
|--------|-------------|
| `--help, -h` | Display help |
| `--version` | Display version information |

## 📋 Main Commands

### 🚀 Project Management

#### `init` - Project Initialization

```bash
celline init [project_name]
```

Initialize a new Celline project.

**Parameters:**
- `project_name` (optional): Project name

**Examples:**
```bash
# Initialize project with current directory name
celline init

# Initialize project with specified name
celline init my-scrna-project
```

#### `info` - System Information Display

```bash
celline info
```

Display Celline system information and available functions.

**Examples:**
```bash
celline info
```

### ⚙️ Configuration Management

#### `config` - Configuration Management

```bash
celline config [options]
```

Manage execution environment configuration.

**Options:**
- `--system {multithreading,PBS}`: Set execution system
- `--nthread N`: Set number of threads
- `--pbs-server SERVER`: Set PBS server name

**Examples:**
```bash
# Interactive configuration
celline config

# Configure multithreaded execution
celline config --system multithreading --nthread 4

# Configure PBS cluster execution
celline config --system PBS --pbs-server my-cluster
```

### 📚 Help System

#### `list` - Function List Display

```bash
celline list
```

Display all available functions in table format.

#### `help` - Detailed Help

```bash
celline help [function_name]
```

**Parameters:**
- `function_name` (optional): Function name for detailed help display

**Examples:**
```bash
# Display general help
celline help

# Display help for specific functions
celline help add
celline help preprocess
```

### 🔧 Function Execution

#### `run` - Function Execution

```bash
celline run <function_name> [function_args] [options]
```

Execute the specified function.

**Options:**
- `--project-dir, -p DIR`: Project directory (default: .)
- `--project-name, -n NAME`: Project name (default: default)

**Examples:**
```bash
# Basic function execution
celline run preprocess

# Function execution with arguments
celline run add GSE123456

# Execute with specified project
celline run preprocess --project-dir /path/to/project
```

### 🌐 Interactive Mode

#### `interactive` - Web Interface Launch

```bash
celline interactive
```

Launch web-based interactive interface and API server.

#### `api` - Launch API Server Only

```bash
celline api
```

Launch only the API server for testing purposes.

## 🔬 Analysis Functions

### 📊 Data Management

#### `add` - Sample Addition

```bash
celline run add <sample_ids> [options]
```

Add new samples to the project.

**Parameters:**
- `sample_ids`: Sample IDs to add (GSE, GSM, SRR, etc.)

**Options:**
- `--title, -t TITLE`: Sample title
- `--from-file, -f FILE`: Load sample IDs from file

**Examples:**
```bash
# Add single sample
celline run add GSE123456

# Add multiple samples
celline run add GSM789012 GSM789013 --title "My samples"

# Add from file
celline run add --from-file samples.txt

# Add from CSV file (requires id, title columns)
celline run add --from-file samples.csv
```

#### `download` - Data Download

```bash
celline run download [options]
```

Download sequencing data for added samples.

**Options:**
- `--nthread, -n N`: Number of download threads (default: 1)
- `--force, -f`: Force re-download of existing files

**Examples:**
```bash
# Basic download
celline run download

# Parallel download
celline run download --nthread 4

# Force re-download
celline run download --force
```

### 🧮 Data Processing

#### `count` - Count Processing

```bash
celline run count [options]
```

Generate expression counts from FASTQ files using Cell Ranger.

**Examples:**
```bash
celline run count
```

#### `create_seurat` - Seurat Object Creation

```bash
celline run create_seurat [options]
```

Create Seurat objects from count data.

**Examples:**
```bash
celline run create_seurat
```

### 🔬 Quality Control & Preprocessing

#### `preprocess` - Preprocessing

```bash
celline run preprocess [options]
```

Execute quality control, filtering, and normalization.

**Options:**
- `--target-celltype, -t TYPES`: Specify target cell types (multiple allowed)

**Examples:**
```bash
# Basic preprocessing
celline run preprocess

# Target specific cell types only
celline run preprocess --target-celltype Neuron Astrocyte
```

### 🔍 Analysis

#### `reduce` - Dimensionality Reduction

```bash
celline run reduce [options]
```

Execute dimensionality reduction using PCA, UMAP, t-SNE.

#### `predict_celltype` - Cell Type Prediction

```bash
celline run predict_celltype [options]
```

Predict cell types using machine learning models.

#### `integrate` - Data Integration

```bash
celline run integrate [options]
```

Execute data integration and batch effect correction for multiple samples.

### 🗄️ Database Management

#### `sync_DB` - Database Synchronization

```bash
celline run sync_DB [options]
```

Synchronize local database with the latest public databases.

#### `set_transcriptome` - Transcriptome Information Setting

```bash
celline run set_transcriptome [options]
```

Configure transcriptome information database for analysis.

## 🎛️ Execution Environment Configuration

### Multithreaded Execution

```bash
# Execute in parallel with 4 threads
celline config --system multithreading --nthread 4

# Function execution
celline run download --nthread 4
```

### PBS Cluster Execution

```bash
# Configure PBS system
celline config --system PBS --pbs-server my-cluster

# Job submission
celline run count  # Executed as PBS job
```

## 📝 Configuration Files

### setting.toml

```toml
[project]
name = "my-project"
version = "1.0.0"

[execution]
system = "multithreading"  # or "PBS"
nthread = 4
pbs_server = "my-cluster"

[R]
r_path = "/usr/bin/R"

[fetch]
wait_time = 4
```

### samples.toml

```toml
# Sample information (auto-generated)
GSM1234567 = "Sample 1 description"
GSM1234568 = "Sample 2 description"
```

## 🔄 Workflow Examples

### Basic Workflow

```bash
#!/bin/bash
# Complete analysis workflow

# 1. Project initialization
celline init my-analysis

# 2. Configuration
celline config --system multithreading --nthread 4

# 3. Sample addition
celline run add GSE123456

# 4. Data processing
celline run download --nthread 4
celline run count
celline run create_seurat

# 5. Quality control & preprocessing
celline run preprocess

# 6. Analysis
celline run predict_celltype
celline run reduce
celline run integrate

# 7. Interactive analysis
celline interactive
```

### Batch Processing

```bash
# Bulk processing of multiple samples
echo "GSE123456\nGSE789012\nGSE345678" > datasets.txt

while read dataset; do
    echo "Processing $dataset"
    celline run add $dataset
    celline run download
    celline run count
    celline run preprocess
done < datasets.txt
```

## 🚨 Error Handling

### Common Errors

```bash
# When function is not found
celline run unknown_function
# Error: Function 'unknown_function' not found.

# When required arguments are missing
celline run add
# Error: Function name is required.

# When configuration error occurs
celline run download
# Error: No samples found. Use 'celline run add' first.
```

### Log Checking

```bash
# Check error logs
cat resources/*/log/*.log

# Real-time log monitoring
tail -f resources/*/log/*.log
```

## 🔧 Debugging and Troubleshooting

### Getting Detailed Information

```bash
# Display detailed system information
celline info

# Check configuration status
celline config

# Check available functions
celline list
```

### Checking Execution Environment

```bash
# Check Python environment
python -c "import celline; print(celline.__version__)"

# Check R environment
which R
R --version
```

---

::alert{type="info"}
For more detailed usage information, please refer to [Configuration](/celline/cli/configuration) and [Advanced Usage](/celline/cli/advanced).
::