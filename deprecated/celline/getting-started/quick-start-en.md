# Quick Start Tutorial

This tutorial will teach you a complete workflow for performing single-cell RNA-seq analysis using Celline with real data.

## 🎯 What You'll Learn in This Tutorial

- Creating and configuring a Celline project
- Obtaining samples from public databases
- Data preprocessing and quality control
- Interactive analysis and visualization

## 📊 Data Used

This tutorial uses single-cell RNA-seq data from mouse brain cells (GSE115189).

- **Dataset**: GSE115189
- **Number of samples**: 2 samples
- **Cell count**: Approximately 3,000 cells
- **Species**: Mus musculus

## 🚀 Step 1: Project Setup

### Creating a Project Directory

```bash
# Create working directory
mkdir celline-tutorial
cd celline-tutorial

# Initialize Celline project
celline init tutorial-project
```

### Configuration Verification

```bash
# Verify installation is correctly completed
celline list

# Display system information
celline info

# Perform basic configuration
celline config
```

In the configuration screen, select the following:
- **Execution system**: `multithreading` (local execution)
- **Number of threads**: `2` (adjust according to your system)

## 🔍 Step 2: Data Acquisition

### Adding Samples

```bash
# Add sample set from GEO
celline run add GSE115189

# To add individual samples
# celline run add GSM3169075 GSM3169076
```

### Verifying Sample Information

```bash
# Check added samples
cat samples.toml
```

Example output:
```toml
GSM3169075 = "Control sample 1"
GSM3169076 = "Control sample 2"
```

## ⬇️ Step 3: Data Download

### Downloading Sequencing Data

```bash
# Parallel download (2 threads)
celline run download --nthread 2

# Check progress
ls -la resources/
```

### Verifying Data Structure

```bash
# Check downloaded files
find resources/ -name "*.fastq.gz" | head -5
```

## 🧮 Step 4: Count Processing

### Counting with Cell Ranger

```bash
# Configure transcriptome database
celline run set_transcriptome

# Execute count processing with Cell Ranger
celline run count
```

This process takes time (several hours~). Upon completion, the following will be generated for each sample:
- Feature-barcode matrix
- Quality metrics
- Cell/gene summary

## 🔬 Step 5: Quality Control and Preprocessing

### Creating Seurat Objects

```bash
# Create Seurat objects
celline run create_seurat
```

### Executing Preprocessing

```bash
# Execute quality control and preprocessing
celline run preprocess

# Check results
ls -la data/*/cell_info.tsv
```

### Checking QC Metrics

```bash
# Check QC statistics for each sample
head data/GSM3169075/cell_info.tsv
```

Example output:
```tsv
barcode	project	sample	cell	cell_type	include	n_genes_by_counts	pct_counts_mt	doublet_score	predicted_doublets
AAACCTGAGAAGGCCT-1	GSE115189	GSM3169075	GSM3169075_1	Neuron	true	2847	1.2	0.05	false
```

## 🌐 Step 6: Interactive Analysis

### Starting the Web Interface

```bash
# Start interactive mode
celline interactive
```

Access http://localhost:8080 in your browser.

### Interface Features

1. **Sample Overview**: Sample list and processing status
2. **Quality Control**: QC metrics visualization
3. **Cell Type Analysis**: Cell type distribution verification
4. **Gene Expression**: Gene expression pattern exploration

## 📈 Step 7: Verifying Analysis Results

### Checking Basic Statistics

```python
# Example verification with Python script
import polars as pl
import scanpy as sc

# Load cell_info.tsv
cell_info = pl.read_csv("data/GSM3169075/cell_info.tsv", separator="\t")

# Display basic statistics
print("Cell count:", len(cell_info))
print("Included cells:", cell_info.filter(pl.col("include")).height)
print("Cell types:", cell_info["cell_type"].unique().to_list())
```

### Visualizing Quality Metrics

```python
# Visualizing QC metrics
import matplotlib.pyplot as plt

# Gene count distribution
plt.figure(figsize=(10, 6))
plt.subplot(1, 2, 1)
plt.hist(cell_info["n_genes_by_counts"], bins=50)
plt.xlabel("Number of genes")
plt.ylabel("Number of cells")

# Mitochondrial gene percentage distribution
plt.subplot(1, 2, 2)
plt.hist(cell_info["pct_counts_mt"], bins=50)
plt.xlabel("Mitochondrial gene percentage")
plt.ylabel("Number of cells")

plt.tight_layout()
plt.show()
```

## 🔬 Step 8: Advanced Analysis (Optional)

### Cell Type Prediction

```bash
# Cell type prediction using pre-trained models
celline run predict_celltype
```

### Dimensionality Reduction and Clustering

```bash
# Dimensionality reduction (PCA, UMAP)
celline run reduce

# Batch effect correction
celline run integrate
```

### Custom Analysis

```python
# Custom analysis example in Python
from celline import Project
from celline.data import Seurat

# Load project
project = Project("./")

# Get Seurat object
seurat = project.seurat("GSE115189", "GSM3169075")

# Execute custom analysis
# ...
```

## 📊 Interpreting Results

### Quality Control Metrics

- **n_genes_by_counts**: Range of 200-5000 is typical
- **pct_counts_mt**: 5% or less is recommended
- **doublet_score**: 0.3 or less is typical

### Cell Filtering

Statistics after filtering:
```bash
# Check filtering statistics
celline run info
```

## 🔄 Workflow Automation

### Batch Processing Script

```bash
#!/bin/bash
# complete_workflow.sh

# Project initialization
celline init

# Add samples
celline run add GSE115189

# Data processing pipeline
celline run download --nthread 4
celline run count
celline run create_seurat
celline run preprocess

# Analysis
celline run predict_celltype
celline run reduce

echo "Analysis complete! Check results with celline interactive."
```

### Python API Usage Example

```python
from celline import Project
from celline.functions import Add, Download, Count, Preprocess

# Create project
project = Project("./tutorial-project")

# Analysis pipeline
project.call(Add([Add.SampleInfo(id="GSE115189")]))
project.call(Download())
project.call(Count())
project.call(Preprocess())

print("Analysis complete!")
```

## 🎉 Congratulations!

In this tutorial, you have completed the following:

- ✅ Creating a Celline project
- ✅ Data acquisition from public databases
- ✅ Data preprocessing and quality control
- ✅ Using the interactive analysis environment

## 📚 Next Steps

- [CLI Reference](/celline/cli) - More detailed command options
- [Interactive Mode](/celline/interactive) - Detailed web interface features
- [API Reference](/celline/api) - Usage with Python programming
- [Developer Guide](/celline/development) - Developing custom features

## 💡 Tips and Best Practices

### Memory Management

```bash
# For large datasets
celline config --nthread 1  # Reduce thread count
```

### Progress Monitoring

```bash
# Check progress with log files
tail -f resources/*/log/*.log
```

### Data Backup

```bash
# Backup important results
tar -czf results_backup.tar.gz data/ results/
```

---

::alert{type="success"}
Tutorial complete! For further analysis, please see [Advanced Usage](/celline/cli/advanced).
::