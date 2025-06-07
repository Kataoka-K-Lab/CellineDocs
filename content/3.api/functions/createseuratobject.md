# CreateSeuratObject Function

Create Seurat objects from Cell Ranger count matrices for downstream analysis.

## Overview

The `CreateSeuratObject` function creates Seurat objects from Cell Ranger output matrices. Seurat objects are the standard data structure for single-cell RNA-seq analysis in R, providing a unified framework for data storage and analysis.

## Class Information

- **Module**: `celline.functions.create_seurat`
- **Class**: `CreateSeuratObject`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `useqc_matrix` | `bool` | Yes | Whether to use quality-controlled filtered matrix |
| `then` | `Optional[Callable[[str], None]]` | No | Callback function executed on successful completion |
| `catch` | `Optional[Callable[[subprocess.CalledProcessError], None]]` | No | Callback function executed on error |

### JobContainer Structure

The `CreateSeuratObject.JobContainer` contains job configuration:

| Field | Type | Description |
|-------|------|-------------|
| `nthread` | `str` | Number of threads (fixed to 1) |
| `cluster_server` | `str` | Cluster server name (if applicable) |
| `jobname` | `str` | Job identifier |
| `logpath` | `str` | Path to log file |
| `r_path` | `str` | Path to R script directory |
| `exec_root` | `str` | Execution root directory |
| `input_h5_path` | `str` | Path to input HDF5 matrix file |
| `data_dir_path` | `str` | Output data directory |
| `proj_name` | `str` | Project name |
| `useqc_matrix` | `str` | QC matrix flag ("true" or "false") |

## Usage Examples

### Python API

#### Basic Usage

```python
from celline import Project
from celline.functions.create_seurat import CreateSeuratObject

# Create project
project = Project("./my-project")

# Create Seurat objects with filtered matrices
seurat_function = CreateSeuratObject(useqc_matrix=True)

# Execute function
result = project.call(seurat_function)
```

#### Using Raw Matrices

```python
from celline import Project
from celline.functions.create_seurat import CreateSeuratObject

# Create project
project = Project("./my-project")

# Create Seurat objects with raw matrices (no Cell Ranger filtering)
seurat_function = CreateSeuratObject(useqc_matrix=False)

# Execute function
result = project.call(seurat_function)
```

#### With Callbacks

```python
from celline import Project
from celline.functions.create_seurat import CreateSeuratObject
import subprocess

def on_success(sample_id: str):
    print(f"Successfully created Seurat object for: {sample_id}")

def on_error(error: subprocess.CalledProcessError):
    print(f"Seurat object creation failed: {error}")

# Create project
project = Project("./my-project")

# Create function with callbacks
seurat_function = CreateSeuratObject(
    useqc_matrix=True,
    then=on_success,
    catch=on_error
)

# Execute function
result = project.call(seurat_function)
```

### CLI Usage

#### Basic Usage

```bash
# Create Seurat objects (filtered matrices)
celline run createseuratobject

# Create with raw matrices
celline run createseuratobject --raw

# Verbose output
celline run createseuratobject --verbose
```

## Implementation Details

### Prerequisites

The function requires samples to be:
1. **Counted**: Cell Ranger count must be completed
2. **Cell Type Predicted**: Cell type prediction must be finished
3. **Preprocessed**: Quality control preprocessing must be done

### R Script Integration

The function executes R scripts to create Seurat objects:

```r
# Load required libraries
library(Seurat)
library(hdf5r)

# Read Cell Ranger output
matrix_path <- "/path/to/filtered_feature_bc_matrix.h5"
expression_matrix <- Read10X_h5(matrix_path)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = expression_matrix,
  project = "MyProject",
  min.cells = 3,
  min.features = 200
)

# Save Seurat object
saveRDS(seurat_obj, file = "seurat.rds")
```

### File Processing

The function processes Cell Ranger outputs:

```
Input:  resources/SAMPLE_ID/counted/outs/filtered_feature_bc_matrix.h5
Output: data/SAMPLE_ID/seurat.rds
```

### Directory Structure

Output organization:

```
project_root/
├── data/
│   └── SAMPLE_ID/
│       ├── seurat.rds           # Seurat object
│       ├── src/
│       │   └── create_seurat.sh # Generated R script
│       └── log/
│           └── create_seurat_*.log
```

## Matrix Selection

### Filtered vs Raw Matrices

| Matrix Type | Description | Use Case |
|-------------|-------------|----------|
| **Filtered** | Cell Ranger filtered cells and features | Standard analysis |
| **Raw** | All detected barcodes and features | Custom filtering workflows |

### Quality Control Impact

Using `useqc_matrix=True`:
- Applies Cell Ranger's cell calling algorithm
- Removes low-quality cells and features
- Reduces computational burden
- Recommended for most analyses

Using `useqc_matrix=False`:
- Preserves all detected barcodes
- Allows custom quality control
- Increases data size and processing time
- Useful for specialized analyses

## Seurat Object Structure

### Standard Components

Created Seurat objects contain:

```r
# Seurat object structure
seurat_obj@assays$RNA@counts        # Raw count matrix
seurat_obj@assays$RNA@data          # Normalized data (initially same as counts)
seurat_obj@meta.data               # Cell metadata
seurat_obj@reductions              # Dimensionality reductions (empty initially)
seurat_obj@graphs                  # Cell-cell graphs (empty initially)
```

### Metadata Integration

The function automatically adds metadata:

| Column | Description |
|--------|-------------|
| `orig.ident` | Sample identifier |
| `nCount_RNA` | Total UMI count per cell |
| `nFeature_RNA` | Number of detected genes per cell |
| `sample_id` | Original sample accession ID |

## Error Handling

### Prerequisite Checking

The function validates prerequisites:

```python
# Sample must be counted, predicted, and preprocessed
if not (sample.path.is_predicted_celltype and sample.path.is_preprocessed):
    print(f"Sample {sample.schema.key} is not ready. Skip")
    continue
```

### Existing Object Detection

```python
# Skip if Seurat object already exists
if os.path.isfile(f"{sample.path.data_sample}/seurat.rds"):
    print(f"Sample {sample.schema.key} already processed. Skip")
    continue
```

### Common Issues

1. **Missing R Dependencies**: Ensure Seurat and hdf5r packages are installed
2. **Insufficient Memory**: Large datasets require substantial RAM
3. **File Permissions**: Check read/write permissions
4. **Corrupted Matrices**: Validate Cell Ranger output integrity

## R Environment Requirements

### Required Packages

```r
# Install required R packages
install.packages(c("Seurat", "hdf5r", "Matrix"))

# For Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

### Memory Configuration

```r
# Increase memory limit if needed
memory.limit(size = 32000)  # 32 GB limit on Windows

# For large datasets
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB
```

## Performance Considerations

### Memory Requirements

| Dataset Size | Memory Needed |
|--------------|---------------|
| <5K cells | 8 GB |
| 5K-20K cells | 16 GB |
| 20K-50K cells | 32 GB |
| >50K cells | 64+ GB |

### Processing Time

Typical processing times:

| Dataset Size | Processing Time |
|--------------|----------------|
| <5K cells | 1-2 minutes |
| 5K-20K cells | 2-5 minutes |
| 20K-50K cells | 5-15 minutes |
| >50K cells | 15+ minutes |

## Cluster Computing

### Job Submission

For cluster environments:

```bash
#!/bin/bash
#PBS -N CreateSeurat
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb
#PBS -l walltime=2:00:00

module load R/4.1.0
Rscript create_seurat.R
```

### Resource Allocation

Recommended cluster resources:
- **CPU**: 1 core (R is primarily single-threaded for this task)
- **Memory**: 32-64 GB depending on dataset size
- **Walltime**: 2-4 hours for large datasets

## Methods

### `call(project: Project) -> Project`

Main execution method that creates Seurat objects for all eligible samples.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Iterates through all samples in the project
2. Checks prerequisites (counted, predicted, preprocessed)
3. Skips samples with existing Seurat objects
4. Generates R scripts for Seurat object creation
5. Executes scripts using ThreadObservable

## Output Validation

### Quality Checks

The function validates successful creation:

```r
# Validate Seurat object
if (class(seurat_obj) == "Seurat") {
  cat("Successfully created Seurat object\n")
  cat("Number of cells:", ncol(seurat_obj), "\n")
  cat("Number of features:", nrow(seurat_obj), "\n")
} else {
  stop("Failed to create valid Seurat object")
}
```

### File Verification

```python
# Check if output file exists and is valid
output_file = f"{sample.path.data_sample}/seurat.rds"
if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
    print("Seurat object created successfully")
else:
    raise FileNotFoundError("Failed to create Seurat object")
```

## Integration with Pipeline

### Typical Workflow

```python
from celline import Project
from celline.functions.count import Count
from celline.functions.preprocess import Preprocess
from celline.functions.predict_celltype import PredictCelltype
from celline.functions.create_seurat import CreateSeuratObject

# Complete pipeline
project = Project("./my-project")

# Process data
project.call(Count(nthread=8))
project.call(Preprocess())
project.call(PredictCelltype())

# Create Seurat objects
project.call(CreateSeuratObject(useqc_matrix=True))
```

## Related Functions

- [Count](count) - Generate count matrices before Seurat object creation
- [Preprocess](preprocess) - Quality control preprocessing
- [PredictCelltype](predict_celltype) - Cell type prediction
- [Integrate](integrate) - Data integration using Seurat objects

## Troubleshooting

### Common Issues

1. **R Package Missing**: Install required R packages
2. **Memory Error**: Increase system memory or use filtered matrices
3. **File Not Found**: Ensure Cell Ranger count completed successfully
4. **Permission Denied**: Check file system permissions

### Debug Mode

Enable detailed R logging:

```r
# Add to R script for debugging
options(error = traceback)
sessionInfo()
```

### Manual Execution

For debugging, run R script manually:

```bash
# Navigate to sample directory
cd data/SAMPLE_ID/src/

# Execute R script manually
Rscript create_seurat.sh
```