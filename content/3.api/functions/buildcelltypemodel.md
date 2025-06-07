# BuildCellTypeModel Function

Build custom cell type prediction models from reference single-cell datasets.

## Overview

The `BuildCellTypeModel` function creates custom scPred-based cell type classification models from user-provided reference datasets. These models can then be used for automated cell type prediction in new datasets with similar biological contexts.

## Class Information

- **Module**: `celline.functions.predict_celltype`
- **Class**: `BuildCellTypeModel`  
- **Base Class**: `CellineFunction`

## Related Classes

### `CellTypeModel`

Data structure for model configuration:

```python
@dataclass
class CellTypeModel:
    species: str
    suffix: Optional[str]
```

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `species` | `str` | Yes | Species name (e.g., "Homo sapiens", "Mus musculus") |
| `suffix` | `str` | Yes | Model identifier/version suffix |
| `nthread` | `int` | Yes | Number of threads for model training |
| `h5matrix_path` | `str` | Yes | Path to reference HDF5 matrix file |
| `celltype_path` | `str` | Yes | Path to cell type annotation TSV file |

### JobContainer Structure

| Field | Type | Description |
|-------|------|-------------|
| `nthread` | `str` | Number of threads for training |
| `cluster_server` | `str` | Cluster server name |
| `jobname` | `str` | Job identifier |
| `logpath` | `str` | Training log file path |
| `h5matrix_path` | `str` | Input matrix file path |
| `celltype_path` | `str` | Cell type annotations path |
| `dist_dir` | `str` | Output directory |
| `r_path` | `str` | R script directory |
| `exec_root` | `str` | Execution root |

## Usage Examples

### Python API

#### Basic Model Building

```python
from celline import Project
from celline.functions.predict_celltype import BuildCellTypeModel

# Create project
project = Project("./my-project")

# Build basic cell type model
build_function = BuildCellTypeModel(
    species="Homo sapiens",
    suffix="brain_atlas",
    nthread=8,
    h5matrix_path="/data/reference_matrix.h5",
    celltype_path="/data/cell_annotations.tsv"
)

# Execute function
result = project.call(build_function)
```

#### Multi-Species Model

```python
from celline import Project
from celline.functions.predict_celltype import BuildCellTypeModel

# Create project
project = Project("./my-project")

# Build models for different species
species_configs = [
    {
        "species": "Homo sapiens",
        "suffix": "pbmc_v1",
        "matrix": "/data/human_pbmc.h5",
        "annotations": "/data/human_celltypes.tsv"
    },
    {
        "species": "Mus musculus", 
        "suffix": "brain_v2",
        "matrix": "/data/mouse_brain.h5",
        "annotations": "/data/mouse_celltypes.tsv"
    }
]

for config in species_configs:
    build_function = BuildCellTypeModel(
        species=config["species"],
        suffix=config["suffix"],
        nthread=16,
        h5matrix_path=config["matrix"],
        celltype_path=config["annotations"]
    )
    
    result = project.call(build_function)
```

#### Tissue-Specific Models

```python
from celline import Project
from celline.functions.predict_celltype import BuildCellTypeModel

# Create project
project = Project("./my-project")

# Build tissue-specific models
tissues = ["brain", "heart", "liver", "lung"]

for tissue in tissues:
    build_function = BuildCellTypeModel(
        species="Homo sapiens",
        suffix=f"{tissue}_atlas_v3",
        nthread=12,
        h5matrix_path=f"/data/{tissue}_reference.h5",
        celltype_path=f"/data/{tissue}_annotations.tsv"
    )
    
    result = project.call(build_function)
```

### CLI Usage

```bash
# Basic model building
celline run buildcelltypemodel \
    --species "Homo sapiens" \
    --suffix brain_v1 \
    --nthread 8 \
    --matrix /data/reference.h5 \
    --celltype /data/annotations.tsv

# High-performance model building
celline run buildcelltypemodel \
    --species "Mus musculus" \
    --suffix cortex_detailed \
    --nthread 16 \
    --matrix /data/mouse_cortex.h5 \
    --celltype /data/cortex_celltypes.tsv
```

## Implementation Details

### Input Data Requirements

#### Matrix File Format

Supported matrix formats:
- **HDF5 (.h5)**: Cell Ranger output format
- **H5Seurat (.h5seurat)**: Seurat HDF5 format  
- **H5SeuratV5 (.h5seuratv5)**: Seurat v5 format
- **Loom (.loom)**: Hierarchical data format

#### Cell Type Annotation Format

Required TSV file structure:

```tsv
cell	celltype
10X82_2_TCTCTCACCAGTTA	Astrocyte
10X82_2_TCTCTCACCAGTTC	Oligodendrocyte  
10X82_2_TCTCTCACCAGTTT	Neuron
PBMC_1_AAACCTGAGAAGGCCT	T_cell
PBMC_1_AAACCTGAGAAGGCCG	B_cell
```

**Column Requirements**:
- `cell`: Unique cell identifier (must match matrix)
- `celltype`: Cell type label (consistent naming)

### Validation Checks

The function performs comprehensive input validation:

```python
# File extension validation
valid_extensions = [".h5", ".loom", ".h5seurat", ".h5seuratv5"]
if not any(h5matrix_path.endswith(ext) for ext in valid_extensions):
    raise ValueError("Invalid matrix file format")

# Annotation file validation
if not celltype_path.endswith(".tsv"):
    raise ValueError("Annotation file must be TSV format")

# Column structure validation
df = pl.read_csv(celltype_path, separator="\t")
if df.columns != ["cell", "celltype"]:
    raise ValueError("Columns must be ['cell', 'celltype']")
```

### Model Training Pipeline

The training process follows these steps:

1. **Data Loading**: Loads reference matrix and cell annotations
2. **Quality Control**: Filters low-quality cells and genes
3. **Normalization**: Applies standard normalization procedures
4. **Feature Selection**: Identifies informative genes for classification
5. **Model Training**: Trains scPred classifier models
6. **Cross-Validation**: Validates model performance
7. **Model Export**: Saves trained models for deployment

### R Script Generation

The function generates optimized R training scripts:

```r
library(Seurat)
library(scPred)
library(Matrix)

# Load reference data
reference_matrix <- Read10X_h5("reference_matrix.h5")
cell_annotations <- read.delim("cell_annotations.tsv")

# Create Seurat object
reference <- CreateSeuratObject(
    counts = reference_matrix,
    min.cells = 3,
    min.features = 200
)

# Add cell type metadata
reference$celltype <- cell_annotations$celltype[
    match(colnames(reference), cell_annotations$cell)
]

# Preprocessing
reference <- NormalizeData(reference)
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)
reference <- RunPCA(reference)

# Train scPred model
reference <- getFeatureSpace(reference, "celltype")
reference <- trainModel(reference)

# Save trained model
saveRDS(reference, "reference.h5seurat")
saveRDS(reference@misc$scPred, "reference.pred")
```

## Model Output Structure

### Generated Files

Models are saved in organized directories:

```
project_root/
├── reference/
│   └── Species_Name/
│       └── model_suffix/
│           ├── reference.h5seurat    # Seurat reference object
│           ├── reference.pred        # scPred model object
│           ├── build.log            # Training log
│           ├── build.sh             # Generated R script
│           └── metrics/
│               ├── training_metrics.csv
│               ├── confusion_matrix.csv
│               └── feature_importance.csv
```

### Model Metadata

Each model includes comprehensive metadata:

```json
{
    "model_info": {
        "species": "Homo sapiens",
        "suffix": "brain_atlas_v2",
        "training_date": "2024-01-15T10:30:00Z",
        "celline_version": "1.0.0"
    },
    "data_info": {
        "n_cells": 50000,
        "n_genes": 20000,
        "n_celltypes": 15,
        "matrix_format": "h5"
    },
    "training_params": {
        "method": "scPred",
        "nthread": 8,
        "cv_folds": 5,
        "feature_selection": "seurat"
    },
    "performance": {
        "accuracy": 0.92,
        "f1_score": 0.89,
        "confusion_matrix": "confusion_matrix.csv"
    }
}
```

## Performance Optimization

### Memory Management

Memory requirements by dataset size:

| Dataset Size | Memory Needed | Recommended Threads |
|--------------|---------------|-------------------|
| 10K cells | 16 GB | 4-8 |
| 50K cells | 64 GB | 8-12 |
| 100K cells | 128 GB | 12-16 |
| 200K+ cells | 256+ GB | 16+ |

### Training Time

Typical training times:

| Dataset Size | Cell Types | Training Time |
|--------------|------------|---------------|
| 10K cells | 5-10 types | 30-60 minutes |
| 50K cells | 10-20 types | 2-4 hours |
| 100K cells | 20+ types | 4-8 hours |

### Cluster Computing

For large-scale model training:

```bash
#!/bin/bash
#PBS -N BuildCellTypeModel
#PBS -l nodes=1:ppn=16
#PBS -l mem=128gb
#PBS -l walltime=8:00:00

module load R/4.1.0
Rscript build_reference.sh
```

## Quality Control

### Training Metrics

The function evaluates model quality using standard metrics:

| Metric | Description | Good Value |
|--------|-------------|------------|
| **Accuracy** | Overall prediction accuracy | > 0.85 |
| **F1 Score** | Harmonic mean of precision/recall | > 0.80 |
| **Per-Class Precision** | Class-specific accuracy | > 0.75 |
| **Cross-Validation Score** | Generalization performance | > 0.80 |

### Feature Importance

```r
# Extract informative genes
important_features <- getFeatureSpace(reference)
top_genes <- head(important_features, 100)

# Visualize feature importance
plot_feature_importance(important_features)
```

## Error Handling

### Input Validation

```python
def validate_inputs(h5matrix_path, celltype_path):
    # Check file existence
    if not os.path.exists(h5matrix_path):
        raise FileNotFoundError(f"Matrix file not found: {h5matrix_path}")
    
    if not os.path.exists(celltype_path):
        raise FileNotFoundError(f"Annotation file not found: {celltype_path}")
    
    # Validate file formats
    if not h5matrix_path.endswith(('.h5', '.h5seurat', '.loom')):
        raise ValueError("Invalid matrix file format")
    
    # Check annotation structure
    df = pd.read_csv(celltype_path, sep='\t')
    if list(df.columns) != ['cell', 'celltype']:
        raise ValueError("Invalid annotation file structure")
```

### Common Issues

1. **Cell ID Mismatch**: Cell IDs in matrix and annotations don't match
2. **Insufficient Data**: Too few cells per cell type for training
3. **Memory Exhaustion**: Dataset too large for available memory
4. **R Dependencies**: Missing required R packages

## Methods

### `call(project: Project) -> Project`

Main execution method that builds the cell type model.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Validates input files and parameters
2. Creates output directory structure
3. Generates R training script with job configuration
4. Executes model training pipeline
5. Validates model output

### `__show_help()`

Displays help information with example data format and usage instructions.

## Best Practices

### Reference Data Selection

1. **High Quality**: Use well-characterized, high-quality reference datasets
2. **Cell Type Balance**: Include sufficient cells per type (≥100 cells)
3. **Biological Relevance**: Choose references relevant to target datasets
4. **Protocol Consistency**: Match sequencing protocols when possible

### Annotation Guidelines

```python
# Good cell type naming conventions
good_names = [
    "T_cell", "B_cell", "NK_cell",           # Immune
    "Neuron", "Astrocyte", "Oligodendrocyte", # Brain
    "Hepatocyte", "Kupffer_cell"              # Liver
]

# Avoid inconsistent naming
avoid_names = [
    "t cell", "T-cell", "T Cell",            # Inconsistent spacing/case
    "Neuron1", "Neuron2",                    # Arbitrary numbering
    "Unknown", "Other", "Unassigned"         # Non-specific labels
]
```

### Model Versioning

```python
# Systematic model versioning
model_versions = {
    "v1.0": "Initial model with basic cell types",
    "v1.1": "Added rare cell types",
    "v2.0": "Updated with improved annotations",
    "v2.1": "Performance optimizations"
}

# Use descriptive suffixes
suffixes = [
    "pbmc_comprehensive",    # Broad tissue coverage
    "brain_detailed",        # High resolution
    "immune_specialized",    # Focused application
    "development_timecourse" # Temporal aspects
]
```

## Related Functions

- [PredictCelltype](predict_celltype) - Use built models for prediction
- [Preprocess](preprocess) - Quality control before model building
- [SetTranscriptome](settranscriptome) - Reference genome setup
- [Count](count) - Generate count matrices for training

## Troubleshooting

### Common Issues

1. **Model Training Failure**: Check R package versions and dependencies
2. **Poor Model Performance**: Evaluate reference data quality and balance
3. **Memory Errors**: Reduce dataset size or use cluster computing
4. **File Format Errors**: Validate input file formats and structure

### Debug Mode

Enable detailed R logging:

```r
# In R script
options(error = traceback)
sessionInfo()

# Debug scPred training
debug(trainModel)
```

### Performance Tuning

```python
# Optimize for large datasets
def optimize_training(n_cells, n_celltypes):
    if n_cells > 100000:
        # Use subsampling for very large datasets
        subsample_size = min(50000, n_cells)
        nthread = min(16, n_celltypes * 2)
    else:
        subsample_size = n_cells
        nthread = min(8, n_celltypes)
    
    return subsample_size, nthread
```

### Manual Model Building

For debugging, build models manually:

```bash
# Navigate to model directory  
cd reference/Homo_sapiens/brain_atlas/

# Execute build script manually
Rscript build.sh

# Check outputs
ls -la *.rds *.pred
```