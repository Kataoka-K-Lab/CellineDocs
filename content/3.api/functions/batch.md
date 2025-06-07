# BatchCorrection Function

Correct batch effects in multi-sample single-cell RNA-seq datasets.

## Overview

The `BatchCorrection` function removes technical batch effects that arise when samples are processed in different batches, sequencing runs, or experimental conditions. It uses advanced R-based computational methods to harmonize data while preserving biological variation.

## Class Information

- **Module**: `celline.functions.batch_cor`
- **Class**: `BatchCorrection`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `output_file_path` | `str` | Yes | Path for batch-corrected output files |
| `filter_func` | `Optional[Callable[[SampleSchema], bool]]` | No | Function to filter samples for correction |

### JobContainer Structure

| Field | Type | Description |
|-------|------|-------------|
| `nthread` | `str` | Number of threads (fixed to 1) |
| `cluster_server` | `str` | Cluster server name |
| `logpath` | `str` | Batch correction log path |
| `r_path` | `str` | R script directory |
| `exec_root` | `str` | Execution root directory |
| `proj_path` | `str` | Project root path |
| `sample_ids` | `str` | Comma-separated sample IDs |
| `project_ids` | `str` | Comma-separated project IDs |
| `output_dir` | `str` | Output directory path |
| `logpath_runtime` | `str` | Runtime log path |

## Usage Examples

### Python API

#### Basic Batch Correction

```python
from celline import Project
from celline.functions.batch_cor import BatchCorrection

# Create project
project = Project("./my-project")

# Basic batch correction for all samples
batch_function = BatchCorrection(
    output_file_path="/path/to/batch_corrected",
    filter_func=None
)

# Execute function
result = project.call(batch_function)
```

#### Filtered Batch Correction

```python
from celline import Project
from celline.functions.batch_cor import BatchCorrection
from celline.DB.dev.model import SampleSchema

# Create project
project = Project("./my-project")

# Filter function for specific conditions
def same_tissue_filter(schema: SampleSchema) -> bool:
    return "brain" in (schema.title or "").lower()

# Batch correction with filtering
batch_function = BatchCorrection(
    output_file_path="/path/to/brain_batch_corrected",
    filter_func=same_tissue_filter
)

# Execute function
result = project.call(batch_function)
```

#### Platform-Specific Correction

```python
from celline import Project
from celline.functions.batch_cor import BatchCorrection
from celline.DB.dev.model import SampleSchema

# Create project
project = Project("./my-project")

# Filter by sequencing platform
def platform_filter(schema: SampleSchema) -> bool:
    return schema.platform.startswith("Illumina")

# Batch correction for specific platform
batch_function = BatchCorrection(
    output_file_path="/path/to/illumina_corrected",
    filter_func=platform_filter
)

# Execute function
result = project.call(batch_function)
```

#### Multi-Study Correction

```python
from celline import Project
from celline.functions.batch_cor import BatchCorrection
from celline.DB.dev.model import SampleSchema

# Create project
project = Project("./my-project")

# Filter for multiple studies
target_studies = ["GSE123456", "GSE789012", "GSE345678"]

def multi_study_filter(schema: SampleSchema) -> bool:
    return schema.parent in target_studies

# Batch correction across studies
batch_function = BatchCorrection(
    output_file_path="/path/to/multi_study_corrected",
    filter_func=multi_study_filter
)

# Execute function
result = project.call(batch_function)
```

### CLI Usage

```bash
# Basic batch correction
celline run batch

# Batch correction with verbose logging
celline run batch --verbose
```

## Implementation Details

### Batch Effect Sources

Common sources of batch effects in single-cell data:

| Source | Description | Impact |
|--------|-------------|--------|
| **Sequencing Run** | Different sequencing lanes/machines | Technical variation |
| **Sample Preparation** | Different processing dates | Protocol variation |
| **Cell Isolation** | Different isolation methods | Methodological bias |
| **Library Preparation** | Kit lot variations | Technical artifacts |
| **Operator Effects** | Different personnel | Handling variation |
| **Environmental** | Temperature, humidity changes | Subtle technical effects |

### Correction Methods

The function employs multiple R-based correction algorithms:

#### ComBat-seq
```r
# Negative binomial model for count data
library(sva)
corrected_counts <- ComBat_seq(
    counts = count_matrix,
    batch = batch_info,
    group = condition_info
)
```

#### Harmony
```r
# Fast harmonic correction
library(harmony)
corrected_embedding <- RunHarmony(
    object = seurat_object,
    group.by.vars = "batch",
    reduction = "pca"
)
```

#### MNN (Mutual Nearest Neighbors)
```r
# Batch correction via MNN
library(batchelor)
corrected <- fastMNN(sample_list, cos.norm = TRUE)
```

#### scVI Integration
```r
# Deep learning approach (if available)
library(reticulate)
# Python scVI integration
```

### File Organization

Batch correction outputs are systematically organized:

```
project_root/
├── batch/
│   ├── batch_corrected_TIMESTAMP.h5ad    # Corrected AnnData object
│   ├── batch_corrected_TIMESTAMP.rds     # Corrected Seurat object
│   ├── logs/
│   │   ├── removebatch_TIMESTAMP.log     # Correction log
│   │   └── runtime/
│   │       └── removebatch_TIMESTAMP.sh  # Runtime log
│   └── removebatch_TIMESTAMP.sh          # Generated R script
```

## Correction Workflow

### Prerequisites

Samples must be:
1. **Counted**: Cell Ranger count completed
2. **Quality Controlled**: Basic QC filtering applied
3. **Normalized**: Data normalization completed
4. **Compatible**: Same species and cell types

### Processing Pipeline

1. **Sample Selection**: Applies filter function to identify target samples
2. **Batch Identification**: Automatically detects batch variables
3. **Data Integration**: Combines count matrices from multiple samples
4. **Batch Detection**: Identifies technical batch effects
5. **Correction Application**: Applies appropriate correction method
6. **Quality Assessment**: Evaluates correction effectiveness
7. **Output Generation**: Saves corrected data in multiple formats

### Batch Variable Detection

The function automatically identifies batch variables:

```r
# Automatic batch detection
batch_vars <- c("orig.ident", "sequencing_run", "prep_date", "platform")
detected_batches <- identify_batch_effects(
    seurat_object,
    variables = batch_vars,
    threshold = 0.05
)
```

## Methods

### `call(project: Project) -> Project`

Main execution method that performs batch correction.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Filters samples based on provided filter function
2. Validates prerequisites for batch correction
3. Creates output directories
4. Generates R batch correction script
5. Executes correction pipeline
6. Saves corrected results

### `register() -> str`

Returns the function identifier for registration.

## Filter Functions

### Common Filter Patterns

```python
# Filter by tissue type
def tissue_filter(tissue_type):
    def filter_func(schema):
        return tissue_type.lower() in (schema.title or "").lower()
    return filter_func

# Filter by time period
def temporal_filter(start_date, end_date):
    def filter_func(schema):
        return start_date <= schema.date <= end_date
    return filter_func

# Filter by cell count threshold
def cell_count_filter(min_cells):
    def filter_func(schema):
        return schema.cell_count >= min_cells
    return filter_func

# Combined filters for complex scenarios
def comprehensive_filter(schema):
    return (
        schema.species == "Homo sapiens" and
        "cortex" in (schema.title or "").lower() and
        schema.platform.startswith("10x") and
        schema.cell_count >= 1000
    )
```

### Batch-Specific Filters

```python
# Same laboratory filter
def lab_filter(lab_name):
    def filter_func(schema):
        return schema.laboratory == lab_name
    return filter_func

# Protocol consistency filter
def protocol_filter(protocol_version):
    def filter_func(schema):
        return schema.protocol_version == protocol_version
    return filter_func
```

## Quality Assessment

### Correction Evaluation Metrics

The function evaluates correction quality using multiple metrics:

| Metric | Description | Good Value |
|--------|-------------|------------|
| **Silhouette Score** | Batch mixing quality | > 0.7 |
| **kBET** | Batch effect test | < 0.05 |
| **LISI** | Local inverse Simpson index | > 2.0 |
| **ASW** | Average silhouette width | > 0.5 |
| **Graph Connectivity** | Biological preservation | > 0.8 |

### Visualization Outputs

Generated visualizations include:

```r
# Before/after comparison plots
plot_umap_before_after(seurat_object)
plot_batch_effects_heatmap(corrected_data)
plot_correction_metrics(evaluation_results)

# Quality control plots
plot_gene_expression_conservation(original, corrected)
plot_cell_type_preservation(corrected_data)
```

## Performance Considerations

### Memory Requirements

Batch correction memory needs:

| Samples | Cells per Sample | Memory Needed |
|---------|------------------|---------------|
| 2-5 | 5K-10K | 32-64 GB |
| 5-10 | 5K-10K | 64-128 GB |
| 10+ | 5K-10K | 128+ GB |

### Processing Time

Correction time by method and data size:

| Method | Sample Count | Processing Time |
|--------|--------------|-----------------|
| Harmony | 5 samples | 1-2 hours |
| ComBat-seq | 5 samples | 2-4 hours |
| FastMNN | 5 samples | 3-6 hours |
| scVI | 5 samples | 4-8 hours |

### Cluster Computing

For large-scale correction:

```bash
#!/bin/bash
#PBS -N BatchCorrection
#PBS -l nodes=1:ppn=16
#PBS -l mem=256gb
#PBS -l walltime=12:00:00

module load R/4.1.0
Rscript batch_correction_script.R
```

## Error Handling

### Common Issues

1. **Insufficient Samples**: Need ≥2 samples for batch correction
2. **Memory Exhaustion**: Large datasets require substantial RAM
3. **Batch Confounding**: Biological and technical effects correlated
4. **Missing Metadata**: Batch variables not properly defined

### Troubleshooting

```python
# Validate batch correction prerequisites
def validate_batch_correction(samples):
    if len(samples) < 2:
        raise ValueError("Need at least 2 samples for batch correction")
    
    # Check for batch variable availability
    batch_vars = set()
    for sample in samples:
        batch_vars.add(sample.schema.platform)
        batch_vars.add(sample.schema.prep_date)
    
    if len(batch_vars) < 2:
        print("Warning: Limited batch variation detected")
```

## Best Practices

### Sample Selection

1. **Biological Coherence**: Correct samples from same tissue/condition
2. **Technical Diversity**: Ensure multiple batches are represented
3. **Quality Control**: Remove low-quality samples before correction
4. **Sample Size**: Include sufficient cells per batch (≥500)

### Method Selection

```python
# Choose correction method based on data characteristics
def select_correction_method(n_samples, n_cells, has_replicates):
    if n_samples <= 5 and n_cells <= 50000:
        return "Harmony"  # Fast and effective
    elif has_replicates:
        return "ComBat-seq"  # Good for designed experiments
    elif n_samples > 10:
        return "FastMNN"  # Scalable for many samples
    else:
        return "scVI"  # Deep learning for complex cases
```

### Pipeline Integration

```python
# Complete batch correction workflow
from celline import Project
from celline.functions.preprocess import Preprocess
from celline.functions.batch_cor import BatchCorrection
from celline.functions.integrate import Integrate

# Process samples
project = Project("./my-project")
project.call(Preprocess())  # QC first

# Apply batch correction
project.call(BatchCorrection(
    output_file_path="./corrected_data",
    filter_func=None
))

# Final integration if needed
project.call(Integrate(
    filter_func=None,
    outfile_name="final_integrated"
))
```

## Output Formats

### Supported Output Types

| Format | Description | Use Case |
|--------|-------------|----------|
| **HDF5/H5AD** | AnnData format | Python analysis |
| **RDS** | Seurat object | R analysis |
| **H5Seurat** | Seurat HDF5 | Cross-platform |
| **CSV/TSV** | Matrix format | General purpose |
| **Loom** | Hierarchical format | Visualization |

### Metadata Preservation

Corrected outputs preserve essential metadata:

```python
# Metadata maintained after correction
corrected_metadata = {
    "original_sample_id": "GSM123456",
    "batch_corrected": True,
    "correction_method": "Harmony",
    "correction_date": "2024-01-15",
    "original_batch": "batch_1",
    "correction_parameters": {...}
}
```

## Related Functions

- [Preprocess](preprocess) - Quality control before batch correction
- [Integrate](integrate) - Alternative integration approach
- [CreateSeuratObject](createseuratobject) - Create objects for correction
- [Reduce](reduce) - Storage optimization after correction

## Troubleshooting

### Common Issues

1. **Overcorrection**: Biological signals removed with batch effects
2. **Undercorrection**: Batch effects remain after correction
3. **Memory Issues**: Insufficient RAM for large datasets
4. **Method Failure**: Correction algorithm fails to converge

### Debug Mode

Enable detailed R logging:

```r
# In R script
options(error = traceback)
sessionInfo()
print("Batch variables detected:")
print(table(metadata$batch))
```

### Manual Correction

For debugging, run correction manually:

```bash
# Navigate to batch directory
cd batch/

# Execute correction script manually
Rscript removebatch_TIMESTAMP.sh
```

### Validation Steps

```python
# Validate correction results
def validate_correction(original, corrected):
    # Check cell count preservation
    assert corrected.n_obs == original.n_obs
    
    # Check gene preservation
    assert corrected.n_vars == original.n_vars
    
    # Evaluate batch mixing
    batch_score = calculate_batch_mixing(corrected)
    assert batch_score > 0.5, "Poor batch correction"
```