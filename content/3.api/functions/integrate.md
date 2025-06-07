# Integrate Function

Integrate multiple single-cell RNA-seq samples using advanced computational methods.

## Overview

The `Integrate` function combines multiple single-cell RNA-seq samples into a unified analysis, correcting for batch effects and technical variations while preserving biological differences. It uses R-based integration methods to create harmonized datasets.

## Class Information

- **Module**: `celline.functions.integrate`
- **Class**: `Integrate`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `filter_func` | `Optional[Callable[[SampleSchema], bool]]` | No | Function to filter samples for integration |
| `outfile_name` | `Optional[str]` | No | Custom output filename (default: timestamped) |

### JobContainer Structure

| Field | Type | Description |
|-------|------|-------------|
| `nthread` | `str` | Number of threads (fixed to 1) |
| `cluster_server` | `str` | Cluster server name |
| `jobname` | `str` | Job identifier |
| `logpath` | `str` | Integration log file path |
| `r_path` | `str` | R script directory |
| `exec_root` | `str` | Execution root |
| `sample_ids` | `str` | Comma-separated sample IDs |
| `project_ids` | `str` | Comma-separated project IDs |
| `all_bcmat_path` | `str` | Comma-separated matrix paths |
| `all_data_sample_dir_path` | `str` | Comma-separated data directories |
| `outfile_path` | `str` | Output file path |
| `logpath_runtime` | `str` | Runtime log path |
| `project_name` | `str` | Project name |

## Usage Examples

### Python API

#### Basic Integration

```python
from celline import Project
from celline.functions.integrate import Integrate

# Create project
project = Project("./my-project")

# Basic integration of all samples
integrate_function = Integrate(filter_func=None, outfile_name=None)

# Execute function
result = project.call(integrate_function)
```

#### Custom Output Name

```python
from celline import Project
from celline.functions.integrate import Integrate

# Create project
project = Project("./my-project")

# Integration with custom output name
integrate_function = Integrate(
    filter_func=None,
    outfile_name="brain_atlas_integrated"
)

# Execute function
result = project.call(integrate_function)
```

#### Filtered Integration

```python
from celline import Project
from celline.functions.integrate import Integrate
from celline.DB.dev.model import SampleSchema

# Create project
project = Project("./my-project")

# Filter function to include only specific samples
def brain_samples_filter(schema: SampleSchema) -> bool:
    return "brain" in schema.title.lower() if schema.title else False

# Integration with sample filtering
integrate_function = Integrate(
    filter_func=brain_samples_filter,
    outfile_name="brain_samples_only"
)

# Execute function
result = project.call(integrate_function)
```

#### Species-Specific Integration

```python
from celline import Project
from celline.functions.integrate import Integrate
from celline.DB.dev.model import SampleSchema

# Create project
project = Project("./my-project")

# Filter by species
def human_samples_filter(schema: SampleSchema) -> bool:
    return schema.species == "Homo sapiens"

# Integration of human samples only
integrate_function = Integrate(
    filter_func=human_samples_filter,
    outfile_name="human_integration"
)

# Execute function
result = project.call(integrate_function)
```

### CLI Usage

```bash
# Basic integration
celline run integrate

# Integration with verbose logging
celline run integrate --verbose
```

## Implementation Details

### Integration Methods

The function uses state-of-the-art R-based integration methods:

- **Seurat Integration**: CCA-based integration
- **Harmony**: Fast integration for large datasets
- **FastMNN**: Mutual nearest neighbors approach
- **scVI**: Deep learning-based integration (if available)

### R Script Integration

The function generates and executes R scripts:

```r
# Load required libraries
library(Seurat)
library(harmony)
library(SeuratWrappers)

# Load sample data
sample_list <- list()
for (i in 1:length(sample_paths)) {
    sample_list[[i]] <- Read10X_h5(sample_paths[i])
}

# Create Seurat objects
seurat_objects <- lapply(sample_list, function(x) {
    CreateSeuratObject(counts = x, min.cells = 3, min.features = 200)
})

# Integrate samples
integrated <- RunHarmony(
    object = merged_object,
    group.by.vars = "sample",
    reduction = "pca"
)

# Save integrated object
saveRDS(integrated, file = output_path)
```

### File Organization

Integration outputs are organized systematically:

```
project_root/
├── integration/
│   ├── integrated_TIMESTAMP.rds    # Integrated Seurat object
│   ├── logs/
│   │   ├── integrate_TIMESTAMP.log     # Integration log
│   │   └── RUNTIME_integrate_TIMESTAMP.log  # Runtime log
│   └── integrate_TIMESTAMP.sh      # Generated R script
```

## Integration Workflow

### Prerequisites

Samples must be:
1. **Counted**: Cell Ranger count completed
2. **Optionally Preprocessed**: QC filtering applied
3. **Compatible**: Same species and similar protocols

### Processing Steps

1. **Sample Discovery**: Identifies all counted samples
2. **Filtering**: Applies user-defined filter function
3. **Validation**: Checks sample compatibility
4. **Matrix Loading**: Loads count matrices from each sample
5. **Preprocessing**: Normalizes and scales data
6. **Integration**: Applies batch correction methods
7. **Dimensionality Reduction**: Computes integrated embeddings
8. **Output Generation**: Saves integrated object

### Quality Control

The integration process includes automatic QC:

```r
# Before integration
pre_integration_metrics <- list(
    n_samples = length(sample_list),
    total_cells = sum(sapply(sample_list, ncol)),
    total_genes = nrow(sample_list[[1]])
)

# After integration
post_integration_metrics <- list(
    integrated_cells = ncol(integrated_object),
    batch_correction_score = calculate_batch_score(integrated_object)
)
```

## Methods

### `call(project: Project) -> Project`

Main execution method that performs sample integration.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Discovers and filters samples
2. Validates prerequisites
3. Generates R integration script
4. Executes integration pipeline
5. Saves integrated results

### `register() -> str`

Returns the function identifier for registration.

## Filter Functions

### Common Filter Patterns

```python
# Filter by tissue type
def tissue_filter(tissue_name):
    def filter_func(schema):
        return tissue_name.lower() in (schema.title or "").lower()
    return filter_func

# Filter by platform
def platform_filter(platform_name):
    def filter_func(schema):
        return schema.platform == platform_name
    return filter_func

# Filter by date range
def date_filter(start_date, end_date):
    def filter_func(schema):
        return start_date <= schema.date <= end_date
    return filter_func

# Combined filters
def combined_filter(schema):
    return (
        schema.species == "Homo sapiens" and
        "brain" in (schema.title or "").lower() and
        schema.platform.startswith("Illumina")
    )
```

### Filter Examples

```python
# Use case: Integrate only cortical samples
cortex_filter = lambda s: "cortex" in (s.title or "").lower()
project.call(Integrate(filter_func=cortex_filter))

# Use case: Integrate samples from specific study
study_filter = lambda s: s.parent == "GSE123456"
project.call(Integrate(filter_func=study_filter))
```

## Performance Considerations

### Memory Requirements

Integration memory needs scale with sample number and size:

| Samples | Cells per Sample | Memory Needed |
|---------|------------------|---------------|
| 2-5 | 5K-10K | 16-32 GB |
| 5-10 | 5K-10K | 32-64 GB |
| 10+ | 5K-10K | 64+ GB |

### Processing Time

Integration time depends on method and data size:

| Method | Sample Count | Processing Time |
|--------|--------------|-----------------|
| Harmony | 5 samples | 30-60 minutes |
| Seurat CCA | 5 samples | 1-2 hours |
| FastMNN | 5 samples | 45-90 minutes |

### Cluster Computing

For large integrations:

```bash
#!/bin/bash
#PBS -N Integration
#PBS -l nodes=1:ppn=8
#PBS -l mem=128gb
#PBS -l walltime=8:00:00

module load R/4.1.0
Rscript integration_script.R
```

## Integration Methods

### Harmony Integration

Fast and effective for most datasets:

```r
# Harmony workflow
merged <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
integrated <- RunHarmony(merged, group.by.vars = "orig.ident")
```

### Seurat CCA Integration

Canonical Correlation Analysis approach:

```r
# CCA workflow
anchors <- FindIntegrationAnchors(object.list = seurat_list)
integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
```

### FastMNN Integration

Mutual nearest neighbors approach:

```r
# FastMNN workflow
library(batchelor)
integrated <- fastMNN(seurat_list, cos.norm = TRUE)
```

## Output Analysis

### Integrated Object Structure

The output Seurat object contains:

```r
# Assays
integrated@assays$RNA          # Original counts
integrated@assays$integrated   # Integrated data (if CCA)

# Reductions
integrated@reductions$pca      # Principal components
integrated@reductions$harmony  # Harmony embeddings
integrated@reductions$umap     # UMAP coordinates

# Metadata
integrated@meta.data$orig.ident     # Original sample identity
integrated@meta.data$sample_id      # Sample identifier
integrated@meta.data$batch          # Batch information
```

### Quality Assessment

Evaluate integration quality:

```r
# Visualize integration
DimPlot(integrated, group.by = "orig.ident")
DimPlot(integrated, group.by = "cell_type")

# Calculate batch metrics
library(harmony)
batch_score <- assess_batch_correction(integrated)
```

## Error Handling

### Common Issues

1. **Memory Exhaustion**: Reduce sample count or use cluster computing
2. **R Package Missing**: Install required packages (Seurat, harmony)
3. **Incompatible Samples**: Ensure same species and protocols
4. **File Not Found**: Verify all samples are counted

### Troubleshooting

```python
# Check sample compatibility
def check_compatibility(samples):
    species_set = set(s.schema.species for s in samples)
    if len(species_set) > 1:
        print(f"Warning: Multiple species found: {species_set}")
    
    platform_set = set(s.schema.platform for s in samples)
    if len(platform_set) > 1:
        print(f"Info: Multiple platforms: {platform_set}")
```

## Best Practices

### Sample Selection

1. **Biological Coherence**: Integrate related samples/conditions
2. **Technical Compatibility**: Use similar protocols and platforms
3. **Quality Control**: Ensure all samples pass QC
4. **Batch Consideration**: Be aware of batch effects

### Integration Strategy

```python
# Progressive integration for large projects
def progressive_integration(project, batch_size=5):
    """Integrate samples in batches to manage memory."""
    samples = list(project.samples.values())
    
    for i in range(0, len(samples), batch_size):
        batch = samples[i:i+batch_size]
        filter_func = lambda s: s.key in [b.schema.key for b in batch]
        
        project.call(Integrate(
            filter_func=filter_func,
            outfile_name=f"batch_{i//batch_size + 1}"
        ))
```

## Related Functions

- [Count](count) - Generate count matrices before integration
- [Preprocess](preprocess) - Quality control before integration
- [CreateSeuratObject](createseuratobject) - Create Seurat objects
- [Reduce](reduce) - Storage optimization after integration

## Troubleshooting

### Common Issues

1. **R Memory Limit**: Increase R memory or use cluster
2. **Package Conflicts**: Use consistent R package versions
3. **File Permissions**: Ensure write access to integration directory
4. **Sample Mismatch**: Verify sample metadata consistency

### Debug Mode

Enable detailed R logging:

```r
# In R script
options(error = traceback)
sessionInfo()
print(ls())  # List objects in environment
```

### Manual Integration

For debugging, run integration manually:

```bash
# Navigate to integration directory
cd integration/

# Execute integration script manually
Rscript integrate_TIMESTAMP.sh
```