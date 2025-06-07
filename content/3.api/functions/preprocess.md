# Preprocess Function

Perform quality control and cell type filtering on counted single-cell data.

## Overview

The `Preprocess` function performs comprehensive quality control on Cell Ranger count output. It calculates quality metrics, detects doublets using Scrublet, and filters cells based on gene count, mitochondrial content, and cell type criteria.

## Class Information

- **Module**: `celline.functions.preprocess`
- **Class**: `Preprocess`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `target_celltype` | `Optional[list[str]]` | No | List of target cell types to include in filtering |

## Usage Examples

### Python API

#### Basic Preprocessing

```python
from celline import Project
from celline.functions.preprocess import Preprocess

# Create project
project = Project("./my-project")

# Basic preprocessing (includes all cell types)
preprocess_function = Preprocess()

# Execute function
result = project.call(preprocess_function)
```

#### Cell Type-Specific Preprocessing

```python
from celline import Project
from celline.functions.preprocess import Preprocess

# Create project
project = Project("./my-project")

# Preprocess with specific cell types
target_celltypes = ["Neuron", "Astrocyte", "Oligodendrocyte"]
preprocess_function = Preprocess(target_celltype=target_celltypes)

# Execute function
result = project.call(preprocess_function)
```

#### Multiple Preprocessing Scenarios

```python
from celline import Project
from celline.functions.preprocess import Preprocess

# Create project
project = Project("./my-project")

# Scenario 1: All cell types
project.call(Preprocess())

# Scenario 2: Neuronal cells only
neuronal_types = ["Neuron", "Interneuron", "Motor_neuron"]
project.call(Preprocess(target_celltype=neuronal_types))

# Scenario 3: Immune cells only
immune_types = ["T_cell", "B_cell", "Macrophage", "NK_cell"]
project.call(Preprocess(target_celltype=immune_types))
```

### CLI Usage

#### Basic Usage

```bash
# Preprocess all cell types
celline run preprocess

# Preprocess specific cell types
celline run preprocess --target-celltype Neuron Astrocyte

# Multiple cell types
celline run preprocess --target-celltype T_cell B_cell Macrophage NK_cell
```

### CLI Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `--target-celltype`, `-t` | `str+` | Target cell types to include in preprocessing |

## Implementation Details

### Prerequisites

The function requires samples to be:
1. **Counted**: Cell Ranger count must be completed
2. **Cell Type Predicted**: Cell type prediction must be finished

### Quality Control Pipeline

The preprocessing pipeline performs the following steps:

1. **Data Loading**: Reads Cell Ranger HDF5 output using Scanpy
2. **Metadata Integration**: Joins count data with cell type predictions
3. **Doublet Detection**: Uses Scrublet to identify potential doublets
4. **QC Metrics Calculation**: Computes standard quality control metrics
5. **Cell Filtering**: Applies multi-criteria filtering
6. **Output Generation**: Saves filtered cell information

### Filtering Criteria

The function applies the following quality control filters:

| Filter | Criteria | Purpose |
|--------|----------|---------|
| **Gene Count (Min)** | ≥ 200 genes | Remove low-quality cells |
| **Gene Count (Max)** | ≤ 5000 genes | Remove potential multiplets |
| **Mitochondrial Content** | ≤ 5% | Remove dying/stressed cells |
| **Doublet Detection** | Not predicted doublet | Remove computational doublets |
| **Cell Type** | In target list | Include only specified types |

### Scanpy Integration

The function uses Scanpy for data processing:

```python
import scanpy as sc

# Read Cell Ranger output
adata = sc.read_10x_h5(path.resources_sample_counted)

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)
```

### Scrublet Doublet Detection

Doublet detection using Scrublet:

```python
import scrublet as scr

# Initialize Scrublet
scrub = scr.Scrublet(adata.X)

# Detect doublets
doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)

# Add to observation metadata
adata.obs["doublet_score"] = doublet_scores
adata.obs["predicted_doublets"] = predicted_doublets
```

## Output Files

### Cell Information File

The function generates a comprehensive cell information file:

```
data/SAMPLE_ID/cell_info.tsv
```

#### File Structure

| Column | Type | Description |
|--------|------|-------------|
| `barcode` | `str` | Cell barcode sequence |
| `project` | `str` | Project identifier |
| `sample` | `str` | Sample identifier |
| `cell` | `str` | Unique cell identifier |
| `cell_type` | `str` | Predicted cell type |
| `include` | `bool` | Pass QC filters |
| `n_genes_by_counts` | `int` | Number of detected genes |
| `total_counts` | `int` | Total UMI count |
| `pct_counts_mt` | `float` | Mitochondrial gene percentage |
| `doublet_score` | `float` | Doublet probability score |
| `predicted_doublets` | `bool` | Doublet prediction |

#### Example Output

```tsv
barcode	project	sample	cell	cell_type	include	n_genes_by_counts	total_counts	pct_counts_mt	doublet_score	predicted_doublets
AAACCTGAGAAGGCCT-1	GSE123456	GSM789012	GSM789012_1	Neuron	true	2450	8924	2.1	0.05	false
AAACCTGAGAAGGCCT-2	GSE123456	GSM789012	GSM789012_2	Astrocyte	true	1876	5234	1.8	0.03	false
AAACCTGAGAAGGCCT-3	GSE123456	GSM789012	GSM789012_3	Doublet	false	5421	15673	8.2	0.92	true
```

## Quality Control Metrics

### Standard Metrics

The function calculates standard single-cell QC metrics:

| Metric | Description | Typical Range |
|--------|-------------|---------------|
| `n_genes_by_counts` | Genes with non-zero counts | 200-5000 |
| `total_counts` | Total UMI count per cell | 1000-50000 |
| `pct_counts_mt` | Mitochondrial gene percentage | 0-20% |
| `doublet_score` | Scrublet doublet score | 0-1 |

### Mitochondrial Gene Detection

Mitochondrial genes are identified by prefix:

```python
# Default mitochondrial gene prefix
mt_prefix = "mt-"  # For mouse
# For human data, use "MT-"

# Mark mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
```

### Cell Type Integration

Cell types from prediction are integrated with QC data:

```python
# Read cell type predictions
celltype_data = pl.read_csv(
    path.data_sample_predicted_celltype,
    separator="\t"
).rename({"scpred_prediction": "cell_type"})

# Join with observation data
obs = obs.join(celltype_data, on="cell")
```

## Error Handling

### Common Issues

1. **Missing Count Data**: Ensure Cell Ranger count completed successfully
2. **Missing Cell Type Predictions**: Run cell type prediction first
3. **Memory Issues**: Large datasets may require substantial RAM
4. **File Permissions**: Check read/write permissions

### Validation Checks

```python
# Check prerequisites
if not (path.is_counted and path.is_predicted_celltype):
    print(f"Sample {sample_id} not ready for preprocessing")
    continue

# Validate cell type predictions exist
if not os.path.exists(path.data_sample_predicted_celltype):
    raise FileNotFoundError("Cell type predictions not found")
```

## Performance Considerations

### Memory Usage

Memory requirements depend on dataset size:

| Dataset Size | Memory Needed |
|--------------|---------------|
| <5K cells | 4 GB |
| 5K-20K cells | 8 GB |
| 20K-50K cells | 16 GB |
| >50K cells | 32+ GB |

### Processing Time

Typical processing times:

| Dataset Size | Processing Time |
|--------------|----------------|
| <5K cells | 30 seconds |
| 5K-20K cells | 1-2 minutes |
| 20K-50K cells | 2-5 minutes |
| >50K cells | 5+ minutes |

### Optimization Tips

1. **Filter Early**: Apply cell type filters before heavy computations
2. **Memory Management**: Process samples sequentially for large datasets
3. **Batch Processing**: Group similar samples for efficiency

## Methods

### `call(project: Project) -> Project`

Main execution method that performs preprocessing on all eligible samples.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Reads samples from `samples.toml`
2. Validates prerequisites for each sample
3. Loads count data and cell type predictions
4. Calculates QC metrics and detects doublets
5. Applies filtering criteria
6. Saves filtered cell information

### `add_cli_args(parser: ArgumentParser) -> None`

Adds CLI-specific arguments to the argument parser.

### `cli(project: Project, args: Namespace) -> Project`

CLI entry point that processes command-line arguments and executes the function.

## Integration with Pipeline

### Typical Workflow

```python
from celline import Project
from celline.functions.count import Count
from celline.functions.predict_celltype import PredictCelltype
from celline.functions.preprocess import Preprocess
from celline.functions.create_seurat import CreateSeuratObject

# Complete pipeline
project = Project("./my-project")

# Count reads
project.call(Count(nthread=8))

# Predict cell types
project.call(PredictCelltype(model))

# Preprocess with QC
project.call(Preprocess(target_celltype=["Neuron", "Astrocyte"]))

# Create Seurat objects
project.call(CreateSeuratObject(useqc_matrix=True))
```

## Cell Type Filtering

### Supported Cell Types

The function supports any cell types from your prediction model:

**Common Neuroscience Types**:
- Neuron, Astrocyte, Oligodendrocyte, Microglia
- Interneuron, Pyramidal_neuron, Motor_neuron
- OPC (Oligodendrocyte Precursor Cell)

**Common Immunology Types**:
- T_cell, B_cell, NK_cell, Macrophage
- CD4_T_cell, CD8_T_cell, Regulatory_T_cell
- Plasma_cell, Dendritic_cell

**Common Development Types**:
- Stem_cell, Progenitor_cell
- Epithelial_cell, Endothelial_cell
- Fibroblast, Smooth_muscle_cell

### Custom Cell Type Lists

```python
# Create custom cell type groups
epithelial_types = [
    "Epithelial_cell", "Basal_cell", "Ciliated_cell",
    "Goblet_cell", "Club_cell"
]

# Preprocess with custom types
project.call(Preprocess(target_celltype=epithelial_types))
```

## Related Functions

- [Count](count) - Generate count matrices before preprocessing
- [PredictCelltype](predict_celltype) - Predict cell types before preprocessing
- [CreateSeuratObject](createseuratobject) - Create Seurat objects after preprocessing
- [Integrate](integrate) - Integrate preprocessed samples

## Troubleshooting

### Common Issues

1. **Memory Error**: Reduce dataset size or increase system memory
2. **Missing Dependencies**: Install required Python packages (scanpy, scrublet, polars)
3. **File Not Found**: Ensure count and prediction steps completed successfully
4. **Empty Output**: Check cell type names match prediction output

### Debug Mode

Enable detailed logging:

```python
import logging
logging.basicConfig(level=logging.DEBUG)

# Also enable scanpy logging
import scanpy as sc
sc.settings.verbosity = 3
```

### Manual Quality Control

For custom QC parameters:

```python
# Custom filtering criteria
custom_filter = {
    "min_genes": 500,      # Higher gene count threshold
    "max_genes": 3000,     # Lower gene count threshold  
    "max_mt_pct": 10,      # Higher mitochondrial threshold
    "min_counts": 1000     # Minimum UMI count
}

# Apply custom filters in your analysis
```