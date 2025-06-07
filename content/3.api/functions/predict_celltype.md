# PredictCelltype Function

Predict cell types in single-cell RNA-seq data using pre-trained scPred models.

## Overview

The `PredictCelltype` function performs automated cell type prediction using scPred-based models. It can work with pre-built reference models or custom-trained models for species-specific cell type classification.

## Class Information

- **Module**: `celline.functions.predict_celltype`
- **Class**: `PredictCelltype`
- **Base Class**: `CellineFunction`

## Related Classes

### `BuildCellTypeModel`

Builds custom cell type prediction models from reference datasets.

- **Module**: `celline.functions.predict_celltype`
- **Class**: `BuildCellTypeModel`
- **Base Class**: `CellineFunction`

### `CellTypeModel`

Data structure for cell type model configuration.

```python
@dataclass
class CellTypeModel:
    species: str
    suffix: Optional[str]
```

## Parameters

### PredictCelltype Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `model` | `CellTypeModel` | Yes | Cell type model configuration |
| `re_predict` | `bool` | No | Re-predict previously processed samples (default: False) |

### BuildCellTypeModel Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `species` | `str` | Yes | Species name (e.g., "Homo sapiens") |
| `suffix` | `str` | Yes | Model identifier suffix |
| `nthread` | `int` | Yes | Number of threads for model building |
| `h5matrix_path` | `str` | Yes | Path to reference HDF5 matrix |
| `celltype_path` | `str` | Yes | Path to cell type annotation TSV |

## Usage Examples

### Python API

#### Basic Cell Type Prediction

```python
from celline import Project
from celline.functions.predict_celltype import PredictCelltype, CellTypeModel

# Create project
project = Project("./my-project")

# Define model
model = CellTypeModel(species="Homo sapiens", suffix="brain")

# Create prediction function
predict_function = PredictCelltype(model)

# Execute function
result = project.call(predict_function)
```

#### Re-prediction of Existing Samples

```python
from celline import Project
from celline.functions.predict_celltype import PredictCelltype, CellTypeModel

# Create project
project = Project("./my-project")

# Define model
model = CellTypeModel(species="Mus musculus", suffix="cortex")

# Re-predict with updated model
predict_function = PredictCelltype(model, re_predict=True)

# Execute function
result = project.call(predict_function)
```

#### Building Custom Cell Type Model

```python
from celline import Project
from celline.functions.predict_celltype import BuildCellTypeModel

# Create project
project = Project("./my-project")

# Build custom model
build_function = BuildCellTypeModel(
    species="Homo sapiens",
    suffix="custom_brain",
    nthread=8,
    h5matrix_path="/path/to/reference_matrix.h5",
    celltype_path="/path/to/cell_annotations.tsv"
)

# Execute function
result = project.call(build_function)
```

### CLI Usage

#### Basic Prediction

```bash
# Predict cell types using default model
celline run predict_celltype

# Predict with specific model
celline run predict_celltype --model brain --species "Homo sapiens"

# Re-predict existing samples
celline run predict_celltype --re-predict
```

#### Model Building

```bash
# Build custom cell type model
celline run buildcelltypemodel \
    --species "Homo sapiens" \
    --suffix brain_v2 \
    --nthread 8 \
    --matrix /path/to/reference.h5 \
    --celltype /path/to/annotations.tsv
```

## Implementation Details

### scPred Integration

The function uses scPred for cell type prediction:

```r
# R script execution
library(scPred)
library(Seurat)

# Load reference model
reference <- readRDS("reference.h5seurat")
predictions <- readRDS("reference.pred")

# Predict cell types
predicted <- scPredict(query_data, reference)
```

### Model Structure

Cell type models are stored in organized directories:

```
project_root/
├── reference/
│   └── Homo_sapiens/
│       ├── brain/
│       │   ├── reference.h5seurat    # Seurat reference object
│       │   ├── reference.pred        # scPred model
│       │   └── build.log            # Build log
│       └── custom_brain/
│           ├── reference.h5seurat
│           ├── reference.pred
│           └── build.log
```

### Species-Specific Processing

The function automatically filters samples by species:

```python
# Filter samples by species
sample_infos = [
    sample for sample in SampleResolver.samples.values()
    if sample.schema.species.replace(" ", "_") == 
       self.model.species.replace(" ", "_")
]
```

### Prediction Output

Results are saved as TSV files:

```
data/SAMPLE_ID/predicted_celltype.tsv
```

#### Output Format

| Column | Type | Description |
|--------|------|-------------|
| `cell` | `str` | Unique cell identifier |
| `scpred_prediction` | `str` | Predicted cell type |
| `scpred_max` | `float` | Maximum prediction probability |
| `scpred_other` | `float` | Alternative prediction probability |

## Model Building Process

### Reference Data Requirements

For building custom models, you need:

1. **HDF5 Matrix**: Cell Ranger output or equivalent
2. **Cell Type Annotations**: TSV file with cell-celltype mapping

#### Annotation File Format

```tsv
cell	celltype
10X82_2_TCTCTCACCAGTTA	Astrocyte
10X82_2_TCTCTCACCAGTTC	Oligodendrocyte
10X82_2_TCTCTCACCAGTTT	Neuron
```

### Validation Checks

The build function validates input data:

```python
# Check file extensions
if not celltype_path.endswith(".tsv"):
    raise ValueError("celltype_path should be .tsv file")

# Validate column structure
df = pl.read_csv(celltype_path, separator="\t")
if df.columns != ["cell", "celltype"]:
    raise ValueError("Columns should be ['cell', 'celltype']")
```

### Model Training Pipeline

1. **Data Loading**: Load reference matrix and annotations
2. **Preprocessing**: Normalize and scale data
3. **Feature Selection**: Identify informative genes
4. **Model Training**: Train scPred classifier
5. **Validation**: Test model performance
6. **Export**: Save model files

## Performance Considerations

### Memory Requirements

| Dataset Size | Memory Needed |
|--------------|---------------|
| Model Building (50K cells) | 32-64 GB |
| Prediction (10K cells) | 8-16 GB |
| Prediction (50K cells) | 16-32 GB |

### Processing Time

| Operation | Dataset Size | Time |
|-----------|--------------|------|
| Model Building | 50K cells | 2-4 hours |
| Prediction | 10K cells | 5-10 minutes |
| Prediction | 50K cells | 15-30 minutes |

### Cluster Computing

For large-scale processing:

```bash
#!/bin/bash
#PBS -N PredictCelltype
#PBS -l nodes=1:ppn=8
#PBS -l mem=32gb
#PBS -l walltime=4:00:00

module load R/4.1.0
Rscript predict_celltype.R
```

## Quality Control

### Prediction Confidence

Monitor prediction quality:

```r
# Check prediction confidence
summary(predictions$scpred_max)

# Identify low-confidence predictions
low_conf <- predictions[predictions$scpred_max < 0.7, ]
```

### Common Cell Types

Typical cell types in different tissues:

**Brain/CNS**:
- Neuron, Astrocyte, Oligodendrocyte, Microglia
- OPC (Oligodendrocyte Precursor), Endothelial

**Immune System**:
- T_cell, B_cell, NK_cell, Macrophage, Dendritic_cell
- CD4_T, CD8_T, Regulatory_T, Plasma_cell

**Development**:
- Stem_cell, Progenitor_cell, Differentiated_cell

## Error Handling

### Common Issues

1. **Model Not Found**: Ensure model is built and registered
2. **Species Mismatch**: Check sample species matches model species
3. **Memory Exhaustion**: Reduce batch size or increase memory
4. **R Package Missing**: Install required R packages (scPred, Seurat)

### Troubleshooting

```python
# Check if model exists
model_dir = f"{Config.PROJ_ROOT}/reference/{species}/{suffix}"
if not os.path.exists(f"{model_dir}/reference.h5seurat"):
    raise FileNotFoundError("Model not found - build first")

# Validate samples ready for prediction
if not sample.path.is_counted:
    print(f"Sample {sample_id} not counted yet")
```

## Methods

### PredictCelltype Methods

#### `call(project: Project) -> Project`

Main execution method for cell type prediction.

#### `register() -> str`

Returns the function identifier for registration.

### BuildCellTypeModel Methods

#### `call(project: Project) -> Project`

Main execution method for model building.

#### `__show_help()`

Displays help information and example data format.

## Integration with Pipeline

### Complete Workflow

```python
from celline import Project
from celline.functions.count import Count
from celline.functions.predict_celltype import (
    BuildCellTypeModel, PredictCelltype, CellTypeModel
)
from celline.functions.preprocess import Preprocess

# Create project
project = Project("./my-project")

# Count reads
project.call(Count(nthread=8))

# Build custom model (optional)
project.call(BuildCellTypeModel(
    species="Homo sapiens",
    suffix="brain_atlas",
    nthread=8,
    h5matrix_path="/data/reference.h5",
    celltype_path="/data/celltypes.tsv"
))

# Predict cell types
model = CellTypeModel(species="Homo sapiens", suffix="brain_atlas")
project.call(PredictCelltype(model))

# Preprocess with predicted types
project.call(Preprocess(target_celltype=["Neuron", "Astrocyte"]))
```

## Model Management

### Available Models

List available models in your project:

```python
import os
from celline.config import Config

reference_dir = f"{Config.PROJ_ROOT}/reference"
for species in os.listdir(reference_dir):
    species_dir = os.path.join(reference_dir, species)
    if os.path.isdir(species_dir):
        print(f"Species: {species.replace('_', ' ')}")
        for model in os.listdir(species_dir):
            print(f"  Model: {model}")
```

### Model Versioning

Use suffixes for model versioning:

```python
# Different model versions
models = [
    CellTypeModel("Homo sapiens", "v1.0"),
    CellTypeModel("Homo sapiens", "v2.0_refined"),
    CellTypeModel("Homo sapiens", "custom_atlas")
]
```

## Related Functions

- [Count](count) - Generate count matrices before prediction
- [Preprocess](preprocess) - Quality control after prediction
- [CreateSeuratObject](createseuratobject) - Create Seurat objects with predictions
- [Integrate](integrate) - Integrate samples with predicted cell types

## Troubleshooting

### Common Issues

1. **R Dependencies**: Install scPred and Seurat packages
2. **Memory Issues**: Use cluster computing for large datasets
3. **Model Quality**: Validate reference data quality
4. **Prediction Accuracy**: Check model-data compatibility

### Debug Mode

Enable R debugging:

```r
# In R script
options(error = traceback)
debug(scPredict)
```

### Manual Prediction

For debugging, run prediction manually:

```bash
# Navigate to reference directory
cd reference/

# Execute prediction script
Rscript prediction.sh
```