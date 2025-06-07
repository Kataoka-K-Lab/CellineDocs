# Functions Reference

This is a detailed reference for each analysis function in Celline.
It explains the parameters, usage examples, and implementation details for all functions.

## 📚 Function Categories

### 🗄️ Data Management
- [Add](#add) - Sample addition and metadata acquisition
- [Download](#download) - Sequencing data download
- [SyncDB](#syncdb) - Database synchronization

### 🧮 Data Processing
- [Count](#count) - Count processing with Cell Ranger
- [CreateSeuratObject](#createseuratobject) - Seurat object creation
- [SetTranscriptome](#settranscriptome) - Transcriptome reference database setup

### 🔬 Quality Control & Preprocessing
- [Preprocess](#preprocess) - Quality control and filtering
- [VCount](#vcount) - Count data verification

### 📊 Analysis & Visualization
- [PredictCelltype](#predictcelltype) - Cell type prediction
- [Reduce](#reduce) - Dimensionality reduction
- [Integrate](#integrate) - Data integration and batch effect correction
- [BatchCor](#batchcor) - Batch correlation analysis

### 🔧 Utilities
- [Info](#info) - System information display
- [Initialize](#initialize) - Initialization and setup
- [Interactive](#interactive) - Interactive mode
- [Job](#job) - Job management
- [Bash](#bash) - Shell command execution

---

## 📊 Data Management Functions

### Add

Retrieves sample IDs from public databases and adds them to the project.

#### Class Definition

```python
class Add(CellineFunction):
    def __init__(self, sample_id: Union[List[SampleInfo], pl.DataFrame]) -> None
```

#### Parameters

**SampleInfo**
```python
class SampleInfo(NamedTuple):
    id: str                    # Sample ID (GSE, GSM, SRR, etc.)
    title: Optional[str] = ""  # Sample title
```

#### Usage Examples

```python
from celline.functions.add import Add
import polars as pl

# Using SampleInfo list
samples = [
    Add.SampleInfo(id="GSE123456", title="Brain tissue study"),
    Add.SampleInfo(id="GSM789012", title="Control sample"),
    Add.SampleInfo(id="GSM789013", title="Treatment sample")
]
add_func = Add(samples)
project.call(add_func)

# Using DataFrame
df = pl.DataFrame({
    "id": ["GSE123456", "GSM789012"],
    "title": ["Study 1", "Sample 1"]
})
add_func = Add(df)
project.call(add_func)

# CLI usage
# celline run add GSE123456
# celline run add GSM789012 GSM789013 --title "My samples"
# celline run add --from-file samples.csv
```

#### Methods

##### `get_samples() -> Dict[str, str]`

Retrieves sample information for the current project.

```python
add_func = Add([])
samples = add_func.get_samples()
print(samples)  # {'GSM123456': 'Sample 1', 'GSM123457': 'Sample 2'}
```

#### Supported Databases

- **GEO (Gene Expression Omnibus)**
  - GSE (Studies)
  - GSM (Samples)
- **SRA (Sequence Read Archive)**
  - SRR (Runs)
- **CNCB (China National Center for Bioinformation)**
  - PRJCA (Projects)
  - CRA (Study)
  - CRR (Runs)

#### Error Handling

```python
try:
    add_func = Add([Add.SampleInfo(id="INVALID_ID")])
    project.call(add_func)
except ValueError as e:
    print(f"Invalid sample ID format: {e}")
except ConnectionError as e:
    print(f"Database connection failed: {e}")
```

---

### Download

Downloads sequencing data for added samples.

#### Class Definition

```python
class Download(CellineFunction):
    def __init__(
        self,
        then: Optional[Callable[[str], None]] = None,
        catch: Optional[Callable[[subprocess.CalledProcessError], None]] = None,
    ) -> None
```

#### Parameters

- `then`: Callback function for download completion
- `catch`: Callback function for error handling

#### Usage Examples

```python
from celline.functions.download import Download

# Basic download
download = Download()
project.call(download)

# Download with callbacks
def on_complete(sample_id):
    print(f"✓ Downloaded: {sample_id}")

def on_error(error):
    print(f"✗ Error downloading: {error}")

download = Download(then=on_complete, catch=on_error)
project.call(download)

# CLI usage
# celline run download
# celline run download --nthread 4
# celline run download --force
```

#### Internal Data Structure

##### `JobContainer`

```python
class JobContainer(NamedTuple):
    filetype: str           # Data format (FASTQ, SRA, etc.)
    nthread: str           # Number of threads
    cluster_server: str    # Cluster server name
    jobname: str          # Job name
    logpath: str          # Log file path
    sample_id: str        # Sample ID
    download_target: str  # Download target directory
    download_source: str  # Download source URL
    run_ids_str: str      # Run ID string
```

#### Download Targets

- **FASTQ files** - Illumina/10x Genomics data
- **SRA files** - Download from SRA archive
- **BAM files** - Aligned data

#### Progress Monitoring

```python
# Check progress with log files
import os
from celline.utils.path import Path

def monitor_download_progress(sample_id):
    path = Path("GSE123456", sample_id)
    log_files = os.listdir(path.resources_sample_log)
    
    for log_file in log_files:
        if "download" in log_file:
            with open(f"{path.resources_sample_log}/{log_file}") as f:
                print(f.read())
```

---

### SyncDB

Synchronizes the local database with the latest public databases.

#### Class Definition

```python
class SyncDB(CellineFunction):
    def __init__(self) -> None
```

#### Usage Examples

```python
from celline.functions.sync_DB import SyncDB

# Database synchronization
sync = SyncDB()
project.call(sync)

# CLI usage
# celline run sync_DB
```

#### Synchronization Targets

- **GEO metadata** - GSE/GSM information
- **SRA metadata** - SRR information  
- **Transcriptome references** - Reference transcriptome information
- **Gene annotations** - Gene annotations

---

## 🧮 Data Processing Functions

### Count

Generates expression count matrices from FASTQ files using Cell Ranger.

#### Class Definition

```python
class Count(CellineFunction):
    def __init__(self) -> None
```

#### Usage Examples

```python
from celline.functions.count import Count

# Count processing
count = Count()
project.call(count)

# CLI usage
# celline run count
```

#### Output Files

- `filtered_feature_bc_matrix.h5` - Filtered expression matrix
- `raw_feature_bc_matrix.h5` - Raw expression matrix
- `metrics_summary.csv` - Quality metrics
- `web_summary.html` - HTML report

#### Cell Ranger Configuration

```python
# Check Cell Ranger parameters
import subprocess

# Check Cell Ranger version
result = subprocess.run(['cellranger', '--version'], 
                       capture_output=True, text=True)
print(f"Cell Ranger version: {result.stdout}")

# Available commands
commands = ['count', 'aggr', 'reanalyze', 'mkfastq']
```

---

### CreateSeuratObject

Creates Seurat objects from count data and performs basic analysis.

#### Class Definition

```python
class CreateSeuratObject(CellineFunction):
    def __init__(self) -> None
```

#### Usage Examples

```python
from celline.functions.create_seurat import CreateSeuratObject

# Seurat object creation
create_seurat = CreateSeuratObject()
project.call(create_seurat)

# CLI usage
# celline run create_seurat
```

#### R Code Example (Internal Execution)

```r
# R code executed internally by Celline
library(Seurat)
library(SeuratDisk)

# Data loading
raw <- Read10X_h5(h5_path)

# Seurat object creation
seurat_obj <- CreateSeuratObject(raw, project = proj) %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(features = VariableFeatures(.)) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(dims = 1:20) %>%
    RunUMAP(dims = 1:20)

# Save
saveRDS(seurat_obj, h5seurat_path)
```

#### Generated Analysis Results

- **Normalization**: LogNormalize
- **Highly Variable Genes**: FindVariableFeatures (top 2000)
- **Scaling**: ScaleData
- **Principal Component Analysis**: PCA (50 components)
- **Neighbor Graph**: FindNeighbors
- **Clustering**: Leiden algorithm
- **UMAP**: 2D embedding

---

### SetTranscriptome

Sets up the transcriptome reference database for analysis.

#### Class Definition

```python
class SetTranscriptome(CellineFunction):
    def __init__(self, species: Optional[str] = None, version: Optional[str] = None) -> None
```

#### Parameters

- `species`: Target species (homo_sapiens, mus_musculus, etc.)
- `version`: Genome version (GRCh38, GRCm39, etc.)

#### Usage Examples

```python
from celline.functions.set_transcriptome import SetTranscriptome

# Auto-detection
set_ref = SetTranscriptome()
project.call(set_ref)

# Manual specification
set_ref = SetTranscriptome(species="homo_sapiens", version="GRCh38")
project.call(set_ref)

# CLI usage
# celline run set_transcriptome
# celline run set_transcriptome --species homo_sapiens --version GRCh38
```

#### Supported Species & Versions

| Species | Available Versions |
|---------|-------------------|
| homo_sapiens | GRCh38, GRCh37 |
| mus_musculus | GRCm39, GRCm38 |
| rattus_norvegicus | Rnor6.0 |
| danio_rerio | GRCz11 |

---

## 🔬 Quality Control & Preprocessing Functions

### Preprocess

Performs quality control, filtering, and doublet detection.

#### Class Definition

```python
class Preprocess(CellineFunction):
    def __init__(self, target_celltype: Optional[List[str]] = None) -> None
```

#### Parameters

- `target_celltype`: List of cell types to filter for

#### Usage Examples

```python
from celline.functions.preprocess import Preprocess

# Basic preprocessing
preprocess = Preprocess()
project.call(preprocess)

# Target specific cell types only
preprocess = Preprocess(target_celltype=["Neuron", "Astrocyte", "Oligodendrocyte"])
project.call(preprocess)

# CLI usage
# celline run preprocess
# celline run preprocess --target-celltype Neuron Astrocyte
```

#### QC Metrics

```python
# Generated QC metrics
qc_metrics = {
    "n_genes_by_counts": "Number of detected genes",
    "total_counts": "Total read count", 
    "pct_counts_mt": "Mitochondrial gene percentage",
    "doublet_score": "Doublet score",
    "predicted_doublets": "Doublet prediction"
}
```

#### Filtering Criteria

| Metric | Lower Limit | Upper Limit |
|---------|-------------|-------------|
| n_genes_by_counts | 200 | 5000 |
| pct_counts_mt | - | 20% |
| doublet_score | - | 0.3 |

#### Output Files

- `cell_info.tsv` - Cell-level QC information
- `qc_plots.pdf` - QC plots
- `filtering_report.txt` - Filtering report

---

### VCount

Validates the integrity of count data.

#### Class Definition

```python
class VCount(CellineFunction):
    def __init__(self) -> None
```

#### Usage Examples

```python
from celline.functions.vcount import VCount

# Count data validation
vcount = VCount()
project.call(vcount)
```

#### Validation Items

- File existence check
- Matrix dimension check
- Count distribution verification
- Quality metrics validation

---

## 📊 Analysis & Visualization Functions

### PredictCelltype

Predicts cell types using machine learning models.

#### Class Definition

```python
class PredictCelltype(CellineFunction):
    def __init__(self, model_path: Optional[str] = None) -> None
```

#### Parameters

- `model_path`: Path to custom trained model

#### Usage Examples

```python
from celline.functions.predict_celltype import PredictCelltype

# Use default model
predict = PredictCelltype()
project.call(predict)

# Use custom model
predict = PredictCelltype(model_path="/path/to/custom_model.pkl")
project.call(predict)

# CLI usage
# celline run predict_celltype
# celline run predict_celltype --model-path /path/to/model.pkl
```

#### Model Information

##### Pre-trained Models

- **Human Brain**: Human brain cell type classification (Allen Brain Atlas)
- **Mouse Brain**: Mouse brain cell type classification (Allen Brain Atlas)
- **PBMC**: Peripheral blood mononuclear cell classification
- **Universal**: General-purpose cell type classification

##### Custom Model Creation

```python
from celline.functions.predict_celltype import BuildCellTypeModel

# Train custom model
model_builder = BuildCellTypeModel(
    training_data="/path/to/training.h5ad",
    features=["CD4", "CD8A", "CD19", "CD14"],
    model_type="random_forest"
)
project.call(model_builder)
```

---

### Reduce

Performs dimensionality reduction using PCA, UMAP, and t-SNE.

#### Class Definition

```python
class Reduce(CellineFunction):
    def __init__(
        self, 
        methods: List[str] = ["pca", "umap"],
        n_components: int = 50,
        n_neighbors: int = 15
    ) -> None
```

#### Parameters

- `methods`: Dimensionality reduction methods to use
- `n_components`: Number of principal components
- `n_neighbors`: Number of neighbors (UMAP)

#### Usage Examples

```python
from celline.functions.reduce import Reduce

# Default settings
reduce = Reduce()
project.call(reduce)

# Custom settings
reduce = Reduce(
    methods=["pca", "umap", "tsne"],
    n_components=30,
    n_neighbors=20
)
project.call(reduce)

# CLI usage
# celline run reduce
# celline run reduce --methods pca umap tsne --n-components 30
```

#### Available Methods

| Method | Description | Parameters |
|--------|-------------|------------|
| PCA | Principal Component Analysis | n_components |
| UMAP | Uniform Manifold Approximation | n_neighbors, min_dist |
| t-SNE | t-distributed Stochastic Neighbor Embedding | perplexity |

---

### Integrate

Performs multi-sample integration and batch effect correction.

#### Class Definition

```python
class Integrate(CellineFunction):
    def __init__(
        self,
        method: str = "harmony",
        batch_key: str = "sample",
        n_components: int = 50
    ) -> None
```

#### Parameters

- `method`: Integration method (harmony, combat, mnn, cca)
- `batch_key`: Column name for batch information
- `n_components`: Number of principal components for integration

#### Usage Examples

```python
from celline.functions.integrate import Integrate

# Use Harmony
integrate = Integrate(method="harmony")
project.call(integrate)

# Use MNN
integrate = Integrate(
    method="mnn",
    batch_key="dataset",
    n_components=30
)
project.call(integrate)

# CLI usage
# celline run integrate
# celline run integrate --method mnn --batch-key dataset
```

#### Integration Method Comparison

| Method | Features | Use Case |
|--------|----------|----------|
| Harmony | Fast, effective | Large-scale datasets |
| Combat | Traditional method | Small to medium datasets |
| MNN | Preserves biological variation | Different protocols |
| CCA | Seurat standard | 10x Genomics data |

---

### BatchCor

Analyzes batch correlations and evaluates batch effects.

#### Class Definition

```python
class BatchCor(CellineFunction):
    def __init__(self, batch_keys: List[str] = ["sample", "dataset"]) -> None
```

#### Usage Examples

```python
from celline.functions.batch_cor import BatchCor

# Batch correlation analysis
batch_cor = BatchCor(batch_keys=["sample", "condition", "batch"])
project.call(batch_cor)
```

---

## 🔧 Utility Functions

### Info

Displays system information and project status.

#### Class Definition

```python
class Info(CellineFunction):
    def __init__(self) -> None
```

#### Usage Examples

```python
from celline.functions.info import Info

# Display system information
info = Info()
project.call(info)

# CLI usage
# celline info
```

### Initialize

Performs project initialization and setup.

#### Class Definition

```python
class Initialize(CellineFunction):
    def __init__(self) -> None
```

#### Usage Examples

```python
from celline.functions.initialize import Initialize

# Initialization
init = Initialize()
project.call(init)

# CLI usage
# celline init
```

### Interactive

Launches the interactive web interface.

#### Class Definition

```python
class Interactive(CellineFunction):
    def __init__(self, port: int = 8080) -> None
```

#### Usage Examples

```python
from celline.functions.interactive import Interactive

# Launch with default port
interactive = Interactive()
project.call(interactive)

# Launch with custom port
interactive = Interactive(port=8090)
project.call(interactive)

# CLI usage
# celline interactive
# celline run interactive --port 8090
```

### Job

Manages and monitors jobs.

#### Class Definition

```python
class Job(CellineFunction):
    def __init__(self, action: str = "status") -> None
```

#### Usage Examples

```python
from celline.functions.job import Job

# Check job status
job = Job(action="status")
project.call(job)

# Cancel job
job = Job(action="cancel")
project.call(job)
```

### Bash

Executes shell commands.

#### Class Definition

```python
class Bash(CellineFunction):
    def __init__(self, command: str, timeout: Optional[int] = None) -> None
```

#### Usage Examples

```python
from celline.functions.bash import Bash

# Execute shell command
bash = Bash("ls -la resources/")
project.call(bash)

# Execute with timeout
bash = Bash("long_running_command.sh", timeout=3600)
project.call(bash)
```

---

> **Info**: For detailed implementation examples and customization methods for each function, refer to the [Developer Guide](/celline/development).