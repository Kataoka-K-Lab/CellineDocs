---
title: API Reference
navigation: false      # ← 見出し専用にしてリンクは非表示
---

This is the complete reference for Celline's Python API.
It provides detailed explanations of classes, functions, and methods for calling Celline's functionality directly from programs.

## 🎯 API Overview

Celline's Python API consists of the following main components:

- **Project**: Project management and workflow control
- **Functions**: Analysis functionality implementation
- **Data Handlers**: Data types and I/O operations
- **Database**: Metadata and sample information management
- **Visualization**: Plotting and interactive features

## 🏗️ Core Classes

### Project

The main class that manages projects and controls analysis workflows.

```python
from celline import Project

# Project creation
project = Project(project_dir="./my-project", proj_name="example")

# Basic usage
project.call(function_instance)
project.parallelize(4)
project.singularize()
```

#### Constructor

```python
def __init__(
    self, 
    project_dir: str, 
    proj_name: str = "", 
    r_path: str = ""
) -> None:
```

**Parameters:**
- `project_dir` (str): Path to the project directory
- `proj_name` (str, optional): Project name (default: directory name)
- `r_path` (str, optional): R path (default: auto-detection)

#### Main Methods

##### `call(func, wait_for_complete=True)`

Executes Celline functions.

```python
from celline.functions.add import Add
from celline.functions.download import Download

# Sample addition
add_func = Add([Add.SampleInfo(id="GSE123456", title="Sample 1")])
project.call(add_func)

# Data download
download_func = Download()
project.call(download_func)
```

##### `parallelize(njobs)` / `singularize()`

Controls parallel execution.

```python
# Execute with 4 parallel processes
project.parallelize(4)

# Return to single thread
project.singularize()
```

##### `call_if_else(condition, true, false)`

Controls execution with conditional branching.

```python
project.call_if_else(
    condition=lambda p: len(p.get_samples()) > 10,
    true=HighThroughputProcessing(),
    false=StandardProcessing()
)
```

##### Execution Environment Control

```python
# Multithreading execution
project.useMultiThreading()

# PBS cluster execution
project.usePBS("cluster-name")
```

### Seurat Data Objects

```python
# Obtaining Seurat objects
seurat = project.seurat(
    project_id="GSE123456",
    sample_id="GSM789012",
    identifier="seurat.seurat",
    via_seurat_disk=False
)

# Direct loading from file path
seurat = project.seurat_from_rawpath("/path/to/seurat.rds")
```

## 🔬 Analysis Function Classes

### Base Class: CellineFunction

The base class for all analysis functions.

```python
from celline.functions._base import CellineFunction
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from celline import Project

class CustomFunction(CellineFunction):
    def __init__(self, param1: str, param2: int):
        super().__init__()
        self.param1 = param1
        self.param2 = param2
    
    def call(self, project: "Project"):
        # Custom processing implementation
        print(f"Processing with {self.param1}, {self.param2}")
        return project
    
    def get_description(self) -> str:
        return "Custom analysis function"
```

### Main Analysis Functions

#### Add - Sample Addition

```python
from celline.functions.add import Add
import polars as pl

# Using SampleInfo
samples = [
    Add.SampleInfo(id="GSE123456", title="Dataset 1"),
    Add.SampleInfo(id="GSM789012", title="Sample 1")
]
add_func = Add(samples)

# Using DataFrame
df = pl.DataFrame({
    "id": ["GSE123456", "GSM789012"],
    "title": ["Dataset 1", "Sample 1"]
})
add_func = Add(df)

project.call(add_func)
```

#### Download - Data Download

```python
from celline.functions.download import Download

# Basic download
download = Download()
project.call(download)

# Download with callbacks
def on_complete(sample_id):
    print(f"Downloaded: {sample_id}")

def on_error(error):
    print(f"Error: {error}")

download = Download(then=on_complete, catch=on_error)
project.call(download)
```

#### Preprocess - Preprocessing

```python
from celline.functions.preprocess import Preprocess

# Basic preprocessing
preprocess = Preprocess()
project.call(preprocess)

# Target specific cell types
preprocess = Preprocess(target_celltype=["Neuron", "Astrocyte"])
project.call(preprocess)
```

#### Other Main Functions

```python
from celline.functions import *

# Count processing
project.call(Count())

# Seurat object creation
project.call(CreateSeuratObject())

# Cell type prediction
project.call(PredictCelltype())

# Dimensionality reduction
project.call(Reduce())

# Data integration
project.call(Integrate())
```

## 🗄️ Database API

### Database Handlers

```python
from celline.DB.dev.handler import HandleResolver
from celline.DB.dev.model import SampleSchema, RunSchema

# Sample ID resolution
resolver = HandleResolver.resolve("GSM123456")
if resolver:
    sample_schema = resolver.sample.search("GSM123456")
    print(f"Sample title: {sample_schema.title}")
    print(f"Parent study: {sample_schema.parent}")
```

### Database Models

#### SampleSchema

```python
from celline.DB.model.sra_gsm import SRA_GSM

# GSM sample search
gsm_model = SRA_GSM()
sample = gsm_model.search("GSM123456")

print(f"Title: {sample.title}")
print(f"Organism: {sample.organism}")
print(f"Library strategy: {sample.library_strategy}")
```

#### Database Synchronization

```python
from celline.functions.sync_DB import SyncDB

# Database synchronization
sync = SyncDB()
project.call(sync)
```

## 📊 Data Handling

### Seurat Data Operations

```python
from celline.data import Seurat

# Seurat object loading
seurat = Seurat("/path/to/seurat.rds", via_seurat_disk=True)

# Data retrieval
expression_matrix = seurat.get_expression()
metadata = seurat.get_metadata()
variable_features = seurat.get_variable_features()
```

### Data Format Conversion

```python
from celline.utils.serialization import NamedTupleAndPolarsStructure
import polars as pl

# NamedTuple and DataFrame interconversion
data_converter = NamedTupleAndPolarsStructure[Add.SampleInfo]

# DataFrame to NamedTuple list
df = pl.DataFrame({"id": ["GSM1", "GSM2"], "title": ["S1", "S2"]})
sample_list = data_converter.deserialize(df, Add.SampleInfo)

# NamedTuple list to DataFrame
samples = [Add.SampleInfo("GSM1", "S1"), Add.SampleInfo("GSM2", "S2")]
df = data_converter.serialize(samples)
```

## 🔧 Utilities

### Path Management

```python
from celline.utils.path import Path

# Project path management
path = Path("GSE123456", "GSM789012")

# Directory preparation
path.prepare()

# Various path retrieval
print(f"Raw data: {path.resources_sample_raw}")
print(f"Counted data: {path.resources_sample_counted}")
print(f"Results: {path.data_sample}")

# Status checking
if path.is_downloaded:
    print("Data has been downloaded")
if path.is_counted:
    print("Data has been counted")
```

### Configuration Management

```python
from celline.config import Config, Setting

# Project configuration
print(f"Project root: {Config.PROJ_ROOT}")
print(f"Execution root: {Config.EXEC_ROOT}")

# Execution settings
print(f"Project name: {Setting.name}")
print(f"R path: {Setting.r_path}")
print(f"Thread count: {Setting.nthread}")
print(f"Execution system: {Setting.system}")

# Configuration saving
Setting.name = "new-project"
Setting.nthread = 8
Setting.flush()  # Save to setting.toml
```

## 🎭 Interactive API

### Web API Server

```python
from celline.api.main import app
import uvicorn

# FastAPI application startup
if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
```

### API Endpoints

```python
import requests

# Get project information
response = requests.get("http://localhost:8000/api/project")
project_info = response.json()

# Add samples
add_request = {
    "sample_ids": ["GSM123456", "GSM789012"]
}
response = requests.post("http://localhost:8000/api/samples/add", 
                        json=add_request)
job_info = response.json()

# Check job status
job_id = job_info["job_id"]
response = requests.get(f"http://localhost:8000/api/jobs/{job_id}")
status = response.json()
```

## 🧵 Parallel Processing and Thread Management

### ThreadObservable

```python
from celline.middleware import ThreadObservable

# Parallel job configuration
ThreadObservable.set_jobs(4)

# Parallel execution of shell commands
shell_files = [
    "script1.sh",
    "script2.sh", 
    "script3.sh"
]
ThreadObservable.call_shell(shell_files)

# Wait for completion
ThreadObservable.watch()
```

### Server System Management

```python
from celline.server import ServerSystem

# Multithreading execution
ServerSystem.useMultiThreading()

# PBS cluster execution
ServerSystem.usePBS("my-cluster")

# Check current configuration
print(f"Current system: {ServerSystem.current_system}")
print(f"Cluster server: {ServerSystem.cluster_server_name}")
```

## 🔍 Logging and Debugging

### Logging

```python
from celline.log.logger import get_logger

# Get logger
logger = get_logger(__name__)

# Log output
logger.info("Processing started")
logger.warning("Low memory warning")
logger.error("Processing failed")

# Contextual logging
logger.info("Processing sample", extra={
    "sample_id": "GSM123456",
    "step": "preprocessing"
})
```

## 📝 Custom Function Creation Example

### Complete Custom Function

```python
from celline.functions._base import CellineFunction
from celline.log.logger import get_logger
import argparse
from typing import Optional, List

class AdvancedQualityControl(CellineFunction):
    """Advanced quality control function"""
    
    def __init__(self, 
                 min_genes: int = 200,
                 max_genes: int = 5000,
                 max_mito_pct: float = 20.0,
                 doublet_threshold: float = 0.3):
        super().__init__()
        self.min_genes = min_genes
        self.max_genes = max_genes
        self.max_mito_pct = max_mito_pct
        self.doublet_threshold = doublet_threshold
        self.logger = get_logger(__name__)
    
    def call(self, project):
        """Main processing"""
        self.logger.info("Starting advanced quality control")
        
        samples = self.get_samples_from_project(project)
        
        for sample_id in samples:
            self.logger.info(f"Processing sample: {sample_id}")
            self.process_sample(project, sample_id)
        
        self.logger.info("Advanced quality control completed")
        return project
    
    def process_sample(self, project, sample_id):
        """Individual sample processing"""
        import scanpy as sc
        import pandas as pd
        
        # Data loading
        adata = self.load_sample_data(project, sample_id)
        
        # Quality metrics calculation
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, percent_top=None, 
                                  log1p=False, inplace=True)
        
        # Filtering
        sc.pp.filter_cells(adata, min_genes=self.min_genes)
        sc.pp.filter_genes(adata, min_cells=3)
        
        # Remove outlier cells
        adata = adata[adata.obs.n_genes_by_counts < self.max_genes, :]
        adata = adata[adata.obs.pct_counts_mt < self.max_mito_pct, :]
        
        # Doublet detection
        doublet_scores = self.detect_doublets(adata)
        adata.obs['doublet_score'] = doublet_scores
        adata = adata[doublet_scores < self.doublet_threshold, :]
        
        # Save results
        self.save_sample_data(project, sample_id, adata)
        
        # Generate QC report
        self.generate_qc_report(project, sample_id, adata)
    
    def detect_doublets(self, adata):
        """Doublet detection"""
        import scrublet as scr
        scrub = scr.Scrublet(adata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        return doublet_scores
    
    def add_cli_args(self, parser: argparse.ArgumentParser):
        """Add CLI arguments"""
        parser.add_argument('--min-genes', type=int, default=200,
                          help='Minimum number of genes per cell')
        parser.add_argument('--max-genes', type=int, default=5000,
                          help='Maximum number of genes per cell')
        parser.add_argument('--max-mito-pct', type=float, default=20.0,
                          help='Maximum mitochondrial gene percentage')
        parser.add_argument('--doublet-threshold', type=float, default=0.3,
                          help='Doublet score threshold')
    
    def cli(self, project, args: Optional[argparse.Namespace] = None):
        """CLI entry point"""
        if args:
            self.min_genes = args.min_genes
            self.max_genes = args.max_genes
            self.max_mito_pct = args.max_mito_pct
            self.doublet_threshold = args.doublet_threshold
        
        return self.call(project)
    
    def get_description(self) -> str:
        return """Advanced quality control with customizable thresholds.
        
        Performs comprehensive QC including gene count filtering,
        mitochondrial gene percentage filtering, and doublet detection."""
    
    def get_usage_examples(self) -> List[str]:
        return [
            "celline run advanced_qc",
            "celline run advanced_qc --min-genes 300 --max-genes 6000",
            "celline run advanced_qc --max-mito-pct 15 --doublet-threshold 0.25"
        ]

# Usage example
qc = AdvancedQualityControl(min_genes=300, max_genes=6000)
project.call(qc)
```

---

> **Info**: For detailed usage examples of the API, see the [Functions Reference](/3.api/functions) section for detailed explanations of individual functions.