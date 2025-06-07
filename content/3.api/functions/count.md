# Count Function

Process downloaded FASTQ files using Cell Ranger to generate feature-barcode matrices.

## Overview

The `Count` function processes downloaded raw sequencing data using Cell Ranger count. It generates feature-barcode matrices from FASTQ files, which are essential for downstream single-cell RNA-seq analysis.

## Class Information

- **Module**: `celline.functions.count`
- **Class**: `Count`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `nthread` | `int` | Yes | Number of threads to use for counting |
| `then` | `Optional[Callable[[str], None]]` | No | Callback function executed on successful completion |
| `catch` | `Optional[Callable[[subprocess.CalledProcessError], None]]` | No | Callback function executed on error |

### JobContainer Structure

The `Count.JobContainer` contains job configuration:

| Field | Type | Description |
|-------|------|-------------|
| `nthread` | `str` | Number of threads to use |
| `cluster_server` | `str` | Cluster server name (if applicable) |
| `jobname` | `str` | Job identifier |
| `logpath` | `str` | Path to log file |
| `sample_id` | `str` | Sample accession ID |
| `dist_dir` | `str` | Output directory |
| `fq_path` | `str` | Path to FASTQ files |
| `transcriptome` | `str` | Path to reference transcriptome |

## Usage Examples

### Python API

#### Basic Usage

```python
from celline import Project
from celline.functions.count import Count

# Create project
project = Project("./my-project")

# Create Count function instance
count_function = Count(nthread=8)

# Execute function
result = project.call(count_function)
```

#### With Callbacks

```python
from celline import Project
from celline.functions.count import Count
import subprocess

def on_success(sample_id: str):
    print(f"Successfully counted: {sample_id}")

def on_error(error: subprocess.CalledProcessError):
    print(f"Count failed: {error}")

# Create project
project = Project("./my-project")

# Create Count function with callbacks
count_function = Count(
    nthread=8,
    then=on_success,
    catch=on_error
)

# Execute function
result = project.call(count_function)
```

### CLI Usage

#### Basic Usage

```bash
# Count with default single thread
celline run count

# Count with multiple threads
celline run count --nthread 8

# Count with maximum performance
celline run count --nthread 16
```

### CLI Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--nthread`, `-n` | `int` | 1 | Number of threads to use for counting |

## Implementation Details

### Cell Ranger Integration

The function generates and executes Cell Ranger count commands:

```bash
cellranger count \
    --id=SAMPLE_ID \
    --transcriptome=/path/to/transcriptome \
    --fastqs=/path/to/fastqs \
    --sample=SAMPLE_ID \
    --localcores=8
```

### Processing Workflow

1. **Sample Discovery**: Reads `samples.toml` to identify samples
2. **Path Validation**: Ensures FASTQ files exist
3. **Transcriptome Resolution**: Finds appropriate reference transcriptome
4. **Script Generation**: Creates Cell Ranger execution scripts
5. **Job Execution**: Runs counting jobs with thread management
6. **Output Validation**: Verifies successful completion

### File Organization

Output files are organized as follows:

```
project_root/
├── resources/
│   └── SAMPLE_ID/
│       ├── raw/
│       │   └── fastqs/           # Input FASTQ files
│       ├── counted/
│       │   └── outs/
│       │       ├── filtered_feature_bc_matrix.h5
│       │       ├── filtered_feature_bc_matrix/
│       │       ├── raw_feature_bc_matrix.h5
│       │       ├── raw_feature_bc_matrix/
│       │       ├── analysis/
│       │       ├── cloupe.cloupe
│       │       ├── metrics_summary.csv
│       │       └── web_summary.html
│       ├── src/
│       │   └── count.sh          # Generated count script
│       └── log/
│           └── count_*.log       # Count logs
```

### Cell Ranger Outputs

The function generates standard Cell Ranger outputs:

| File | Description |
|------|-------------|
| `filtered_feature_bc_matrix.h5` | Filtered feature-barcode matrix in HDF5 format |
| `filtered_feature_bc_matrix/` | Filtered matrix in Matrix Market format |
| `raw_feature_bc_matrix.h5` | Raw feature-barcode matrix in HDF5 format |
| `raw_feature_bc_matrix/` | Raw matrix in Matrix Market format |
| `analysis/` | Secondary analysis results (PCA, t-SNE, UMAP, clustering) |
| `cloupe.cloupe` | Loupe Cell Browser file |
| `metrics_summary.csv` | Run statistics and metrics |
| `web_summary.html` | Summary report |

## Transcriptome Requirements

### Reference Preparation

Before running count, ensure transcriptomes are registered:

```python
from celline.DB.model.transcriptome import Transcriptome

# Add human reference
Transcriptome().add_path(
    species="Homo sapiens",
    built_path="/path/to/refdata-gex-GRCh38-2020-A"
)

# Add mouse reference
Transcriptome().add_path(
    species="Mus musculus", 
    built_path="/path/to/refdata-gex-mm10-2020-A"
)
```

### Supported Species

The function automatically resolves transcriptomes for:
- Homo sapiens (human)
- Mus musculus (mouse)
- Custom species (user-defined)

## Performance Optimization

### Thread Configuration

Optimal thread counts by system type:

| System Type | Recommended Threads |
|-------------|-------------------|
| Laptop/Desktop | 4-8 |
| Workstation | 8-16 |
| HPC Node | 16-32 |
| Cloud Instance | Instance-dependent |

### Memory Requirements

Cell Ranger memory requirements:

| Dataset Size | Memory Needed |
|--------------|---------------|
| <10K cells | 32 GB |
| 10K-50K cells | 64 GB |
| 50K-100K cells | 128 GB |
| >100K cells | 256+ GB |

### Storage Requirements

Typical storage needs:

| Input Size | Output Size |
|------------|-------------|
| 1 GB FASTQ | 2-3 GB |
| 5 GB FASTQ | 10-15 GB |
| 10 GB FASTQ | 20-30 GB |

## Error Handling

### Common Issues

1. **Insufficient Memory**: Cell Ranger requires substantial RAM
2. **Missing Transcriptome**: Reference must be registered first
3. **Invalid FASTQ**: Files must be properly formatted
4. **Disk Space**: Ensure adequate storage for outputs

### Recovery Strategies

```python
# Retry with reduced threads if memory issues
try:
    count_function = Count(nthread=16)
    project.call(count_function)
except MemoryError:
    count_function = Count(nthread=8)
    project.call(count_function)
```

## Cluster Computing

### PBS/Slurm Integration

For cluster environments:

```python
from celline.server import ServerSystem

# Configure cluster settings
ServerSystem.job_system = ServerSystem.JobType.PBS
ServerSystem.cluster_server_name = "my-cluster"

# Count jobs will be submitted to cluster
project.call(Count(nthread=16))
```

### Job Templates

The function uses customizable job templates:

```bash
#!/bin/bash
#PBS -N Count
#PBS -l nodes=1:ppn=16
#PBS -l mem=64gb
#PBS -l walltime=24:00:00

cd $PBS_O_WORKDIR
cellranger count --id=SAMPLE_ID ...
```

## Methods

### `call(project: Project) -> Project`

Main execution method that orchestrates the counting process.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

### `add_cli_args(parser: ArgumentParser) -> None`

Adds CLI-specific arguments to the argument parser.

### `cli(project: Project, args: Namespace) -> Project`

CLI entry point that processes command-line arguments and executes the function.

## Quality Control

### Metrics Validation

The function automatically validates key metrics:

```python
# Example metrics from Cell Ranger output
metrics = {
    "estimated_number_of_cells": 5000,
    "mean_reads_per_cell": 50000,
    "median_genes_per_cell": 2500,
    "fraction_reads_in_cells": 0.85,
    "valid_barcodes": 0.98
}
```

### Automatic Warnings

The function warns about potential issues:
- Low cell recovery rates
- High mitochondrial gene expression
- Poor sequencing quality
- Insufficient read depth

## Integration with Pipeline

### Typical Workflow

```python
from celline import Project
from celline.functions.add import Add
from celline.functions.download import Download
from celline.functions.count import Count
from celline.functions.preprocess import Preprocess

# Complete pipeline
project = Project("./my-project")

# Add samples
project.call(Add([Add.SampleInfo(id="GSE123456")]))

# Download data
project.call(Download())

# Count reads
project.call(Count(nthread=8))

# Preprocess results
project.call(Preprocess())
```

## Related Functions

- [Download](download) - Download FASTQ files before counting
- [SetTranscriptome](settranscriptome) - Register reference transcriptomes
- [Preprocess](preprocess) - Quality control after counting
- [CreateSeuratObject](createseuratobject) - Create Seurat objects from count data

## Troubleshooting

### Common Issues

1. **Cell Ranger Not Found**: Ensure Cell Ranger is installed and in PATH
2. **Permission Errors**: Check write permissions in output directories
3. **Transcriptome Missing**: Register required reference genomes
4. **Memory Exhaustion**: Reduce thread count or increase system memory

### Debug Mode

Enable detailed logging:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

### Manual Execution

For debugging, run Cell Ranger manually:

```bash
# Navigate to sample directory
cd resources/SAMPLE_ID/src/

# Execute count script manually
bash count.sh
```