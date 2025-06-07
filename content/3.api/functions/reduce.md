# Reduce Function

Clean up and reduce storage space by removing unnecessary Cell Ranger output files.

## Overview

The `Reduce` function optimizes storage usage by removing large, redundant files from Cell Ranger output while preserving essential data files needed for downstream analysis. This is particularly useful for large projects with many samples.

## Class Information

- **Module**: `celline.functions.reduce`
- **Class**: `Reduce`
- **Base Class**: `CellineFunction`

## Parameters

The `Reduce` function takes no constructor parameters.

## Usage Examples

### Python API

#### Basic Usage

```python
from celline import Project
from celline.functions.reduce import Reduce

# Create project
project = Project("./my-project")

# Create Reduce function instance
reduce_function = Reduce()

# Execute function
result = project.call(reduce_function)
```

#### Pipeline Integration

```python
from celline import Project
from celline.functions.count import Count
from celline.functions.reduce import Reduce
from celline.functions.preprocess import Preprocess

# Create project
project = Project("./my-project")

# Complete counting
project.call(Count(nthread=8))

# Reduce storage usage
project.call(Reduce())

# Continue with preprocessing
project.call(Preprocess())
```

### CLI Usage

#### Basic Usage

```bash
# Reduce storage for all samples
celline run reduce

# Typically used after counting
celline run count --nthread 8
celline run reduce
```

## Implementation Details

### Files Preserved

The function preserves essential Cell Ranger output files:

| File/Directory | Description | Size Impact |
|----------------|-------------|-------------|
| `outs/filtered_feature_bc_matrix.h5` | Filtered count matrix (HDF5) | Essential |
| `outs/molecule_info.h5` | Molecule-level information | Large |
| `outs/web_summary.html` | Summary report | Small |
| `outs/filtered_feature_bc_matrix/` | Filtered matrix (Matrix Market) | Medium |
| `├── matrix.mtx.gz` | Count matrix | |
| `├── features.tsv.gz` | Gene information | |
| `└── barcodes.tsv.gz` | Cell barcodes | |
| `_log/` | Log files | Small |

### Files Removed

The function removes space-consuming, non-essential files:

| File/Directory | Description | Typical Size |
|----------------|-------------|--------------|
| `outs/raw_feature_bc_matrix.h5` | Raw (unfiltered) matrix | Very Large |
| `outs/raw_feature_bc_matrix/` | Raw matrix (Matrix Market) | Very Large |
| `outs/analysis/` | Secondary analysis results | Large |
| `├── pca/` | PCA results | |
| `├── tsne/` | t-SNE results | |
| `├── umap/` | UMAP results | |
| `├── clustering/` | Clustering results | |
| `└── diffexp/` | Differential expression | |
| `outs/cloupe.cloupe` | Loupe Cell Browser file | Medium |
| Temporary BAM files | Alignment files | Very Large |
| STAR index cache | Alignment index | Large |

### Storage Savings

Typical storage reduction by dataset size:

| Dataset Size | Original Size | After Reduction | Savings |
|--------------|---------------|-----------------|---------|
| 5K cells | 15 GB | 3 GB | 80% |
| 20K cells | 45 GB | 8 GB | 82% |
| 50K cells | 120 GB | 20 GB | 83% |
| 100K cells | 250 GB | 40 GB | 84% |

## Directory Structure

### Before Reduction

```
resources/SAMPLE_ID/counted/
├── outs/
│   ├── filtered_feature_bc_matrix.h5      # KEEP
│   ├── filtered_feature_bc_matrix/        # KEEP
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   ├── raw_feature_bc_matrix.h5          # REMOVE
│   ├── raw_feature_bc_matrix/            # REMOVE
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   ├── molecule_info.h5                  # KEEP
│   ├── web_summary.html                  # KEEP
│   ├── cloupe.cloupe                     # REMOVE
│   └── analysis/                         # REMOVE
│       ├── pca/
│       ├── tsne/
│       ├── umap/
│       ├── clustering/
│       └── diffexp/
├── _log/                                 # KEEP
└── tmp/                                  # REMOVE
```

### After Reduction

```
resources/SAMPLE_ID/counted/
├── outs/
│   ├── filtered_feature_bc_matrix.h5
│   ├── filtered_feature_bc_matrix/
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   ├── molecule_info.h5
│   └── web_summary.html
└── _log/
```

## Safety Considerations

### Data Preservation

The function preserves all data necessary for:
- Downstream analysis workflows
- Quality control and metrics
- Matrix manipulation and filtering
- Integration with other tools

### Files Safely Removed

Removed files can be regenerated if needed:
- **Raw matrices**: Can be recreated from FASTQ files
- **Secondary analysis**: Can be recomputed from filtered matrices
- **Loupe files**: Can be regenerated from count matrices

### Recovery Options

If removed files are needed later:

```bash
# Re-run Cell Ranger count to regenerate all files
cellranger count \
    --id=SAMPLE_ID \
    --transcriptome=/path/to/transcriptome \
    --fastqs=/path/to/fastqs \
    --sample=SAMPLE_ID
```

## Performance Impact

### Processing Speed

The reduction process is fast:

| Number of Samples | Processing Time |
|-------------------|-----------------|
| 1-5 samples | 10-30 seconds |
| 5-20 samples | 30-60 seconds |
| 20+ samples | 1-2 minutes |

### I/O Considerations

- Uses efficient file system operations
- Processes samples sequentially to avoid I/O contention
- Removes files in order (files before directories)
- Handles symbolic links and empty directories

## Implementation Details

### File Walking Algorithm

```python
# Process files bottom-up to handle directories properly
for foldername, subfolders, filenames in os.walk(
    target_path, topdown=False, followlinks=True
):
    # Remove files not in keep list
    for filename in filenames:
        rel_path = os.path.relpath(
            os.path.join(foldername, filename), target_path
        )
        if rel_path not in FILES_TO_KEEP:
            os.remove(os.path.join(foldername, filename))
    
    # Remove empty directories and symbolic links
    for subfolder in subfolders:
        full_path = os.path.join(foldername, subfolder)
        if os.path.islink(full_path):
            os.remove(full_path)
        elif os.path.isdir(full_path) and not os.listdir(full_path):
            os.rmdir(full_path)
```

### Error Handling

The function handles various edge cases:
- Permission errors (skips protected files)
- Symbolic links (removes safely)
- Empty directories (removes after contents)
- Missing files (continues processing)

## Methods

### `call(project: Project) -> Project`

Main execution method that reduces storage for all counted samples.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Iterates through all samples in the project
2. Checks if counting output exists
3. Walks directory structure bottom-up
4. Removes files not in preservation list
5. Cleans up empty directories and links

## Best Practices

### When to Use

**Recommended scenarios**:
- After completing Cell Ranger count
- Before long-term storage or archival
- When disk space is limited
- For production pipelines

**Avoid when**:
- You need Loupe Cell Browser files
- Secondary analysis results are required
- Disk space is not a concern
- You may need to regenerate raw matrices

### Pipeline Timing

```python
# Optimal pipeline timing
project.call(Count(nthread=8))           # Count reads
# ... extract any needed secondary analysis
project.call(Reduce())                   # Reduce storage
project.call(Preprocess())               # Continue pipeline
```

### Backup Strategy

Consider backing up before reduction:

```bash
# Create backup of essential samples
tar -czf backup_sample_X.tar.gz resources/SAMPLE_X/counted/

# Run reduction
celline run reduce

# Remove backup after confirming success
rm backup_sample_X.tar.gz
```

## Integration with Cluster Storage

### Shared File Systems

For shared cluster storage:

```python
# Check available space before reduction
import shutil

total, used, free = shutil.disk_usage("/shared/project/")
print(f"Free space: {free // (2**30)} GB")

# Run reduction to free space
project.call(Reduce())

# Check space after reduction
total, used, free = shutil.disk_usage("/shared/project/")
print(f"Free space after reduction: {free // (2**30)} GB")
```

### Automated Cleanup

For automated pipelines:

```python
# Automated storage management
def smart_reduce(project, threshold_gb=100):
    """Reduce storage if free space below threshold."""
    import shutil
    
    total, used, free = shutil.disk_usage(Config.PROJ_ROOT)
    free_gb = free // (2**30)
    
    if free_gb < threshold_gb:
        print(f"Low disk space ({free_gb} GB), running reduction...")
        project.call(Reduce())
    else:
        print(f"Sufficient disk space ({free_gb} GB), skipping reduction")

# Use in pipeline
smart_reduce(project, threshold_gb=100)
```

## Monitoring and Reporting

### Storage Reports

Generate storage usage reports:

```python
import os

def storage_report(project_root):
    """Generate storage usage report."""
    total_size = 0
    sample_sizes = {}
    
    for sample_dir in os.listdir(f"{project_root}/resources"):
        sample_path = f"{project_root}/resources/{sample_dir}"
        if os.path.isdir(sample_path):
            size = sum(
                os.path.getsize(os.path.join(dirpath, filename))
                for dirpath, dirnames, filenames in os.walk(sample_path)
                for filename in filenames
            )
            sample_sizes[sample_dir] = size // (2**30)  # GB
            total_size += size
    
    print(f"Total storage: {total_size // (2**30)} GB")
    for sample, size in sorted(sample_sizes.items()):
        print(f"  {sample}: {size} GB")
```

## Related Functions

- [Count](count) - Generate output files before reduction
- [Info](info) - Check storage usage and file status
- [Download](download) - Download source data (kept separate from count output)

## Troubleshooting

### Common Issues

1. **Permission Errors**: Ensure write permissions on count directories
2. **Files in Use**: Close any applications accessing count files
3. **Symbolic Links**: Function handles links safely
4. **Network Storage**: May be slower on network file systems

### Recovery from Errors

If reduction is interrupted:

```python
# Check which samples were processed
processed_samples = []
for sample in project.samples:
    count_dir = f"resources/{sample}/counted"
    if os.path.exists(f"{count_dir}/outs/raw_feature_bc_matrix.h5"):
        print(f"{sample}: Not reduced")
    else:
        print(f"{sample}: Already reduced")
        processed_samples.append(sample)
```

### Manual File Removal

For manual cleanup:

```bash
# Remove specific file types manually
find resources/*/counted -name "raw_feature_bc_matrix.h5" -delete
find resources/*/counted -name "cloupe.cloupe" -delete
find resources/*/counted -name "analysis" -type d -exec rm -rf {} +
```