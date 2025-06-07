# Download Function

Download raw sequencing data files for samples in your project.

## Overview

The `Download` function automatically downloads raw sequencing data files (FASTQ, BAM, etc.) for all samples that have been added to your project but not yet downloaded. It supports multi-threading and provides comprehensive progress tracking.

## Class Information

- **Module**: `celline.functions.download`
- **Class**: `Download`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `then` | `Optional[Callable[[str], None]]` | No | Callback function executed on successful download |
| `catch` | `Optional[Callable[[subprocess.CalledProcessError], None]]` | No | Callback function executed on download error |

### JobContainer Structure

The `Download.JobContainer` is a NamedTuple containing job information:

| Field | Type | Description |
|-------|------|-------------|
| `filetype` | `str` | File type strategy for download |
| `nthread` | `str` | Number of threads to use |
| `cluster_server` | `str` | Cluster server name (if applicable) |
| `jobname` | `str` | Job identifier |
| `logpath` | `str` | Path to log file |
| `sample_id` | `str` | Sample accession ID |
| `download_target` | `str` | Target directory for downloads |
| `download_source` | `str` | Source URL for download |
| `run_ids_str` | `str` | Comma-separated run IDs |

## Usage Examples

### Python API

#### Basic Usage

```python
from celline import Project
from celline.functions.download import Download

# Create project
project = Project("./my-project")

# Create Download function instance
download_function = Download()

# Execute function
result = project.call(download_function)
```

#### With Threading

```python
from celline import Project
from celline.functions.download import Download

# Create project
project = Project("./my-project")

# Create Download function instance with threading
download_function = Download()
download_function.nthread = 4  # Use 4 threads

# Execute function
result = project.call(download_function)
```

#### With Callbacks

```python
from celline import Project
from celline.functions.download import Download
import subprocess

def on_success(sample_id: str):
    print(f"Successfully downloaded: {sample_id}")

def on_error(error: subprocess.CalledProcessError):
    print(f"Download failed: {error}")

# Create project
project = Project("./my-project")

# Create Download function instance with callbacks
download_function = Download(then=on_success, catch=on_error)

# Execute function
result = project.call(download_function)
```

### CLI Usage

#### Basic Download

```bash
# Download data for all samples
celline run download

# Download with multiple threads for faster processing
celline run download --nthread 4

# Force re-download existing files
celline run download --force

# Combine options
celline run download --nthread 2 --force
```

#### CLI Arguments

| Argument | Short | Type | Default | Description |
|----------|-------|------|---------|-------------|
| `--nthread` | `-n` | `int` | `1` | Number of threads for parallel download |
| `--force` | `-f` | `flag` | `False` | Force re-download existing files |

**Understanding CLI options:**
- **Threading**: The `--nthread` parameter enables parallel downloads, significantly reducing download time for multiple samples
- **Force mode**: Use `--force` when you need to re-download corrupted files or update existing data
- **Progress monitoring**: Check `resources/[SAMPLE_ID]/log/download_*.log` files for detailed progress and error information

## Implementation Details

### Download Process

1. **Sample Resolution**: Identifies all samples in the project
2. **Handler Resolution**: Determines appropriate database handler for each sample
3. **Schema Retrieval**: Fetches sample and run metadata
4. **Path Preparation**: Creates necessary directory structure
5. **Template Generation**: Creates download scripts from templates
6. **Execution**: Runs download jobs with thread management

### File Organization

Downloaded files are organized in the following structure:

```
project_root/
├── resources/
│   └── SAMPLE_ID/
│       ├── raw/
│       │   └── fastqs/          # Downloaded FASTQ files
│       ├── src/
│       │   └── download.sh      # Generated download script
│       └── log/
│           └── download_*.log   # Download logs
```

### Download Scripts

The function generates shell scripts for each sample using templates:

```bash
#!/bin/bash
# Generated download script for SAMPLE_ID

# Download configuration
NTHREAD=4
SAMPLE_ID="GSM123456"
DOWNLOAD_TARGET="/path/to/raw/data"
DOWNLOAD_SOURCE="ftp://ftp.ncbi.nlm.nih.gov/..."

# Execute download
fastq-dump --split-files --outdir $DOWNLOAD_TARGET $RUN_IDS
```

### Thread Management

The function uses `ThreadObservable` for concurrent downloads:
- Each sample gets its own download job
- Jobs are executed in parallel up to the specified thread limit
- Progress is tracked and reported in real-time

### Status Checking

The function checks download status before processing:
- Skips samples already downloaded (`path.is_downloaded`)
- Skips samples already processed (`path.is_counted`)
- Optionally force re-download with CLI flag

## File Types Supported

The function supports various sequencing file formats:

- **FASTQ**: Single-cell RNA-seq data
- **BAM/SAM**: Aligned sequencing data
- **SRA**: Sequence Read Archive format
- **CRAM**: Compressed reference-aligned format

## Error Handling

### Common Issues

1. **Network Connectivity**: Handles timeout and connection errors
2. **Storage Space**: Checks available disk space before download
3. **File Corruption**: Validates downloaded file integrity
4. **Permission Issues**: Handles file system permission errors

### Recovery Mechanisms

- Automatic retry for failed downloads
- Partial download resumption
- Corrupted file detection and re-download
- Detailed error logging for debugging

## Methods

### `call(project: Project) -> Project`

Main execution method that orchestrates the download process.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

### `add_cli_args(parser: ArgumentParser) -> None`

Adds CLI-specific arguments to the argument parser.

### `cli(project: Project, args: Namespace) -> Project`

CLI entry point that processes command-line arguments and executes the function.

## Performance Considerations

### Threading

- **Optimal Thread Count**: 2-4 threads for most systems
- **Network Bandwidth**: Consider your connection speed
- **Storage I/O**: Balance threads with storage performance

### Memory Usage

- Large datasets may require significant memory
- Monitor system resources during download
- Consider batch processing for very large projects

### Storage Requirements

Example storage requirements for common datasets:

| Dataset Type | File Size Range | Storage Needed |
|--------------|----------------|----------------|
| Single-cell RNA-seq | 1-10 GB | 10-100 GB |
| Bulk RNA-seq | 5-50 GB | 50-500 GB |
| ATAC-seq | 2-20 GB | 20-200 GB |

## Integration

### Pipeline Integration

The Download function integrates with other Celline functions:

```python
from celline import Project
from celline.functions.add import Add
from celline.functions.download import Download
from celline.functions.count import Count

# Create project and add samples
project = Project("./my-project")
project.call(Add([Add.SampleInfo(id="GSE123456")]))

# Download data
project.call(Download())

# Process data
project.call(Count())
```

### Cluster Computing

For cluster environments:

```python
# Configure cluster settings before download
from celline.server import ServerSystem
ServerSystem.cluster_server_name = "slurm"

# Downloads will be submitted as cluster jobs
project.call(Download())
```

## Related Functions

- [Add](add) - Add samples to project before downloading
- [Count](count) - Process downloaded data
- [Info](info) - Check download status
- [Job](job) - Monitor download jobs

## Troubleshooting

### Common Issues

1. **Slow Downloads**: Check network connection and reduce thread count
2. **Disk Space**: Ensure sufficient storage space
3. **Permission Errors**: Check file system permissions
4. **Network Timeouts**: Increase timeout settings or retry

### Debug Mode

Enable detailed logging:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

### Manual Download

For problematic samples, you can manually download:

```bash
# Navigate to sample directory
cd resources/SAMPLE_ID/src/

# Execute download script manually
bash download.sh
```