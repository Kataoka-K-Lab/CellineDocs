# Add Function

Add accession IDs to your project and automatically fetch their metadata from public databases.

## Overview

The `Add` function allows you to add sample accession IDs (such as GSE, GSM, SRR) to your Celline project. It automatically resolves the appropriate database handler and retrieves metadata for the specified samples.

## Class Information

- **Module**: `celline.functions.add`
- **Class**: `Add`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sample_id` | `Union[list[Add.SampleInfo], pl.DataFrame]` | Yes | Accession IDs to add to the project |

### SampleInfo Structure

The `Add.SampleInfo` is a NamedTuple with the following fields:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | `str` | Yes | Sample accession ID (e.g., GSE123456, GSM789012) |
| `title` | `Optional[str]` | No | Custom title for the sample (defaults to empty string) |

## Usage Examples

### Python API

#### Adding Single Sample

```python
from celline import Project
from celline.functions.add import Add

# Create project
project = Project("./my-project")

# Create sample info
sample_info = Add.SampleInfo(id="GSE123456", title="My RNA-seq study")

# Create Add function instance
add_function = Add([sample_info])

# Execute function
result = project.call(add_function)
```

#### Adding Multiple Samples

```python
from celline import Project
from celline.functions.add import Add

# Create project
project = Project("./my-project")

# Create multiple sample infos
sample_infos = [
    Add.SampleInfo(id="GSM789012", title="Control sample"),
    Add.SampleInfo(id="GSM789013", title="Treatment sample"),
    Add.SampleInfo(id="GSM789014", title="Replicate sample")
]

# Create Add function instance
add_function = Add(sample_infos)

# Execute function
result = project.call(add_function)
```

#### Using Polars DataFrame

```python
import polars as pl
from celline import Project
from celline.functions.add import Add

# Create project
project = Project("./my-project")

# Create DataFrame with sample information
df = pl.DataFrame({
    "id": ["GSE123456", "GSM789012", "GSM789013"],
    "title": ["Study 1", "Control", "Treatment"]
})

# Create Add function instance
add_function = Add(df)

# Execute function
result = project.call(add_function)
```

### CLI Usage

#### Basic Usage

```bash
# Add single accession ID
celline run add GSE123456

# Add multiple accession IDs
celline run add GSM789012 GSM789013 GSM789014

# Add with custom title
celline run add GSM789012 --title "My sample"
```

#### File Input

```bash
# Read sample IDs from text file (one per line)
celline run add --from-file samples.txt

# Read from CSV file (requires 'id' and 'title' columns)
celline run add --from-file samples.csv

# Read from TSV file
celline run add --from-file samples.tsv
```

### CLI Arguments

| Argument | Type | Required | Description |
|----------|------|----------|-------------|
| `sample_ids` | `str+` | Yes* | One or more sample IDs to add |
| `--title`, `-t` | `str` | No | Optional title for the samples |
| `--from-file`, `-f` | `str` | No | Read sample IDs from file |

*Required if `--from-file` is not used

### File Formats

#### Text File Format
```
GSE123456
GSM789012
GSM789013
```

#### CSV/TSV Format
```csv
id,title
GSE123456,RNA-seq study
GSM789012,Control sample
GSM789013,Treatment sample
```

## Implementation Details

### Database Resolution

The function uses `HandleResolver` to automatically determine the appropriate database handler based on the accession ID format:

- **GSE IDs**: Gene Expression Omnibus (GEO) series
- **GSM IDs**: GEO samples  
- **SRR IDs**: Sequence Read Archive (SRA) runs
- **Other formats**: Automatically detected

### Project Integration

Added samples are stored in the project's `samples.toml` file with the following structure:

```toml
[samples]
"GSM789012" = "Control sample"
"GSM789013" = "Treatment sample"
```

### Progress Tracking

The function provides real-time progress feedback:
- Rich progress bar for multiple samples
- Color-coded status messages (success, warning, error)
- Detailed logging for debugging

### Error Handling

The function handles various error scenarios:
- Invalid accession ID formats
- Network connectivity issues
- Database access problems
- File reading errors (for CLI file input)

## Methods

### `get_samples() -> dict[str, str]`

Retrieves existing sample information from the project's `samples.toml` file.

**Returns**: Dictionary mapping sample IDs to titles

### `call(project: Project) -> Project`

Main execution method that adds the specified samples to the project.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

### `add_cli_args(parser: ArgumentParser) -> None`

Adds CLI-specific arguments to the argument parser.

### `cli(project: Project, args: Namespace) -> Project`

CLI entry point that processes command-line arguments and executes the function.

## Notes

- Parallel calculations are not supported for this function
- Duplicate sample IDs are automatically skipped
- The function requires internet connectivity to fetch metadata
- Large datasets may take considerable time to process

## Related Functions

- [Download](download) - Download the added samples
- [SyncDB](syncdb) - Synchronize database information
- [Info](info) - View project information

## Troubleshooting

### Common Issues

1. **Network Connection**: Ensure stable internet connection for metadata retrieval
2. **Invalid IDs**: Verify accession ID format (GSE, GSM, SRR, etc.)
3. **File Permissions**: Check read permissions for input files in CLI mode
4. **Memory Usage**: Large sample lists may require adequate memory

### Debug Mode

Enable detailed logging for troubleshooting:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```