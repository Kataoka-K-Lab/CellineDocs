# SyncDB Function

Synchronize database information for samples in your project.

## Overview

The `SyncDB` function synchronizes metadata from external databases (GEO, SRA, etc.) for all samples in your project. It ensures that your local database cache is up-to-date with the latest information from public repositories.

## Class Information

- **Module**: `celline.functions.sync_DB`
- **Class**: `SyncDB`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `force_update_target` | `Optional[List[str]]` | No | List of specific sample IDs to force update |

## Usage Examples

### Python API

#### Basic Synchronization

```python
from celline import Project
from celline.functions.sync_DB import SyncDB

# Create project
project = Project("./my-project")

# Create SyncDB function instance
sync_function = SyncDB()

# Execute function
result = project.call(sync_function)
```

#### Force Update Specific Samples

```python
from celline import Project
from celline.functions.sync_DB import SyncDB

# Create project
project = Project("./my-project")

# Force update specific samples
force_update_samples = ["GSM123456", "GSM789012"]
sync_function = SyncDB(force_update_target=force_update_samples)

# Execute function
result = project.call(sync_function)
```

#### Full Project Synchronization

```python
from celline import Project
from celline.functions.sync_DB import SyncDB

# Create project
project = Project("./my-project")

# Synchronize all samples in project
sync_function = SyncDB()
result = project.call(sync_function)

print(f"Synchronized {len(result)} samples")
```

### CLI Usage

#### Basic Usage

```bash
# Synchronize all samples in project
celline run syncdb

# Synchronize with verbose output
celline run syncdb --verbose
```

## Implementation Details

### Synchronization Process

1. **Sample Discovery**: Reads `samples.toml` to identify all project samples
2. **Handler Resolution**: Determines appropriate database handler for each sample
3. **Metadata Fetching**: Retrieves latest metadata from external databases
4. **Local Update**: Updates local database cache with new information
5. **Progress Tracking**: Provides real-time progress feedback

### Database Sources

The function synchronizes with multiple external databases:

- **GEO (Gene Expression Omnibus)**: Sample and series metadata
- **SRA (Sequence Read Archive)**: Run and experiment information
- **CNCB (China National Center for Bioinformation)**: Additional metadata sources
- **Custom Sources**: Project-specific database endpoints

### File Dependencies

The function requires the following project files:

```
project_root/
├── samples.toml          # Sample registry (required)
├── setting.toml          # Project configuration
└── resources/
    └── samples/          # Sample-specific directories
```

### Sample Registry Format

The `samples.toml` file format:

```toml
[samples]
"GSM123456" = "Control sample"
"GSM789012" = "Treatment sample"
"SRR567890" = "RNA-seq run"
```

## Synchronization Behavior

### Normal Mode

In normal mode, the function:
- Checks if metadata already exists locally
- Skips samples with recent metadata (cached)
- Only fetches missing or outdated information
- Optimizes network requests

### Force Update Mode

When `force_update_target` is specified:
- Ignores local cache for specified samples
- Forces fresh metadata retrieval
- Useful for correcting outdated information
- Slower but ensures accuracy

## Error Handling

### Common Issues

1. **Network Connectivity**: Handles timeout and connection errors
2. **Invalid Sample IDs**: Gracefully handles unrecognized accession formats
3. **Database Unavailability**: Retries failed requests with exponential backoff
4. **File Access**: Handles permission and file system errors

### Recovery Strategies

The function implements several recovery mechanisms:

```python
# Automatic retry for network errors
try:
    handler.add(sample, force_search)
except NetworkError:
    # Retry with exponential backoff
    retry_with_backoff(handler.add, sample, force_search)
```

## Methods

### `call(project: Project) -> SyncDB`

Main execution method that performs database synchronization.

**Parameters**:
- `project`: The Celline project instance

**Returns**: SyncDB instance (for method chaining)

**Raises**:
- `FileNotFoundError`: If `samples.toml` file is missing
- `NotImplementedError`: If handler cannot be resolved for a sample

## Performance Considerations

### Batch Processing

For large projects with many samples:

```python
# Process samples in batches to avoid memory issues
batch_size = 10
sample_batches = [samples[i:i+batch_size] for i in range(0, len(samples), batch_size)]

for batch in sample_batches:
    sync_function = SyncDB(force_update_target=batch)
    project.call(sync_function)
```

### Network Optimization

- Implements request pooling for efficiency
- Uses compression for large metadata transfers
- Caches frequently accessed information
- Respects database rate limits

### Memory Usage

- Processes samples sequentially to minimize memory footprint
- Clears metadata cache periodically
- Optimizes data structures for large datasets

## Database Schema Updates

The function automatically handles database schema changes:

```python
# Metadata structure is versioned and backward-compatible
sample_metadata = {
    "version": "2.1",
    "accession": "GSM123456",
    "title": "Sample title",
    "organism": "Homo sapiens",
    "platform": "Illumina HiSeq 2500",
    "last_updated": "2024-01-15T10:30:00Z"
}
```

## Integration with Other Functions

### Workflow Integration

Typical usage in a complete workflow:

```python
from celline import Project
from celline.functions.add import Add
from celline.functions.sync_DB import SyncDB
from celline.functions.download import Download

# Create project and add samples
project = Project("./my-project")
project.call(Add([Add.SampleInfo(id="GSE123456")]))

# Synchronize metadata
project.call(SyncDB())

# Download data (benefits from updated metadata)
project.call(Download())
```

### Dependency Chain

The SyncDB function is typically used:
1. After adding new samples with [Add](add)
2. Before downloading data with [Download](download)
3. Periodically to refresh metadata
4. When troubleshooting data issues

## Monitoring and Logging

### Progress Tracking

The function provides detailed progress information:

```
Fetching... ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 25/25 samples
✓ GSM123456: Updated metadata
✓ GSM789012: Cached (no update needed)
✗ GSM999999: Failed to resolve
```

### Logging Output

Enable detailed logging for debugging:

```python
import logging
logging.basicConfig(level=logging.INFO)

# Detailed sync information will be logged
project.call(SyncDB())
```

## CLI Arguments

While the current implementation doesn't expose CLI arguments, future versions may include:

| Argument | Type | Description |
|----------|------|-------------|
| `--force` | `flag` | Force update all samples |
| `--samples` | `str+` | Specific samples to sync |
| `--timeout` | `int` | Network timeout in seconds |
| `--retry` | `int` | Number of retry attempts |

## Related Functions

- [Add](add) - Add samples before synchronization
- [Download](download) - Download data after synchronization
- [Info](info) - View synchronized metadata
- [Initialize](init) - Initialize project structure

## Troubleshooting

### Common Issues

1. **Missing samples.toml**: Ensure project is properly initialized
2. **Network timeout**: Check internet connection and database availability
3. **Permission errors**: Verify read/write permissions in project directory
4. **Memory issues**: Process large projects in batches

### Manual Database Update

For debugging, you can manually update specific samples:

```python
from celline.DB.dev.handler import HandleResolver

# Manually resolve and update a sample
handler = HandleResolver.resolve("GSM123456")
if handler:
    handler.add("GSM123456", force_search=True)
```

### Cache Management

Clear local database cache if needed:

```python
import shutil
from celline.config import Config

# Clear all cached metadata (use with caution)
cache_dir = f"{Config.PROJ_ROOT}/.cache"
if os.path.exists(cache_dir):
    shutil.rmtree(cache_dir)
```

## Best Practices

1. **Regular Synchronization**: Run sync periodically for active projects
2. **Force Updates**: Use force mode sparingly to avoid unnecessary network load
3. **Error Handling**: Always check return values and handle exceptions
4. **Batch Processing**: For large projects, process samples in batches
5. **Network Considerations**: Run during off-peak hours for large synchronizations