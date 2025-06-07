# SetTranscriptome Function

Register reference transcriptome paths for Cell Ranger count operations.

## Overview

The `SetTranscriptome` function registers reference transcriptome paths in the Celline database. These references are required for Cell Ranger count operations and must be properly formatted Cell Ranger-compatible reference packages.

## Class Information

- **Module**: `celline.functions.set_transcriptome`
- **Class**: `SetTranscriptome`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `species` | `str` | Yes | Species name (e.g., "Homo sapiens", "Mus musculus") |
| `built_path` | `str` | Yes | Path to Cell Ranger reference directory |
| `force_update` | `bool` | No | Force update existing entry (default: True) |

## Usage Examples

### Python API

#### Register Human Reference

```python
from celline import Project
from celline.functions.set_transcriptome import SetTranscriptome

# Create project
project = Project("./my-project")

# Register human reference transcriptome
transcriptome_function = SetTranscriptome(
    species="Homo sapiens",
    built_path="/path/to/refdata-gex-GRCh38-2020-A"
)

# Execute function
result = project.call(transcriptome_function)
```

#### Register Mouse Reference

```python
from celline import Project
from celline.functions.set_transcriptome import SetTranscriptome

# Create project
project = Project("./my-project")

# Register mouse reference transcriptome
transcriptome_function = SetTranscriptome(
    species="Mus musculus",
    built_path="/path/to/refdata-gex-mm10-2020-A"
)

# Execute function
result = project.call(transcriptome_function)
```

#### Custom Species Reference

```python
from celline import Project
from celline.functions.set_transcriptome import SetTranscriptome

# Create project
project = Project("./my-project")

# Register custom species reference
transcriptome_function = SetTranscriptome(
    species="Danio rerio",
    built_path="/path/to/custom-zebrafish-reference",
    force_update=True
)

# Execute function
result = project.call(transcriptome_function)
```

#### Multiple References

```python
from celline import Project
from celline.functions.set_transcriptome import SetTranscriptome

# Create project
project = Project("./my-project")

# Register multiple species references
references = [
    ("Homo sapiens", "/path/to/human-ref"),
    ("Mus musculus", "/path/to/mouse-ref"),
    ("Rattus norvegicus", "/path/to/rat-ref")
]

for species, path in references:
    transcriptome_function = SetTranscriptome(
        species=species,
        built_path=path
    )
    project.call(transcriptome_function)
```

### CLI Usage

While the current implementation doesn't expose CLI arguments, you can register transcriptomes programmatically:

```bash
# Example workflow with transcriptome registration
python -c "
from celline import Project
from celline.functions.set_transcriptome import SetTranscriptome

project = Project('.')
project.call(SetTranscriptome(
    species='Homo sapiens',
    built_path='/data/references/GRCh38'
))
"
```

## Reference Requirements

### Cell Ranger Compatibility

Transcriptome references must be Cell Ranger-compatible packages containing:

```
reference_directory/
├── fasta/
│   ├── genome.fa          # Reference genome FASTA
│   └── genome.fa.fai      # FASTA index
├── genes/
│   └── genes.gtf          # Gene annotation GTF
├── pickle/
│   └── ...                # Pickled reference data
├── star/                  # STAR index files
│   ├── Genome
│   ├── SA
│   ├── SAindex
│   └── ...
└── reference.json         # Reference metadata
```

### Pre-built References

10x Genomics provides pre-built references:

| Species | Reference Package | Version |
|---------|------------------|---------|
| Human | `refdata-gex-GRCh38-2020-A` | GRCh38 |
| Mouse | `refdata-gex-mm10-2020-A` | mm10 |
| Human (v3.1.0) | `refdata-gex-GRCh38-2024-A` | GRCh38 |
| Mouse (v3.1.0) | `refdata-gex-mm39-2024-A` | mm39 |

### Download Instructions

```bash
# Download human reference (GRCh38)
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz

# Download mouse reference (mm10)
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
tar -xzf refdata-gex-mm10-2020-A.tar.gz
```

## Implementation Details

### Database Integration

The function updates the Transcriptome database table:

```python
# Internal database operation
Transcriptome().add_path(
    species=self.species,
    built_path=self.built_path,
    force_update=self.force_update
)
```

### Database Schema

The transcriptome database stores:

| Field | Type | Description |
|-------|------|-------------|
| `species` | `str` | Species name (primary key) |
| `built_path` | `str` | Path to reference directory |
| `created_at` | `datetime` | Registration timestamp |
| `updated_at` | `datetime` | Last update timestamp |

### Species Name Conventions

Use standard binomial nomenclature:
- Human: `"Homo sapiens"`
- Mouse: `"Mus musculus"`
- Rat: `"Rattus norvegicus"`
- Zebrafish: `"Danio rerio"`
- Fruit fly: `"Drosophila melanogaster"`

## Force Update Behavior

### Default Behavior (`force_update=True`)

- Overwrites existing entries
- Updates path to new location
- Refreshes metadata timestamps
- Recommended for most use cases

### Preserve Existing (`force_update=False`)

- Skips if species already registered
- Prevents accidental overwrites
- Useful for initialization scripts

## Integration with Count Function

### Automatic Resolution

The Count function automatically resolves transcriptomes:

```python
from celline.functions.count import Count
from celline.functions.set_transcriptome import SetTranscriptome

# Register transcriptome first
project.call(SetTranscriptome(
    species="Homo sapiens",
    built_path="/path/to/human-ref"
))

# Count function will automatically use registered transcriptome
project.call(Count(nthread=8))
```

### Species Matching

Transcriptomes are matched by sample species metadata:
1. Sample species determined from database metadata
2. Corresponding transcriptome path retrieved
3. Used in Cell Ranger count command

## Custom Reference Creation

### Building Custom References

For non-standard organisms, create custom references:

```bash
# Create custom reference with cellranger mkref
cellranger mkref \
    --genome=custom_genome \
    --fasta=genome.fa \
    --genes=annotations.gtf \
    --ref-version=1.0.0
```

### Quality Validation

Validate custom references:

```bash
# Test reference with cellranger count
cellranger count \
    --id=test_run \
    --transcriptome=/path/to/custom-reference \
    --fastqs=/path/to/test-fastqs \
    --sample=test_sample \
    --expect-cells=1000
```

## Methods

### `call(project: Project) -> Project`

Main execution method that registers the transcriptome in the database.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Validates input parameters
2. Calls `Transcriptome().add_path()` with provided parameters
3. Updates database with new/updated entry

## Error Handling

### Common Issues

1. **Invalid Path**: Ensure reference directory exists and is accessible
2. **Permission Errors**: Check read permissions for reference directory
3. **Corrupted Reference**: Validate Cell Ranger reference structure
4. **Disk Space**: Ensure adequate space for reference storage

### Validation Checks

```python
import os

def validate_reference(built_path: str) -> bool:
    """Validate Cell Ranger reference structure."""
    required_files = [
        "reference.json",
        "fasta/genome.fa",
        "genes/genes.gtf"
    ]
    
    for file_path in required_files:
        if not os.path.exists(os.path.join(built_path, file_path)):
            return False
    return True
```

## Storage Considerations

### Reference Sizes

Typical reference package sizes:

| Species | Reference Size | Disk Space |
|---------|---------------|------------|
| Human (GRCh38) | 8-12 GB | 15-20 GB (with index) |
| Mouse (mm10) | 6-8 GB | 12-15 GB (with index) |
| Custom genome | Variable | 2-3x FASTA size |

### Network Storage

For shared environments:

```python
# Use network-accessible paths
transcriptome_function = SetTranscriptome(
    species="Homo sapiens",
    built_path="/shared/references/GRCh38-2020-A"
)
```

## Best Practices

### Reference Management

1. **Version Control**: Include reference version in path names
2. **Centralized Storage**: Use shared locations for multi-user environments
3. **Backup Strategy**: Maintain backups of custom references
4. **Documentation**: Document custom reference build procedures

### Path Conventions

```python
# Recommended path structure
base_path = "/data/cellranger_references"
references = {
    "Homo sapiens": f"{base_path}/GRCh38-2020-A",
    "Mus musculus": f"{base_path}/mm10-2020-A",
    "Custom species": f"{base_path}/custom-species-v1.0"
}
```

## Related Functions

- [Count](count) - Uses registered transcriptomes for counting
- [Info](info) - Display registered transcriptomes
- [Initialize](init) - Set up project with default transcriptomes

## Troubleshooting

### Common Issues

1. **Reference Not Found**: Verify path exists and is accessible
2. **Permission Denied**: Check file system permissions
3. **Invalid Reference**: Validate Cell Ranger reference structure
4. **Database Error**: Check database file permissions

### Verification

Check registered transcriptomes:

```python
from celline.DB.model.transcriptome import Transcriptome

# List all registered transcriptomes
transcriptomes = Transcriptome().list_all()
for species, path in transcriptomes.items():
    print(f"{species}: {path}")
```

### Manual Database Update

For direct database manipulation:

```python
from celline.DB.model.transcriptome import Transcriptome

# Get transcriptome database instance
db = Transcriptome()

# Check if species exists
existing = db.search("Homo sapiens")
if existing:
    print(f"Found: {existing}")
else:
    print("Species not registered")
```