# Info Function

Display comprehensive project information, status, and system details.

## Overview

The `Info` function provides a detailed overview of your Celline project, including project paths, configuration settings, sample status, system information, and analysis progress. It serves as a quick diagnostic tool for project management and troubleshooting.

## Class Information

- **Module**: `celline.functions.info`
- **Class**: `Info`
- **Base Class**: `CellineFunction`

## Parameters

The `Info` function takes no constructor parameters.

## Usage Examples

### Python API

#### Basic Information Display

```python
from celline import Project
from celline.functions.info import Info

# Create project
project = Project("./my-project")

# Display project information
info_function = Info()

# Execute function
result = project.call(info_function)
```

#### Programmatic Information Access

```python
from celline import Project
from celline.functions.info import Info

# Create project
project = Project("./my-project")

# Get project information
info_function = Info()
result = info_function.call(project)

# Access project details programmatically
print(f"Project path: {project.PROJ_PATH}")
print(f"Execution path: {project.EXEC_PATH}")
```

### CLI Usage

#### Basic Usage

```bash
# Display project information
celline run info

# Use info for project diagnostics
celline info  # Alternative command
```

#### Pipeline Integration

```bash
# Check project status before analysis
celline run info
celline run count --nthread 8
celline run info  # Check again after processing
```

## Implementation Details

### Displayed Information

The function provides comprehensive project information including:

#### Project Paths
- **Project Root**: Main project directory location
- **Execution Path**: Current execution context
- **Data Directory**: Location of processed data
- **Resources Directory**: Raw data and intermediate files
- **Logs Directory**: Analysis logs and outputs

#### Project Configuration
- **Project Name**: From settings configuration
- **Version**: Project version information
- **Settings File**: Configuration file location
- **Samples File**: Sample registry location

#### Sample Status
- **Total Samples**: Number of samples in project
- **Downloaded**: Samples with raw data
- **Counted**: Samples processed with Cell Ranger
- **Preprocessed**: Samples with quality control
- **Predicted**: Samples with cell type predictions

#### System Information
- **Celline Version**: Current software version
- **Python Version**: Runtime environment
- **R Version**: R environment status
- **Dependencies**: Required software availability

### Extended Information Display

For comprehensive diagnostics:

```python
from celline import Project
from celline.functions.info import Info
from celline.config import Config, Setting
from celline.sample import SampleResolver
import os

class ExtendedInfo(Info):
    def call(self, project):
        # Basic project info
        super().call(project)
        
        # Extended system information
        console.print("\n[cyan]Extended Information:[/cyan]")
        
        # Storage usage
        self._show_storage_info()
        
        # Sample details
        self._show_sample_details()
        
        # Dependency status
        self._show_dependency_status()
        
        # Recent activity
        self._show_recent_activity()
        
        return project
    
    def _show_storage_info(self):
        """Display storage usage information."""
        import shutil
        
        total, used, free = shutil.disk_usage(Config.PROJ_ROOT)
        console.print(f"Storage usage:")
        console.print(f"  Total: {total // (2**30)} GB")
        console.print(f"  Used: {used // (2**30)} GB") 
        console.print(f"  Free: {free // (2**30)} GB")
    
    def _show_sample_details(self):
        """Display detailed sample information."""
        console.print(f"\nSample Status:")
        for sample_id, sample_info in SampleResolver.samples.items():
            status = []
            if sample_info.path.is_downloaded:
                status.append("downloaded")
            if sample_info.path.is_counted:
                status.append("counted")
            if sample_info.path.is_preprocessed:
                status.append("preprocessed")
            if sample_info.path.is_predicted_celltype:
                status.append("predicted")
            
            status_str = ", ".join(status) if status else "new"
            console.print(f"  {sample_id}: {status_str}")
```

## Information Categories

### Project Structure

```
Project Information Display:
├── Project Paths
│   ├── Root Directory: /path/to/project
│   ├── Data Directory: /path/to/project/data
│   ├── Resources Directory: /path/to/project/resources
│   └── Results Directory: /path/to/project/results
├── Configuration
│   ├── Project Name: My Single Cell Project
│   ├── Version: 1.0.0
│   ├── Settings File: setting.toml
│   └── Samples File: samples.toml
├── Sample Status
│   ├── Total Samples: 15
│   ├── Downloaded: 12
│   ├── Counted: 10
│   ├── Preprocessed: 8
│   └── Predicted: 6
└── System Information
    ├── Celline Version: 1.0.0
    ├── Python Version: 3.9.7
    ├── R Version: 4.1.0
    └── Cell Ranger: Available
```

### Sample Status Indicators

| Status | Description | Files Present |
|--------|-------------|---------------|
| **New** | Sample added but not processed | `samples.toml` entry only |
| **Downloaded** | Raw data downloaded | `resources/SAMPLE/raw/` files |
| **Counted** | Cell Ranger processing completed | `resources/SAMPLE/counted/outs/` |
| **Preprocessed** | Quality control applied | `data/SAMPLE/cell_info.tsv` |
| **Predicted** | Cell type prediction completed | `data/SAMPLE/predicted_celltype.tsv` |
| **Integrated** | Included in integration analysis | `integration/` outputs |

### System Dependencies

The function checks availability of required tools:

| Tool | Purpose | Status Check |
|------|---------|-------------|
| **Python** | Runtime environment | Version and packages |
| **R** | Statistical computing | Version and libraries |
| **Cell Ranger** | Count processing | Installation path |
| **FastQ-dump** | SRA data download | Command availability |
| **Git** | Version control | Repository status |

## Methods

### `register() -> str`

Returns the function identifier for registration.

### `call(project: Project) -> Project`

Main execution method that displays project information.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Output**: Rich formatted project information display

### `cli(project: Project, args: Namespace) -> Project`

CLI entry point that displays project information.

### `get_description() -> str`

Returns function description for help documentation.

### `get_usage_examples() -> list[str]`

Returns usage examples for CLI help.

## Advanced Usage

### Custom Information Display

```python
from celline import Project
from celline.functions.info import Info
from rich.table import Table
from rich.console import Console

class CustomInfo(Info):
    def call(self, project):
        console = Console()
        
        # Create custom information table
        table = Table(title="Project Dashboard")
        table.add_column("Category", style="cyan")
        table.add_column("Status", style="green")
        table.add_column("Details", style="white")
        
        # Add project information rows
        table.add_row("Project", "Active", project.PROJ_PATH)
        table.add_row("Samples", "15 total", "12 processed")
        table.add_row("Storage", "45 GB", "120 GB available")
        
        console.print(table)
        return project
```

### Automated Health Checks

```python
from celline import Project
from celline.functions.info import Info

class HealthCheck(Info):
    def call(self, project):
        console = Console()
        issues = []
        
        # Check for common issues
        if not os.path.exists(f"{project.PROJ_PATH}/samples.toml"):
            issues.append("Missing samples.toml file")
        
        if not os.path.exists(f"{project.PROJ_PATH}/setting.toml"):
            issues.append("Missing setting.toml file")
        
        # Check disk space
        import shutil
        _, _, free = shutil.disk_usage(project.PROJ_PATH)
        if free < 10 * (2**30):  # Less than 10 GB
            issues.append("Low disk space")
        
        # Display health status
        if issues:
            console.print("[red]Health Check Issues:[/red]")
            for issue in issues:
                console.print(f"  ⚠ {issue}")
        else:
            console.print("[green]✓ All health checks passed[/green]")
        
        return project
```

### Integration with Monitoring

```python
import json
import datetime
from celline.functions.info import Info

class MonitoringInfo(Info):
    def call(self, project):
        # Collect project metrics
        metrics = {
            "timestamp": datetime.datetime.now().isoformat(),
            "project_path": project.PROJ_PATH,
            "sample_count": len(SampleResolver.samples),
            "processed_samples": sum(
                1 for s in SampleResolver.samples.values() 
                if s.path.is_counted
            ),
            "storage_usage": self._get_storage_usage(),
            "system_status": self._check_system_status()
        }
        
        # Save metrics for monitoring
        with open(f"{project.PROJ_PATH}/metrics.json", "w") as f:
            json.dump(metrics, f, indent=2)
        
        # Display summary
        console.print(f"[cyan]Monitoring metrics saved to metrics.json[/cyan]")
        
        return project
```

## Output Formats

### Console Display

The default output provides formatted console display:

```
Project Information:
Project path: /Users/researcher/my-project
Exec path: /Users/researcher/my-project

Storage Information:
├── Total space: 500 GB
├── Used space: 120 GB
└── Available: 380 GB

Sample Status:
├── Total samples: 15
├── Downloaded: 12
├── Counted: 10
├── Preprocessed: 8
└── Predicted: 6

System Dependencies:
├── Python 3.9.7: ✓ Available
├── R 4.1.0: ✓ Available
├── Cell Ranger: ✓ Available
└── FastQ-dump: ✓ Available
```

### JSON Export

For programmatic access:

```python
def export_info_json(project):
    info_data = {
        "project": {
            "name": Setting.name,
            "path": project.PROJ_PATH,
            "version": Setting.version
        },
        "samples": {
            "total": len(SampleResolver.samples),
            "status_counts": count_sample_status()
        },
        "system": get_system_info(),
        "timestamp": datetime.datetime.now().isoformat()
    }
    
    with open("project_info.json", "w") as f:
        json.dump(info_data, f, indent=2)
```

## Performance Considerations

### Fast Information Retrieval

The function is optimized for quick information display:

- **Minimal I/O**: Only reads necessary configuration files
- **Cached Status**: Uses cached sample status when available
- **Lazy Loading**: Loads detailed information only when requested

### Large Project Handling

For projects with many samples:

```python
class OptimizedInfo(Info):
    def call(self, project):
        # Use sampling for large projects
        sample_count = len(SampleResolver.samples)
        
        if sample_count > 100:
            # Sample-based status estimation
            sample_subset = list(SampleResolver.samples.values())[:10]
            estimated_status = self._estimate_status(sample_subset)
            console.print(f"Status (estimated from {len(sample_subset)} samples):")
        else:
            # Full status check for smaller projects
            actual_status = self._get_full_status()
            console.print("Status (complete):")
```

## Related Functions

- [Initialize](init) - Set up project structure
- [Job](job) - Monitor running jobs
- [Bash](bash) - Execute diagnostic commands
- All analysis functions - Check their prerequisites and outputs

## Troubleshooting

### Common Information Issues

1. **Missing Configuration**: Project not properly initialized
2. **Permission Errors**: Cannot read project files
3. **Corrupted Status**: Sample status inconsistent
4. **Large Projects**: Slow information retrieval

### Debug Information

Enable detailed diagnostic output:

```python
class DebugInfo(Info):
    def call(self, project):
        console.print("[yellow]Debug Information:[/yellow]")
        
        # File system permissions
        self._check_permissions()
        
        # Configuration validation
        self._validate_config()
        
        # Sample integrity
        self._check_sample_integrity()
        
        return project
```

### Manual Information Gathering

For troubleshooting, gather information manually:

```bash
# Project structure
ls -la
cat setting.toml
cat samples.toml

# Sample status
find resources/ -name "*.h5" | wc -l
find data/ -name "cell_info.tsv" | wc -l

# System dependencies
which cellranger
R --version
python --version
```