# Initialize Function

Set up a new Celline project with proper directory structure and configuration.

## Overview

The `Initialize` function creates a new Celline project from scratch, setting up the necessary directory structure, configuration files, and validating system dependencies. It provides an interactive setup process that guides users through project configuration.

## Class Information

- **Module**: `celline.functions.initialize`
- **Class**: `Initialize`
- **Base Class**: `CellineFunction`

## Parameters

The `Initialize` function takes no constructor parameters. All configuration is handled through interactive prompts.

## Usage Examples

### Python API

#### Basic Project Initialization

```python
from celline import Project
from celline.functions.initialize import Initialize

# Create project (can be empty directory)
project = Project("./my-new-project")

# Initialize project
init_function = Initialize()

# Execute function (will prompt for configuration)
result = project.call(init_function)
```

#### Automated Initialization

```python
from celline import Project
from celline.functions.initialize import Initialize
from celline.config import Setting

# Pre-configure settings
Setting.name = "My Research Project"
Setting.r_path = "/usr/local/bin/R"
Setting.version = "1.0.0"

# Initialize with existing settings
project = Project("./automated-project")
init_function = Initialize()
result = project.call(init_function)
```

### CLI Usage

#### Interactive Initialization

```bash
# Navigate to project directory
mkdir my-project
cd my-project

# Initialize interactively
celline init

# Alternative: Initialize from anywhere
celline init --project-dir ./my-project
```

#### Quick Setup

```bash
# Initialize with minimal prompts
celline init --quick

# Initialize with specific settings
celline init --name "Brain Atlas Project" --r-path "/opt/R/bin/R"
```

## Implementation Details

### Initialization Process

The function follows a comprehensive setup workflow:

1. **Dependency Validation**: Checks system requirements
2. **R Environment Setup**: Configures R installation
3. **Interactive Configuration**: Prompts for project details
4. **Directory Creation**: Sets up project structure
5. **File Generation**: Creates configuration files
6. **Validation**: Confirms successful setup

### System Dependencies

The function validates required software:

| Dependency | Purpose | Validation |
|------------|---------|------------|
| **Cell Ranger** | Single-cell processing | Command availability |
| **R** | Statistical computing | Version and packages |
| **Python** | Runtime environment | Version compatibility |
| **FastQ-dump** | SRA data access | Installation check |

### Dependency Validation

```python
from celline.utils.dependencies import DependencyValidator

# Comprehensive dependency check
dependencies_ok = DependencyValidator.validate_dependencies(
    show_details=True,
    check_r_packages=False
)

if not dependencies_ok:
    console.print("[red]Missing required dependencies[/red]")
    console.print("[yellow]Please install and retry initialization[/yellow]")
    return project
```

### R Environment Configuration

The function provides intelligent R installation selection:

```python
# Automatic R detection
selected_r_path = DependencyValidator.select_r_installation()

if selected_r_path is None:
    console.print("[red]R installation required[/red]")
    return project

# Configure R path in settings
settings.r_path = selected_r_path
```

## Project Structure

### Created Directory Structure

```
project_root/
├── data/                    # Processed data files
├── results/                 # Analysis results and outputs
├── scripts/                 # Custom analysis scripts
├── resources/              # Raw data and intermediate files (auto-created)
├── logs/                   # Process logs (auto-created)
├── integration/            # Integration outputs (auto-created)
├── batch/                  # Batch correction outputs (auto-created)
├── samples.toml            # Sample configuration
└── setting.toml            # Project settings
```

### Configuration Files

#### samples.toml

Sample registry with examples:

```toml
# Sample configuration for My Research Project
# Add your samples here following this format:
#
# [samples.sample1]
# name = "Sample 1"
# path = "data/sample1"
#
# [samples.sample2] 
# name = "Sample 2"
# path = "data/sample2"

# Example:
# [samples.GSM123456]
# name = "Control sample"
# path = "data/GSM123456"
```

#### setting.toml

Project configuration:

```toml
[project]
name = "My Research Project"
version = "1.0.0"
description = "Single cell analysis project"

[analysis]
# Analysis parameters go here
# threads = 4
# memory_limit = "8G"

[paths]
data_dir = "data"
results_dir = "results"
scripts_dir = "scripts"
```

### Internal Configuration

The function also configures internal Celline settings:

```python
settings = Setting()
settings.name = project_name
settings.r_path = selected_r_path
settings.version = "0.1"
settings.wait_time = 4
settings.flush()  # Save to internal config
```

## Interactive Configuration

### Project Setup Prompts

The function uses inquirer for interactive configuration:

```python
import inquirer

questions = [
    inquirer.Text(
        name="projname",
        message="What is a name of your project?"
    ),
    inquirer.List(
        name="analysis_type",
        message="What type of analysis?",
        choices=["Single-cell RNA-seq", "ATAC-seq", "Multi-omics"]
    ),
    inquirer.Confirm(
        name="use_cluster",
        message="Will you use cluster computing?",
        default=False
    )
]

result = inquirer.prompt(questions, raise_keyboard_interrupt=True)
```

### Advanced Configuration Options

```python
# Extended configuration prompts
advanced_questions = [
    inquirer.Text(
        name="description",
        message="Project description (optional):"
    ),
    inquirer.List(
        name="species",
        message="Primary species:",
        choices=["Homo sapiens", "Mus musculus", "Other"]
    ),
    inquirer.Text(
        name="contact",
        message="Contact email (optional):"
    ),
    inquirer.Path(
        name="data_location",
        message="Default data location:",
        path_type=inquirer.Path.DIRECTORY
    )
]
```

## Methods

### `register() -> str`

Returns the function identifier for registration.

### `call(project: Project) -> Project`

Main execution method that performs project initialization.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Validates system dependencies
2. Configures R environment
3. Prompts for project configuration
4. Creates directory structure
5. Generates configuration files
6. Provides setup guidance

## Error Handling

### Common Issues

1. **Missing Dependencies**: Required software not installed
2. **Permission Errors**: Cannot create directories/files
3. **R Configuration**: R not found or incompatible version
4. **Existing Project**: Directory already contains project files

### Dependency Resolution

```python
def resolve_dependencies():
    missing = []
    
    # Check Cell Ranger
    if not shutil.which("cellranger"):
        missing.append("Cell Ranger")
    
    # Check R
    if not shutil.which("R"):
        missing.append("R")
    
    # Check Python packages
    try:
        import scanpy, pandas, numpy
    except ImportError as e:
        missing.append(f"Python package: {e.name}")
    
    if missing:
        console.print(f"[red]Missing dependencies: {', '.join(missing)}[/red]")
        return False
    
    return True
```

### Recovery Strategies

```python
def handle_init_failure(error_type, project_path):
    if error_type == "permission":
        console.print("[yellow]Try running with sudo or check directory permissions[/yellow]")
    
    elif error_type == "existing_project":
        console.print("[yellow]Directory already contains project files[/yellow]")
        console.print("Use --force to overwrite or choose different directory")
    
    elif error_type == "dependencies":
        console.print("[yellow]Install missing dependencies and retry[/yellow]")
        show_installation_guide()
```

## Validation and Testing

### Post-Initialization Validation

```python
def validate_initialization(project_path):
    required_files = [
        "samples.toml",
        "setting.toml"
    ]
    
    required_dirs = [
        "data",
        "results", 
        "scripts"
    ]
    
    # Check files
    for file in required_files:
        if not os.path.exists(os.path.join(project_path, file)):
            return False, f"Missing file: {file}"
    
    # Check directories
    for directory in required_dirs:
        if not os.path.exists(os.path.join(project_path, directory)):
            return False, f"Missing directory: {directory}"
    
    return True, "Initialization successful"
```

### Integration Testing

```python
def test_project_setup(project_path):
    """Test basic project functionality after initialization."""
    
    # Test configuration loading
    try:
        project = Project(project_path)
        console.print("[green]✓ Project loading successful[/green]")
    except Exception as e:
        console.print(f"[red]✗ Project loading failed: {e}[/red]")
        return False
    
    # Test basic function availability
    try:
        from celline.functions.info import Info
        project.call(Info())
        console.print("[green]✓ Basic functions accessible[/green]")
    except Exception as e:
        console.print(f"[red]✗ Function access failed: {e}[/red]")
        return False
    
    return True
```

## Advanced Features

### Template Projects

```python
def create_from_template(template_name, project_path):
    """Create project from predefined template."""
    
    templates = {
        "single_cell": {
            "description": "Standard single-cell RNA-seq analysis",
            "samples": ["control", "treatment"],
            "analysis_steps": ["count", "preprocess", "integrate"]
        },
        "multi_omics": {
            "description": "Multi-modal single-cell analysis",
            "samples": ["rna_sample", "atac_sample"],
            "analysis_steps": ["count", "preprocess", "multimodal_integrate"]
        },
        "developmental": {
            "description": "Developmental trajectory analysis",
            "samples": ["timepoint_1", "timepoint_2", "timepoint_3"],
            "analysis_steps": ["count", "preprocess", "trajectory"]
        }
    }
    
    if template_name in templates:
        template = templates[template_name]
        setup_template_project(project_path, template)
```

### Custom Configuration

```python
def custom_configuration():
    """Advanced project configuration options."""
    
    config_options = {
        "analysis": {
            "default_threads": 8,
            "memory_limit": "32G",
            "temp_dir": "/tmp/celline"
        },
        "cluster": {
            "scheduler": "slurm",
            "partition": "compute",
            "time_limit": "24:00:00"
        },
        "output": {
            "compression": True,
            "format": "h5ad",
            "backup": True
        }
    }
    
    return config_options
```

## Best Practices

### Project Organization

1. **Descriptive Names**: Use clear, descriptive project names
2. **Version Control**: Initialize git repository if needed
3. **Documentation**: Include README with project description
4. **Backup Strategy**: Plan for data backup and versioning

### Configuration Management

```python
# Recommended project naming convention
project_names = [
    "2024_brain_atlas_v1",           # Date_tissue_version
    "covid_pbmc_longitudinal",       # Condition_tissue_design
    "development_cortex_timecourse"  # Process_tissue_design
]

# Environment-specific configuration
environments = {
    "development": {
        "debug": True,
        "threads": 2,
        "memory": "8G"
    },
    "production": {
        "debug": False,
        "threads": 16,
        "memory": "64G"
    }
}
```

### Security Considerations

```python
def secure_initialization():
    """Apply security best practices during initialization."""
    
    # Set appropriate file permissions
    os.chmod("setting.toml", 0o600)  # Owner read/write only
    os.chmod("samples.toml", 0o644)  # Standard permissions
    
    # Create .gitignore for sensitive files
    gitignore_content = """
# Sensitive configuration
.env
*.key
*.pem

# Large data files
*.h5
*.h5ad
*.rds

# Temporary files
tmp/
.cache/
"""
    
    with open(".gitignore", "w") as f:
        f.write(gitignore_content)
```

## Related Functions

- [Info](info) - Display project information after initialization
- [Add](add) - Add samples to initialized project
- [Job](job) - Monitor analysis jobs
- All analysis functions - Require proper project initialization

## Troubleshooting

### Common Initialization Issues

1. **Dependency Installation**: Use package managers for missing software
2. **Permission Problems**: Check directory write permissions
3. **R Configuration**: Ensure R is in system PATH
4. **Network Issues**: May affect package installations

### Installation Guides

#### Cell Ranger Installation

```bash
# Download Cell Ranger
curl -o cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=..."

# Extract and add to PATH
tar -xzf cellranger-7.1.0.tar.gz
export PATH=/path/to/cellranger-7.1.0:$PATH
```

#### R Package Installation

```r
# Install required R packages
install.packages(c("Seurat", "scPred", "harmony"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scran")
```

### Manual Initialization

For debugging, manually create project structure:

```bash
# Create directories
mkdir -p data results scripts

# Create basic configuration
cat > setting.toml << EOF
[project]
name = "Manual Project"
version = "1.0.0"
EOF

cat > samples.toml << EOF
# Add samples here
EOF
```

### Recovery from Failed Initialization

```python
def recover_failed_init(project_path):
    """Recover from partially failed initialization."""
    
    # Clean up partial files
    partial_files = ["setting.toml.tmp", "samples.toml.tmp"]
    for file in partial_files:
        if os.path.exists(os.path.join(project_path, file)):
            os.remove(os.path.join(project_path, file))
    
    # Restart initialization
    console.print("[yellow]Cleaning up and retrying initialization...[/yellow]")
    initialize_function = Initialize()
    return initialize_function.call(project)
```