---
title: Functions Reference
navigation: false      # ← 見出し専用にしてリンクは非表示
---

This is a detailed reference for each analysis function in Celline. Each function has its own dedicated page with comprehensive documentation.

## 📚 Function Categories

### 🗄️ Data Management
- [Add](add) - Add accession ID to your project
- [Download](download) - Download data into your project
- [SyncDB](syncdb) - Database synchronization

### 🧮 Data Processing  
- [Count](count) - Count processing with Cell Ranger
- [CreateSeuratObject](createseuratobject) - Create Seurat objects
- [SetTranscriptome](settranscriptome) - Set transcriptome reference

### 🔬 Quality Control & Preprocessing
- [Preprocess](preprocess) - Quality control and filtering

### 📊 Analysis & Visualization
- [PredictCelltype](predict_celltype) - Cell type prediction
- [BuildCellTypeModel](buildcelltypemodel) - Build cell type prediction model
- [Reduce](reduce) - Dimensionality reduction
- [Integrate](integrate) - Data integration
- [BatchCorrection](batch) - Batch effect correction

### 🔧 Utilities
- [Info](info) - System information display
- [Initialize](init) - Project initialization
- [Interactive](interactive) - Launch interactive web interface
- [Job](job) - Job management
- [Bash](bash) - Execute bash commands

---

## Quick Reference

| Command | Class | Description |
|---------|-------|-------------|
| `add` | Add | Add accession ID to your project |
| `download` | Download | Download data into your project |
| `syncdb` | SyncDB | Database synchronization |
| `count` | Count | Count processing with Cell Ranger |
| `createseuratobject` | CreateSeuratObject | Create Seurat objects |
| `settranscriptome` | SetTranscriptome | Set transcriptome reference |
| `preprocess` | Preprocess | Quality control and filtering |
| `predict_celltype` | PredictCelltype | Cell type prediction |
| `buildcelltypemodel` | BuildCellTypeModel | Build cell type prediction model |
| `reduce` | Reduce | Dimensionality reduction |
| `integrate` | Integrate | Data integration |
| `batch` | BatchCorrection | Batch effect correction |
| `info` | Info | System information display |
| `init` | Initialize | Project initialization |
| `interactive` | Interactive | Launch interactive web interface |
| `job` | Job | Job management |
| `bash` | Bash | Execute bash commands |

## Usage Pattern

All functions follow the same pattern:

```python
from celline import Project
from celline.functions.module_name import ClassName

# Create project
project = Project("./my-project")

# Create function instance
function = ClassName(parameters)

# Execute function
result = project.call(function)
```

Click on any function name above to see detailed documentation with parameters, examples, and usage instructions.