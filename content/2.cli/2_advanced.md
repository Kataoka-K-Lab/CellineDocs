# Advanced Usage

This document explains advanced features and usage methods of Celline.
It provides detailed explanations of advanced functions such as custom workflows, large-scale data processing, and integrated analysis.

## 🚀 Large-Scale Data Processing

### Execution in Cluster Environments

#### PBS/Torque Cluster

```bash
# PBS system configuration
celline config --system PBS --pbs-server my-cluster

# Job resource configuration
cat > pbs_config.toml << EOF
[pbs]
queue = "normal"
walltime = "48:00:00"
memory = "128GB"
ncpus = 32
nodes = 1
EOF

# Large dataset processing
celline run count --config pbs_config.toml
```

#### Slurm Cluster (Custom Implementation)

```bash
#!/bin/bash
#SBATCH --job-name=celline-analysis
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16

module load R/4.3.0 cellranger/7.0.0

# Process multiple samples in parallel
celline run download --nthread 16
celline run count --nthread 16
```

### Memory-Efficient Processing

```python
# Memory management for large datasets
from celline import Project
from celline.utils.memory import MemoryManager

project = Project("./large-dataset")

# Limit memory usage
with MemoryManager(max_memory="32GB") as mem:
    # Process in chunks
    for chunk in project.iter_samples(chunk_size=10):
        project.process_chunk(chunk)
        mem.clear_cache()
```

## 🔬 Custom Workflows

### Pipeline Design

```python
from celline import Project, Pipeline
from celline.functions import *

# Custom pipeline definition
class CustomSingleCellPipeline(Pipeline):
    def __init__(self, project: Project):
        super().__init__(project)
        
    def run(self):
        # Step-by-step processing
        self.stage1_data_acquisition()
        self.stage2_quality_control()
        self.stage3_analysis()
        self.stage4_integration()
    
    def stage1_data_acquisition(self):
        """Data acquisition stage"""
        self.project.call(Add(self.sample_list))
        self.project.call(Download())
        
    def stage2_quality_control(self):
        """Quality control stage"""
        self.project.call(Count())
        self.project.call(Preprocess())
        
    def stage3_analysis(self):
        """Analysis stage"""
        self.project.call(PredictCelltype())
        self.project.call(Reduce())
        
    def stage4_integration(self):
        """Integration stage"""
        self.project.call(Integrate())

# Pipeline execution
project = Project("./my-analysis")
pipeline = CustomSingleCellPipeline(project)
pipeline.run()
```

### Conditional Processing

```python
from celline.functions import *

project = Project("./conditional-analysis")

# Processing branching based on conditions
project.call_if_else(
    condition=lambda p: p.has_spatial_data(),
    true=SpatialAnalysis(),
    false=StandardAnalysis()
)

# Processing based on data type
if project.is_10x_data():
    project.call(Count())
elif project.is_smartseq_data():
    project.call(SmartSeqProcessing())
else:
    project.call(CustomProcessing())
```

## 📊 Batch Processing and Automation

### Batch Processing of Multiple Projects

```bash
#!/bin/bash
# batch_analysis.sh

# Project list file
PROJECTS="project_list.txt"

while IFS=',' read -r project_name gse_id description; do
    echo "Processing $project_name ($gse_id)"
    
    # Create project directory
    mkdir -p "$project_name"
    cd "$project_name"
    
    # Execute analysis pipeline
    celline init "$project_name"
    celline run add "$gse_id"
    celline run download --nthread 8
    celline run count
    celline run preprocess
    celline run predict_celltype
    
    # Archive results
    tar -czf "../${project_name}_results.tar.gz" results/
    
    cd ..
done < "$PROJECTS"
```

### Automation Using Makefile

```makefile
# Makefile for Celline analysis

# Variable definitions
PROJECT_NAME := my-scrna-analysis
GSE_ID := GSE123456
NTHREAD := 8

# Default target
all: init download count preprocess analyze

# Project initialization
init:
	celline init $(PROJECT_NAME)

# Sample addition
add:
	celline run add $(GSE_ID)

# Data download
download: add
	celline run download --nthread $(NTHREAD)

# Count processing
count: download
	celline run count

# Preprocessing
preprocess: count
	celline run preprocess

# Analysis
analyze: preprocess
	celline run predict_celltype
	celline run reduce
	celline run integrate

# Check results
check:
	celline info
	find results/ -name "*.csv" -o -name "*.png" | head -10

# Cleanup
clean:
	rm -rf resources/*/raw resources/*/tmp
	
# Complete cleanup
distclean:
	rm -rf resources/ data/ results/

.PHONY: all init add download count preprocess analyze check clean distclean
```

## 🔄 Integrated Analysis

### Multi-Dataset Integration

```python
from celline.functions.integrate import MultiDatasetIntegration

# Multi-dataset integration analysis
integration = MultiDatasetIntegration([
    "dataset1_GSE123456",
    "dataset2_GSE789012", 
    "dataset3_GSE345678"
])

project = Project("./integrated-analysis")
project.call(integration)

# Comparison of batch effect correction methods
methods = ["harmony", "combat", "mnn", "cca"]
for method in methods:
    integration.method = method
    project.call(integration)
    project.save_results(f"integration_{method}")
```

### Meta-Analysis

```python
from celline.functions.meta import MetaAnalysis

# Meta-analysis execution
meta = MetaAnalysis([
    {"gse": "GSE123456", "condition": "control", "tissue": "brain"},
    {"gse": "GSE789012", "condition": "disease", "tissue": "brain"},
    {"gse": "GSE345678", "condition": "treatment", "tissue": "brain"}
])

project = Project("./meta-analysis")
project.call(meta)

# Statistical analysis of results
meta.perform_differential_analysis()
meta.generate_forest_plots()
```

## 🧬 Custom Analysis Functions

### Creating New Analysis Functions

```python
from celline.functions._base import CellineFunction
import argparse

class CustomCellTypeClassification(CellineFunction):
    """Custom cell type classification function"""
    
    def __init__(self, model_path: str, confidence_threshold: float = 0.8):
        super().__init__()
        self.model_path = model_path
        self.confidence_threshold = confidence_threshold
    
    def call(self, project):
        """Main processing"""
        import pickle
        import scanpy as sc
        
        # Load custom model
        with open(self.model_path, 'rb') as f:
            model = pickle.load(f)
        
        # Process each sample
        for sample_id in project.get_samples():
            adata = project.load_sample_data(sample_id)
            
            # Feature extraction
            features = self.extract_features(adata)
            
            # Execute prediction
            predictions = model.predict(features)
            confidence = model.predict_proba(features)
            
            # Confidence filtering
            high_conf_mask = confidence.max(axis=1) >= self.confidence_threshold
            
            # Save results
            adata.obs['custom_celltype'] = predictions
            adata.obs['celltype_confidence'] = confidence.max(axis=1)
            adata.obs['high_confidence'] = high_conf_mask
            
            project.save_sample_data(sample_id, adata)
        
        return project
    
    def extract_features(self, adata):
        """Feature extraction"""
        # Custom feature calculation
        import numpy as np
        
        # Highly variable gene expression
        hvg_expr = adata[:, adata.var.highly_variable].X.toarray()
        
        # Pathway scores
        pathway_scores = self.calculate_pathway_scores(adata)
        
        # Concatenation
        features = np.concatenate([hvg_expr, pathway_scores], axis=1)
        return features
    
    def add_cli_args(self, parser: argparse.ArgumentParser):
        """Add CLI arguments"""
        parser.add_argument('--model-path', required=True,
                          help='Path to the trained classification model')
        parser.add_argument('--confidence-threshold', type=float, default=0.8,
                          help='Confidence threshold for predictions')
    
    def get_description(self):
        return "Custom cell type classification using trained models"
```

### Registering Analysis Functions

```python
# Register custom function to Celline
from celline.cli.registry import get_registry

registry = get_registry()
registry.register_function(
    name="custom_classify",
    class_ref=CustomCellTypeClassification,
    module_path="my_custom_functions.classification"
)

# Use via CLI
# celline run custom_classify --model-path my_model.pkl
```

## 📊 Advanced Data Visualization Features

### Custom Plot Creation

```python
from celline.visualization import CellinePlotter

class AdvancedVisualization(CellineFunction):
    def call(self, project):
        plotter = CellinePlotter(project)
        
        # Integrated UMAP for multiple samples
        plotter.integrated_umap(
            samples=project.get_samples(),
            color_by="celltype",
            split_by="condition"
        )
        
        # Expression heatmap
        plotter.expression_heatmap(
            genes=["CD4", "CD8A", "IL2", "IFNG"],
            group_by="celltype",
            save_path="results/expression_heatmap.pdf"
        )
        
        # Trajectory analysis
        plotter.trajectory_plot(
            root_cell="stem_cell",
            save_path="results/trajectory.pdf"
        )
        
        return project
```

### Interactive Visualization

```python
import plotly.graph_objects as go
from celline.interactive import InteractivePlotter

def create_interactive_plots(project):
    plotter = InteractivePlotter(project)
    
    # 3D UMAP
    fig = plotter.plot_3d_umap(
        color_by="celltype",
        hover_data=["sample", "condition", "batch"]
    )
    
    # Interactive gene expression
    expr_fig = plotter.plot_gene_expression(
        genes=["CD4", "CD8A"],
        plot_type="violin"
    )
    
    # Dashboard creation
    dashboard = plotter.create_dashboard([fig, expr_fig])
    dashboard.serve(port=8050)
```

## 🔌 API Integration

### External Database Integration

```python
from celline.database import ExternalDBConnector

class CustomDBIntegration(CellineFunction):
    def __init__(self, db_config: dict):
        self.db_config = db_config
    
    def call(self, project):
        # Fetch additional information from external database
        connector = ExternalDBConnector(self.db_config)
        
        for sample_id in project.get_samples():
            # Extend metadata
            extended_meta = connector.fetch_extended_metadata(sample_id)
            
            # Update gene annotations
            gene_annotations = connector.fetch_gene_annotations(
                species=project.get_species()
            )
            
            # Integrate into project
            project.update_sample_metadata(sample_id, extended_meta)
            project.update_gene_annotations(gene_annotations)
        
        return project
```

### Machine Learning Pipeline Integration

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier

class MLPipeline(CellineFunction):
    def __init__(self, model_config: dict):
        self.model_config = model_config
    
    def call(self, project):
        # Build machine learning pipeline
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', RandomForestClassifier(**self.model_config))
        ])
        
        # Prepare training data
        X_train, y_train = project.prepare_training_data()
        
        # Train model
        pipeline.fit(X_train, y_train)
        
        # Execute predictions
        for sample_id in project.get_samples():
            X_test = project.prepare_test_data(sample_id)
            predictions = pipeline.predict(X_test)
            project.save_predictions(sample_id, predictions)
        
        # Save model
        project.save_model(pipeline, "trained_classifier.pkl")
        
        return project
```

## 🔧 Performance Optimization

### Parallel Processing Optimization

```python
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from celline.utils import ProcessManager

class OptimizedProcessing(CellineFunction):
    def call(self, project):
        samples = project.get_samples()
        
        # CPU-intensive tasks use process parallelization
        with ProcessPoolExecutor(max_workers=8) as executor:
            cpu_futures = [
                executor.submit(self.cpu_intensive_task, sample)
                for sample in samples
            ]
        
        # I/O-intensive tasks use thread parallelization
        with ThreadPoolExecutor(max_workers=16) as executor:
            io_futures = [
                executor.submit(self.io_intensive_task, sample)
                for sample in samples
            ]
        
        # Collect results
        cpu_results = [f.result() for f in cpu_futures]
        io_results = [f.result() for f in io_futures]
        
        return project
```

### Memory Usage Optimization

```python
from celline.utils.memory import MemoryProfiler, LazyLoader

class MemoryOptimizedAnalysis(CellineFunction):
    def call(self, project):
        with MemoryProfiler() as profiler:
            # Use lazy loading
            lazy_data = LazyLoader(project.data_path)
            
            # Process in chunks
            for chunk in lazy_data.iter_chunks(chunk_size=1000):
                self.process_chunk(chunk)
                
                # Monitor memory usage
                if profiler.get_memory_usage() > profiler.memory_limit:
                    profiler.clear_cache()
                    gc.collect()
        
        return project
```

## 🚨 Error Handling and Recovery

### Robust Error Handling

```python
from celline.utils.exceptions import CellineException
from celline.utils.recovery import RecoveryManager

class RobustAnalysis(CellineFunction):
    def call(self, project):
        recovery = RecoveryManager(project)
        
        try:
            # Execute analysis steps
            for step in self.analysis_steps:
                checkpoint = recovery.create_checkpoint(step.name)
                
                try:
                    step.execute(project)
                    recovery.mark_success(checkpoint)
                except Exception as e:
                    recovery.mark_failure(checkpoint, e)
                    
                    # Attempt automatic recovery
                    if recovery.can_recover(step):
                        recovery.recover_from_checkpoint(checkpoint)
                    else:
                        raise CellineException(f"Step {step.name} failed: {e}")
        
        except CellineException as e:
            # Generate detailed error report
            recovery.generate_error_report(e)
            raise
        
        return project
```

---

> **Warning**: When using advanced features, conduct thorough testing and verify operation with small datasets before running in production environments.