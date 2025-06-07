# Developer Guide

This is the developer guide for Celline. It provides detailed explanations of creating custom functions, contributing to the codebase, and plugin development.

## 🎯 Development Overview

Celline is designed with extensibility in mind, and functionality can be extended through the following methods:

- **Creating Custom Functions** - Adding new analysis capabilities
- **Adding Database Handlers** - Supporting new data sources
- **Extending Execution Backends** - Supporting new computing environments
- **Developing Web UI Components** - Extending the interface

## 🛠️ Development Environment Setup

### 1. Getting the Source Code

```bash
# Clone the repository
git clone https://github.com/YUYA556223/celline.git
cd celline

# Checkout the development branch
git checkout develop
```

### 2. Installing Development Dependencies

```bash
# Development with UV environment
uv sync --all-extras

# Or development with pip
pip install -e ".[dev]"

# Install development tools
pip install ruff black pytest pytest-cov sphinx
```

### 3. Verifying Development Environment

```bash
# Run tests
pytest tests/

# Check code formatting
ruff check src/

# Type checking
mypy src/celline/
```

## 🔧 Custom Function Development

### Basic Function Class

```python
from celline.functions._base import CellineFunction
from celline.log.logger import get_logger
import argparse
from typing import Optional, List, TYPE_CHECKING

if TYPE_CHECKING:
    from celline import Project

class MyCustomFunction(CellineFunction):
    """Example custom analysis function"""
    
    def __init__(self, parameter1: str, parameter2: int = 10):
        """
        Initialization method
        
        Args:
            parameter1: Required parameter
            parameter2: Optional parameter
        """
        super().__init__()
        self.parameter1 = parameter1
        self.parameter2 = parameter2
        self.logger = get_logger(__name__)
    
    def call(self, project: "Project") -> "Project":
        """
        Main processing - required implementation method
        
        Args:
            project: Celline project instance
            
        Returns:
            project: Processed project
        """
        self.logger.info(f"Starting custom analysis with {self.parameter1}")
        
        try:
            # Get samples
            samples = self._get_samples_from_project(project)
            
            # Process each sample
            for sample_id in samples:
                self.logger.info(f"Processing sample: {sample_id}")
                self._process_sample(project, sample_id)
            
            self.logger.info("Custom analysis completed successfully")
            
        except Exception as e:
            self.logger.error(f"Error in custom analysis: {e}")
            raise
        
        return project
    
    def _get_samples_from_project(self, project: "Project") -> List[str]:
        """Get sample list from project"""
        import toml
        from celline.config import Config
        
        samples_file = f"{Config.PROJ_ROOT}/samples.toml"
        if not os.path.exists(samples_file):
            return []
        
        with open(samples_file, 'r') as f:
            samples = toml.load(f)
        
        return list(samples.keys())
    
    def _process_sample(self, project: "Project", sample_id: str):
        """Process individual sample"""
        from celline.utils.path import Path
        import scanpy as sc
        
        # Path setup
        path = Path("dummy_project", sample_id)  # Use actual project ID
        
        # Load data
        if path.is_counted:
            adata = sc.read_10x_h5(f"{path.resources_sample_counted}/outs/filtered_feature_bc_matrix.h5")
            
            # Custom processing
            result = self._custom_analysis(adata)
            
            # Save results
            self._save_results(path, result)
    
    def _custom_analysis(self, adata):
        """Custom analysis logic"""
        import numpy as np
        
        # Example: high expression gene detection
        mean_expression = np.mean(adata.X.toarray(), axis=0)
        high_expr_genes = adata.var_names[mean_expression > np.percentile(mean_expression, 90)]
        
        return {
            "high_expression_genes": high_expr_genes.tolist(),
            "total_cells": adata.n_obs,
            "total_genes": adata.n_vars,
            "parameter1_used": self.parameter1,
            "parameter2_used": self.parameter2
        }
    
    def _save_results(self, path: "Path", result: dict):
        """Save results"""
        import json
        
        result_file = f"{path.data_sample}/custom_analysis_results.json"
        with open(result_file, 'w') as f:
            json.dump(result, f, indent=2)
    
    def add_cli_args(self, parser: argparse.ArgumentParser) -> None:
        """Define CLI arguments"""
        parser.add_argument(
            '--parameter1', '-p1',
            required=True,
            help='Required parameter 1'
        )
        parser.add_argument(
            '--parameter2', '-p2',
            type=int,
            default=10,
            help='Optional parameter 2 (default: 10)'
        )
    
    def cli(self, project: "Project", args: Optional[argparse.Namespace] = None) -> "Project":
        """CLI execution entry point"""
        if args:
            self.parameter1 = args.parameter1
            self.parameter2 = args.parameter2
        
        return self.call(project)
    
    def get_description(self) -> str:
        """Function description"""
        return """Custom analysis function example.
        
        This function demonstrates how to create custom analysis
        functions that integrate with the Celline framework."""
    
    def get_usage_examples(self) -> List[str]:
        """Usage examples"""
        return [
            "celline run my_custom --parameter1 value1",
            "celline run my_custom --parameter1 value1 --parameter2 20"
        ]
```

### Function Registration

```python
# Register function to Celline
from celline.cli.registry import get_registry

def register_custom_functions():
    """Register custom functions"""
    registry = get_registry()
    
    registry.register_function(
        name="my_custom",
        class_ref=MyCustomFunction,
        module_path="my_package.custom_functions"
    )

# Execute during package initialization
register_custom_functions()
```

### Creating Tests

```python
# tests/test_custom_function.py
import pytest
import tempfile
import os
from celline import Project
from my_package.custom_functions import MyCustomFunction

class TestMyCustomFunction:
    
    @pytest.fixture
    def temp_project(self):
        """Temporary project for testing"""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create configuration file
            setting_content = """
[project]
name = "test_project"
version = "1.0.0"

[execution]
system = "multithreading"
nthread = 1

[R]
r_path = "/usr/bin/R"

[fetch]
wait_time = 1
"""
            with open(f"{temp_dir}/setting.toml", 'w') as f:
                f.write(setting_content)
            
            # Create sample file
            with open(f"{temp_dir}/samples.toml", 'w') as f:
                f.write('GSM123456 = "Test Sample"\n')
            
            yield Project(temp_dir, "test_project")
    
    def test_custom_function_creation(self):
        """Test custom function creation"""
        func = MyCustomFunction("test_param", 5)
        assert func.parameter1 == "test_param"
        assert func.parameter2 == 5
    
    def test_custom_function_execution(self, temp_project):
        """Test custom function execution"""
        func = MyCustomFunction("test_param", 5)
        
        # Prepare mock data if needed
        # ...
        
        result_project = func.call(temp_project)
        assert result_project is not None
    
    def test_cli_args_parsing(self):
        """Test CLI argument parsing"""
        import argparse
        
        func = MyCustomFunction("dummy", 1)
        parser = argparse.ArgumentParser()
        func.add_cli_args(parser)
        
        args = parser.parse_args(['--parameter1', 'test_value', '--parameter2', '15'])
        assert args.parameter1 == 'test_value'
        assert args.parameter2 == 15
```

## 🗄️ Database Handler Development

### Custom Database Handler

```python
from celline.DB.dev.handler import DatabaseHandler
from celline.DB.dev.model import SampleSchema
from typing import Optional
import requests

class CustomDatabaseHandler(DatabaseHandler):
    """Example custom database handler"""
    
    def __init__(self):
        self.api_base = "https://api.customdb.org/v1"
        self.headers = {"Accept": "application/json"}
    
    def can_handle(self, sample_id: str) -> bool:
        """Check if this handler can process the sample ID"""
        return sample_id.startswith("CUSTOM")
    
    def fetch_sample_metadata(self, sample_id: str) -> Optional[SampleSchema]:
        """Fetch sample metadata"""
        try:
            response = requests.get(
                f"{self.api_base}/samples/{sample_id}",
                headers=self.headers,
                timeout=30
            )
            response.raise_for_status()
            
            data = response.json()
            
            # Convert to SampleSchema
            schema = SampleSchema(
                key=data["id"],
                title=data["title"],
                organism=data["organism"],
                library_strategy=data["library_strategy"],
                parent=data["study_id"],
                children=data.get("run_ids", ""),
                # Other required fields
            )
            
            return schema
            
        except requests.RequestException as e:
            self.logger.error(f"Failed to fetch metadata for {sample_id}: {e}")
            return None
    
    def add(self, sample_id: str) -> bool:
        """Add sample to local database"""
        metadata = self.fetch_sample_metadata(sample_id)
        if metadata:
            # Save to local database
            self._save_to_local_db(metadata)
            return True
        return False
    
    def _save_to_local_db(self, schema: SampleSchema):
        """Save to local database"""
        import polars as pl
        from celline.config import Config
        
        # Load existing data
        db_file = f"{Config.EXEC_ROOT}/DB/CUSTOM_SAMPLES.parquet"
        
        if os.path.exists(db_file):
            df = pl.read_parquet(db_file)
        else:
            df = pl.DataFrame()
        
        # Add new record
        new_record = pl.DataFrame([{
            "key": schema.key,
            "title": schema.title,
            "organism": schema.organism,
            "library_strategy": schema.library_strategy,
            "parent": schema.parent,
            "children": schema.children
        }])
        
        # Combine data and save
        updated_df = pl.concat([df, new_record]) if not df.is_empty() else new_record
        updated_df.write_parquet(db_file)

# Register handler
from celline.DB.dev.handler import HandleResolver

def register_custom_handler():
    """Register custom handler"""
    custom_handler = CustomDatabaseHandler()
    HandleResolver.register_handler(custom_handler)

register_custom_handler()
```

## 🎨 Web UI Component Development

### Vue.js Component

```vue
<!-- CustomAnalysisPanel.vue -->
<template>
  <div class="custom-analysis-panel">
    <h3>Custom Analysis</h3>
    
    <div class="parameter-section">
      <h4>Parameters</h4>
      
      <div class="form-group">
        <label for="parameter1">Parameter 1:</label>
        <input 
          id="parameter1"
          v-model="parameters.parameter1" 
          type="text"
          placeholder="Enter parameter 1"
          required
        />
      </div>
      
      <div class="form-group">
        <label for="parameter2">Parameter 2:</label>
        <input 
          id="parameter2"
          v-model="parameters.parameter2" 
          type="number"
          min="1"
          max="100"
        />
      </div>
    </div>
    
    <div class="action-section">
      <button 
        @click="executeAnalysis" 
        :disabled="isRunning || !canExecute"
        class="btn-primary"
      >
        {{ isRunning ? 'Running...' : 'Execute Analysis' }}
      </button>
    </div>
    
    <div v-if="result" class="result-section">
      <h4>Results</h4>
      <div class="result-content">
        <p>High expression genes found: {{ result.high_expression_genes.length }}</p>
        <p>Total cells: {{ result.total_cells }}</p>
        <p>Total genes: {{ result.total_genes }}</p>
      </div>
    </div>
    
    <div v-if="error" class="error-section">
      <p class="error-message">{{ error }}</p>
    </div>
  </div>
</template>

<script>
export default {
  name: 'CustomAnalysisPanel',
  
  data() {
    return {
      parameters: {
        parameter1: '',
        parameter2: 10
      },
      isRunning: false,
      result: null,
      error: null
    };
  },
  
  computed: {
    canExecute() {
      return this.parameters.parameter1.trim() !== '';
    }
  },
  
  methods: {
    async executeAnalysis() {
      this.isRunning = true;
      this.error = null;
      
      try {
        const response = await fetch('/api/functions/my_custom/execute', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json'
          },
          body: JSON.stringify(this.parameters)
        });
        
        if (!response.ok) {
          throw new Error(`HTTP error! status: ${response.status}`);
        }
        
        const data = await response.json();
        
        if (data.job_id) {
          // Monitor job
          await this.monitorJob(data.job_id);
        } else {
          this.result = data.result;
        }
        
      } catch (error) {
        this.error = `Analysis failed: ${error.message}`;
      } finally {
        this.isRunning = false;
      }
    },
    
    async monitorJob(jobId) {
      const checkStatus = async () => {
        try {
          const response = await fetch(`/api/jobs/${jobId}`);
          const status = await response.json();
          
          if (status.status === 'completed') {
            this.result = status.result;
            return true;
          } else if (status.status === 'failed') {
            this.error = status.message;
            return true;
          }
          
          return false;
        } catch (error) {
          this.error = `Job monitoring failed: ${error.message}`;
          return true;
        }
      };
      
      // Polling for job status
      while (!(await checkStatus())) {
        await new Promise(resolve => setTimeout(resolve, 2000));
      }
    }
  }
};
</script>

<style scoped>
.custom-analysis-panel {
  padding: 20px;
  border: 1px solid #ddd;
  border-radius: 8px;
  margin: 10px 0;
}

.form-group {
  margin: 10px 0;
}

.form-group label {
  display: block;
  margin-bottom: 5px;
  font-weight: bold;
}

.form-group input {
  width: 100%;
  padding: 8px;
  border: 1px solid #ccc;
  border-radius: 4px;
}

.btn-primary {
  background-color: #007bff;
  color: white;
  border: none;
  padding: 10px 20px;
  border-radius: 4px;
  cursor: pointer;
  font-size: 16px;
}

.btn-primary:disabled {
  background-color: #6c757d;
  cursor: not-allowed;
}

.result-section {
  margin-top: 20px;
  padding: 15px;
  background-color: #f8f9fa;
  border-radius: 4px;
}

.error-section {
  margin-top: 20px;
  padding: 15px;
  background-color: #f8d7da;
  border: 1px solid #f5c6cb;
  border-radius: 4px;
}

.error-message {
  color: #721c24;
  margin: 0;
}
</style>
```

### Adding FastAPI Endpoints

```python
# celline/api/custom_endpoints.py
from fastapi import APIRouter, HTTPException, BackgroundTasks
from pydantic import BaseModel
from typing import Dict, Any
import uuid
from datetime import datetime

from celline.functions.my_custom import MyCustomFunction
from celline import Project

router = APIRouter(prefix="/api/functions/my_custom", tags=["custom"])

class CustomAnalysisRequest(BaseModel):
    parameter1: str
    parameter2: int = 10

class CustomAnalysisResponse(BaseModel):
    job_id: str
    status: str

# Job storage (use Redis etc. in actual applications)
active_jobs: Dict[str, Dict[str, Any]] = {}

@router.post("/execute", response_model=CustomAnalysisResponse)
async def execute_custom_analysis(
    request: CustomAnalysisRequest,
    background_tasks: BackgroundTasks
):
    """Execute custom analysis"""
    job_id = str(uuid.uuid4())
    
    # Initialize job information
    active_jobs[job_id] = {
        "status": "pending",
        "created_at": datetime.now(),
        "progress": 0.0,
        "message": "Analysis queued",
        "result": None
    }
    
    # Execute in background task
    background_tasks.add_task(
        run_custom_analysis,
        job_id,
        request.parameter1,
        request.parameter2
    )
    
    return CustomAnalysisResponse(job_id=job_id, status="started")

async def run_custom_analysis(job_id: str, param1: str, param2: int):
    """Execute custom analysis in background"""
    try:
        # Update job status
        active_jobs[job_id]["status"] = "running"
        active_jobs[job_id]["message"] = "Starting custom analysis"
        active_jobs[job_id]["progress"] = 10.0
        
        # Get project
        project = Project("./")  # Current project
        
        # Execute custom function
        custom_func = MyCustomFunction(param1, param2)
        
        active_jobs[job_id]["progress"] = 50.0
        active_jobs[job_id]["message"] = "Running analysis"
        
        # Execute analysis (if run synchronously)
        result_project = custom_func.call(project)
        
        # Load results
        result = load_analysis_result(project)
        
        # Complete
        active_jobs[job_id]["status"] = "completed"
        active_jobs[job_id]["progress"] = 100.0
        active_jobs[job_id]["message"] = "Analysis completed successfully"
        active_jobs[job_id]["result"] = result
        
    except Exception as e:
        active_jobs[job_id]["status"] = "failed"
        active_jobs[job_id]["message"] = f"Analysis failed: {str(e)}"

def load_analysis_result(project: Project) -> Dict[str, Any]:
    """Load analysis results"""
    import json
    import os
    from celline.config import Config
    
    # Load result file
    result_file = f"{Config.PROJ_ROOT}/data/*/custom_analysis_results.json"
    
    # Use proper file path in actual implementation
    # Simplified here
    return {
        "high_expression_genes": ["GENE1", "GENE2", "GENE3"],
        "total_cells": 1000,
        "total_genes": 20000
    }

@router.get("/status/{job_id}")
async def get_job_status(job_id: str):
    """Get job status"""
    if job_id not in active_jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    return active_jobs[job_id]

# Add router to main app
# celline/api/main.py
from .custom_endpoints import router as custom_router

app.include_router(custom_router)
```

## 📊 Testing and Debugging

### Writing Unit Tests

```python
# tests/test_custom_analysis.py
import pytest
import tempfile
import json
from unittest.mock import Mock, patch
from celline import Project
from my_package.custom_functions import MyCustomFunction

class TestCustomAnalysisIntegration:
    
    @pytest.fixture
    def mock_scanpy_data(self):
        """Mock of Scanpy AnnData object"""
        import numpy as np
        
        mock_adata = Mock()
        mock_adata.n_obs = 1000
        mock_adata.n_vars = 20000
        mock_adata.X.toarray.return_value = np.random.rand(1000, 20000)
        mock_adata.var_names = [f"GENE_{i}" for i in range(20000)]
        
        return mock_adata
    
    @patch('scanpy.read_10x_h5')
    def test_full_analysis_workflow(self, mock_read_h5, mock_scanpy_data):
        """Test full analysis workflow"""
        # Setup mock
        mock_read_h5.return_value = mock_scanpy_data
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Prepare test environment
            project = self._setup_test_project(temp_dir)
            
            # Execute custom function
            func = MyCustomFunction("test_param", 5)
            result_project = func.call(project)
            
            # Verify results
            assert result_project is not None
            
            # Check result files
            result_files = self._find_result_files(temp_dir)
            assert len(result_files) > 0
            
            # Verify result content
            with open(result_files[0], 'r') as f:
                result = json.load(f)
            
            assert "high_expression_genes" in result
            assert "total_cells" in result
            assert result["parameter1_used"] == "test_param"
            assert result["parameter2_used"] == 5
    
    def _setup_test_project(self, temp_dir):
        """Setup test project"""
        # Create configuration file
        setting_content = """
[project]
name = "test_project"
version = "1.0.0"

[execution]
system = "multithreading"
nthread = 1

[R]
r_path = "/usr/bin/R"

[fetch]
wait_time = 1
"""
        with open(f"{temp_dir}/setting.toml", 'w') as f:
            f.write(setting_content)
        
        # Create sample file
        with open(f"{temp_dir}/samples.toml", 'w') as f:
            f.write('GSM123456 = "Test Sample"\n')
        
        # Create data directory
        data_dir = f"{temp_dir}/data/GSM123456"
        os.makedirs(data_dir, exist_ok=True)
        
        # Create resource directory
        resource_dir = f"{temp_dir}/resources/GSM123456/counted/outs"
        os.makedirs(resource_dir, exist_ok=True)
        
        return Project(temp_dir, "test_project")
    
    def _find_result_files(self, temp_dir):
        """Find result files"""
        import glob
        return glob.glob(f"{temp_dir}/data/*/custom_analysis_results.json")
```

### Debug Logging

```python
# Setup detailed logging for debugging
import logging
from celline.log.logger import get_logger

def setup_debug_logging():
    """Setup debug logging"""
    
    # Root logger configuration
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('celline_debug.log'),
            logging.StreamHandler()
        ]
    )
    
    # Celline specific logger configuration
    celline_logger = get_logger("celline")
    celline_logger.setLevel(logging.DEBUG)
    
    # Adjust external library log levels
    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)

# Debug mode execution
if __name__ == "__main__":
    setup_debug_logging()
    
    # Test execution of custom function
    logger = get_logger(__name__)
    logger.debug("Starting debug session")
    
    try:
        project = Project("./test_project")
        func = MyCustomFunction("debug_param", 1)
        func.call(project)
        
    except Exception as e:
        logger.exception("Error occurred during debug session")
        raise
```

## 🚀 Contributing

### Creating Pull Requests

```bash
# 1. Clone forked repository
git clone https://github.com/yourusername/celline.git
cd celline

# 2. Create new branch
git checkout -b feature/my-new-feature

# 3. Make changes
# Code changes, add tests, etc.

# 4. Run tests
pytest tests/

# 5. Code formatting
ruff format src/
ruff check src/ --fix

# 6. Commit
git add .
git commit -m "Add new custom analysis feature"

# 7. Push
git push origin feature/my-new-feature

# 8. Create pull request on GitHub
```

### Code Quality Standards

```python
# Coding convention examples

# 1. Use type hints
def process_data(data: List[Dict[str, Any]]) -> pd.DataFrame:
    """Use clear type hints"""
    pass

# 2. Documentation
def complex_function(param1: str, param2: int) -> bool:
    """
    Add detailed documentation for complex functions
    
    Args:
        param1: Description of parameter 1
        param2: Description of parameter 2
    
    Returns:
        Boolean result of processing
    
    Raises:
        ValueError: When invalid parameters are provided
    
    Example:
        >>> result = complex_function("test", 42)
        >>> assert result is True
    """
    pass

# 3. Error handling
def safe_function():
    """Proper error handling"""
    try:
        # Processing
        pass
    except SpecificException as e:
        logger.error(f"Specific error occurred: {e}")
        raise
    except Exception as e:
        logger.exception("Unexpected error occurred")
        raise RuntimeError(f"Function failed: {e}") from e
```

### Documentation

```python
# Docstring writing examples
class NewAnalysisFunction(CellineFunction):
    """
    New analysis function class
    
    This function performs XX analysis and generates YY results.
    
    Attributes:
        param1: Analysis parameter 1
        param2: Analysis parameter 2
    
    Example:
        >>> func = NewAnalysisFunction("value1", 10)
        >>> result = func.call(project)
        >>> print(result.status)
    """
    
    def __init__(self, param1: str, param2: int):
        """
        Initialization
        
        Args:
            param1: Required parameter
            param2: Optional parameter
        """
        pass
    
    def call(self, project: "Project") -> "Project":
        """
        Main processing
        
        Args:
            project: Celline project
            
        Returns:
            Processed project
            
        Raises:
            CellineException: When processing fails
        """
        pass
```

## 🔧 Debugging and Profiling

### Performance Measurement

```python
import cProfile
import pstats
from functools import wraps

def profile_function(func):
    """Function profiling decorator"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        profiler = cProfile.Profile()
        profiler.enable()
        
        try:
            result = func(*args, **kwargs)
        finally:
            profiler.disable()
            
            # Output results
            stats = pstats.Stats(profiler)
            stats.sort_stats('cumulative')
            stats.print_stats(20)  # Top 20
            
        return result
    return wrapper

# Usage example
@profile_function
def heavy_analysis_function(data):
    """Heavy processing function"""
    # Processing content
    pass
```

### Memory Usage Monitoring

```python
import psutil
import os
from contextlib import contextmanager

@contextmanager
def memory_monitor(description: str = ""):
    """Monitor memory usage"""
    process = psutil.Process(os.getpid())
    start_memory = process.memory_info().rss / 1024 / 1024  # MB
    
    print(f"Memory before {description}: {start_memory:.2f} MB")
    
    try:
        yield
    finally:
        end_memory = process.memory_info().rss / 1024 / 1024  # MB
        diff = end_memory - start_memory
        print(f"Memory after {description}: {end_memory:.2f} MB (diff: {diff:+.2f} MB)")

# Usage example
with memory_monitor("custom analysis"):
    result = heavy_function(large_data)
```

---

> **Success**: We welcome contributions to the Celline ecosystem using this developer guide!
> If you have questions or suggestions, please use GitHub Issues or Discussions.