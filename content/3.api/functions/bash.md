# Bash Function

Execute shell commands and scripts with promise-based asynchronous handling.

## Overview

The `Bash` function provides a powerful interface for executing shell commands within Celline workflows. It features a Promise-based execution model with success and error callbacks, enabling seamless integration of custom shell scripts and system commands into analysis pipelines.

## Class Information

- **Module**: `celline.functions.bash`
- **Class**: `Bash`
- **Base Class**: `CellineFunction`

## Related Classes

### `Promise`

Asynchronous execution handler for shell commands:

```python
class Promise:
    def __init__(self, command: str) -> None
    def execute(self) -> "Promise"
    def then(self, callback: Callable[[str], None]) -> "Promise"
    def error(self, callback: Callable[[str], None]) -> "Promise"
    def is_done(self) -> bool
```

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `cmd` | `str` | Yes | Shell command to execute |
| `then` | `Optional[Callable[[str], None]]` | No | Success callback function |
| `catch` | `Optional[Callable[[str], None]]` | No | Error callback function |

### Promise Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `command` | `str` | Yes | Shell command string |

## Usage Examples

### Python API

#### Basic Command Execution

```python
from celline import Project
from celline.functions.bash import Bash

# Create project
project = Project("./my-project")

# Execute simple command
bash_function = Bash(
    cmd="echo 'Hello from Celline'",
    then=None,
    catch=None
)

# Execute function
result = project.call(bash_function)
```

#### Command with Callbacks

```python
from celline import Project
from celline.functions.bash import Bash

def on_success(output: str):
    print(f"Command succeeded: {output.strip()}")

def on_error(error: str):
    print(f"Command failed: {error.strip()}")

# Create project
project = Project("./my-project")

# Execute with callbacks
bash_function = Bash(
    cmd="ls -la data/",
    then=on_success,
    catch=on_error
)

# Execute function
result = project.call(bash_function)
```

#### File Processing Commands

```python
from celline import Project
from celline.functions.bash import Bash

def process_success(output: str):
    lines = output.strip().split('\n')
    print(f"Processed {len(lines)} files")

def process_error(error: str):
    print(f"Processing failed: {error}")

# Create project
project = Project("./my-project")

# Process files with find command
bash_function = Bash(
    cmd="find data/ -name '*.h5' -type f",
    then=process_success,
    catch=process_error
)

result = project.call(bash_function)
```

#### Complex Pipeline Commands

```python
from celline import Project
from celline.functions.bash import Bash

# Multi-step data processing
def pipeline_success(output: str):
    print("Data pipeline completed successfully")
    print(f"Output: {output}")

def pipeline_error(error: str):
    print(f"Pipeline failed at: {error}")

# Create project
project = Project("./my-project")

# Complex shell pipeline
bash_function = Bash(
    cmd="""
    cd data && \
    find . -name '*.tsv' | \
    xargs wc -l | \
    sort -nr | \
    head -10
    """,
    then=pipeline_success,
    catch=pipeline_error
)

result = project.call(bash_function)
```

### Direct Promise Usage

#### Promise Chain Execution

```python
from celline.functions.bash import Promise

# Create and execute promise
promise = Promise("python --version")

# Chain operations
result = promise.then(
    lambda output: print(f"Python version: {output.strip()}")
).error(
    lambda error: print(f"Error: {error.strip()}")
).execute()

# Check completion
if promise.is_done():
    print("Command completed")
```

#### Multiple Command Execution

```python
from celline.functions.bash import Promise

commands = [
    "python --version",
    "R --version", 
    "cellranger --version"
]

def check_version(cmd: str):
    def success_handler(output: str):
        print(f"{cmd}: {output.strip()}")
    
    def error_handler(error: str):
        print(f"{cmd}: Not available")
    
    Promise(cmd).then(success_handler).error(error_handler).execute()

# Check all versions
for cmd in commands:
    check_version(cmd)
```

## CLI Usage

The Bash function is primarily designed for programmatic use. For direct CLI execution, use standard shell commands or create custom scripts.

### Integration with CLI Functions

```python
# Custom CLI function using Bash
class CustomScript(CellineFunction):
    def __init__(self, script_path: str):
        self.script_path = script_path
    
    def call(self, project):
        def on_success(output):
            console.print(f"[green]Script completed successfully[/green]")
            console.print(output)
        
        def on_error(error):
            console.print(f"[red]Script failed[/red]")
            console.print(error)
        
        bash_function = Bash(
            cmd=f"bash {self.script_path}",
            then=on_success,
            catch=on_error
        )
        
        return bash_function.call(project)
```

## Implementation Details

### Promise-Based Execution

The function uses a Promise pattern for asynchronous command handling:

```python
class Promise:
    def __init__(self, command: str) -> None:
        self._is_done: bool = False
        self._result: Optional[str] = None
        self._error: Optional[str] = None
        self._callbacks: List[Callable[[str], None]] = []
        self._error_callbacks: List[Callable[[str], None]] = []
        self._command: str = command
    
    def execute(self) -> "Promise":
        try:
            result: bytes = subprocess.check_output(
                self._command,
                shell=True,
                stderr=subprocess.STDOUT,
                executable=os.environ.get("SHELL"),
            )
            self.resolve(result.decode())
        except subprocess.CalledProcessError as e:
            self.catch_error(e.output.decode())
        return self
```

### Shell Environment

Commands execute in the system shell environment:

- **Shell Detection**: Uses `$SHELL` environment variable
- **Path Resolution**: Inherits system `$PATH`
- **Environment Variables**: Access to all environment variables
- **Working Directory**: Executes in current working directory

### Error Handling

Comprehensive error handling for command execution:

```python
def execute_with_error_handling(self) -> "Promise":
    try:
        # Execute command
        result = subprocess.check_output(
            self._command,
            shell=True,
            stderr=subprocess.STDOUT,
            executable=os.environ.get("SHELL"),
            timeout=300  # 5 minute timeout
        )
        self.resolve(result.decode())
        
    except subprocess.CalledProcessError as e:
        # Command returned non-zero exit code
        self.catch_error(e.output.decode())
        
    except subprocess.TimeoutExpired:
        # Command exceeded timeout
        self.catch_error("Command timed out")
        
    except FileNotFoundError:
        # Shell or command not found
        self.catch_error("Command or shell not found")
        
    except Exception as e:
        # Other execution errors
        self.catch_error(str(e))
    
    return self
```

## Advanced Usage Patterns

### Conditional Execution

```python
from celline.functions.bash import Bash

def conditional_execution(project):
    # Check if file exists before processing
    def file_exists(output: str):
        if "No such file" not in output:
            # File exists, proceed with processing
            process_function = Bash(
                cmd="python analyze_data.py data/input.csv",
                then=lambda out: print("Analysis complete"),
                catch=lambda err: print("Analysis failed")
            )
            project.call(process_function)
        else:
            print("Input file not found, skipping analysis")
    
    # Check file existence
    check_function = Bash(
        cmd="ls data/input.csv 2>&1",
        then=file_exists,
        catch=lambda err: print("File check failed")
    )
    
    project.call(check_function)
```

### Pipeline Orchestration

```python
from celline.functions.bash import Bash

class DataPipeline:
    def __init__(self, project):
        self.project = project
        self.step = 0
        self.steps = [
            "python preprocess.py",
            "Rscript normalize.R", 
            "python integrate.py",
            "Rscript visualize.R"
        ]
    
    def run_pipeline(self):
        self.execute_step(self.step)
    
    def execute_step(self, step_idx):
        if step_idx >= len(self.steps):
            print("Pipeline completed successfully")
            return
        
        def on_success(output):
            print(f"Step {step_idx + 1} completed")
            self.execute_step(step_idx + 1)
        
        def on_error(error):
            print(f"Pipeline failed at step {step_idx + 1}: {error}")
        
        bash_function = Bash(
            cmd=self.steps[step_idx],
            then=on_success,
            catch=on_error
        )
        
        self.project.call(bash_function)

# Usage
pipeline = DataPipeline(project)
pipeline.run_pipeline()
```

### Resource Monitoring

```python
from celline.functions.bash import Bash
import time

class ResourceMonitor:
    def __init__(self, project):
        self.project = project
        self.monitoring = True
    
    def start_monitoring(self):
        self.check_resources()
    
    def check_resources(self):
        def process_stats(output: str):
            lines = output.strip().split('\n')
            # Parse CPU, memory, disk usage
            for line in lines:
                if 'Cpu' in line:
                    print(f"CPU: {line}")
                elif 'MiB Mem' in line:
                    print(f"Memory: {line}")
            
            # Continue monitoring
            if self.monitoring:
                time.sleep(10)
                self.check_resources()
        
        bash_function = Bash(
            cmd="top -bn1 | head -5",
            then=process_stats,
            catch=lambda err: print("Monitoring failed")
        )
        
        self.project.call(bash_function)
    
    def stop_monitoring(self):
        self.monitoring = False

# Usage
monitor = ResourceMonitor(project)
monitor.start_monitoring()
```

## Security Considerations

### Command Injection Prevention

```python
import shlex

def safe_command_execution(user_input: str):
    # Sanitize user input
    safe_input = shlex.quote(user_input)
    
    # Use safe command construction
    command = f"ls {safe_input}"
    
    bash_function = Bash(
        cmd=command,
        then=lambda out: print("Safe execution"),
        catch=lambda err: print("Error occurred")
    )
    
    return bash_function
```

### Environment Isolation

```python
def isolated_execution(command: str, env_vars: dict):
    # Create isolated environment
    import os
    isolated_env = os.environ.copy()
    isolated_env.update(env_vars)
    
    # Execute with custom environment
    result = subprocess.run(
        command,
        shell=True,
        env=isolated_env,
        capture_output=True,
        text=True
    )
    
    return result
```

## Performance Optimization

### Parallel Execution

```python
import threading
from concurrent.futures import ThreadPoolExecutor

def parallel_bash_execution(commands: list):
    def execute_command(cmd: str):
        promise = Promise(cmd)
        return promise.execute()
    
    # Execute commands in parallel
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(execute_command, cmd) for cmd in commands]
        
        # Wait for completion
        for future in futures:
            result = future.result()
            print(f"Command completed: {result.is_done()}")
```

### Memory-Efficient Streaming

```python
def stream_output(command: str):
    """Stream command output for large outputs."""
    process = subprocess.Popen(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )
    
    # Stream output line by line
    for line in iter(process.stdout.readline, ''):
        yield line.strip()
    
    process.stdout.close()
    return_code = process.wait()
    
    if return_code != 0:
        error = process.stderr.read()
        raise subprocess.CalledProcessError(return_code, command, error)
```

## Methods

### Bash Methods

#### `call(project: Project) -> Project`

Main execution method that runs the shell command.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Creates Promise with configured command
2. Attaches success and error callbacks
3. Executes command through Promise.execute()
4. Returns project instance

### Promise Methods

#### `execute() -> Promise`

Executes the shell command and triggers callbacks.

**Returns**: Self for method chaining

#### `then(callback: Callable[[str], None]) -> Promise`

Registers success callback function.

**Parameters**:
- `callback`: Function to call on successful execution

**Returns**: Self for method chaining

#### `error(callback: Callable[[str], None]) -> Promise`

Registers error callback function.

**Parameters**:
- `callback`: Function to call on execution error

**Returns**: Self for method chaining

#### `is_done() -> bool`

Checks if command execution is complete.

**Returns**: True if execution finished

## Common Use Cases

### System Diagnostics

```python
def system_diagnostics(project):
    """Run comprehensive system diagnostics."""
    
    diagnostics = [
        ("Disk Space", "df -h"),
        ("Memory Usage", "free -h"),
        ("CPU Info", "lscpu | head -10"),
        ("Process Count", "ps aux | wc -l"),
        ("Network Status", "ping -c 1 google.com")
    ]
    
    for name, cmd in diagnostics:
        bash_function = Bash(
            cmd=cmd,
            then=lambda out, n=name: print(f"{n}: {out.strip()}"),
            catch=lambda err, n=name: print(f"{n}: Failed")
        )
        project.call(bash_function)
```

### Data Validation

```python
def validate_data_files(project):
    """Validate data file integrity."""
    
    def check_file_integrity(output: str):
        files = output.strip().split('\n')
        print(f"Found {len(files)} data files")
        
        for file in files:
            if file.strip():
                # Check each file
                check_cmd = f"file {file.strip()}"
                file_check = Bash(
                    cmd=check_cmd,
                    then=lambda out: print(f"✓ {file}: Valid"),
                    catch=lambda err: print(f"✗ {file}: Invalid")
                )
                project.call(file_check)
    
    # Find all data files
    find_function = Bash(
        cmd="find data/ -name '*.h5' -o -name '*.csv' -o -name '*.tsv'",
        then=check_file_integrity,
        catch=lambda err: print("File search failed")
    )
    
    project.call(find_function)
```

### Custom Tool Integration

```python
def integrate_custom_tool(project, tool_path: str, input_file: str):
    """Integrate external analysis tool."""
    
    def tool_success(output: str):
        print("Custom tool completed successfully")
        # Parse tool output
        results = parse_tool_output(output)
        save_results(results)
    
    def tool_error(error: str):
        print(f"Custom tool failed: {error}")
        # Handle tool-specific errors
        if "memory" in error.lower():
            print("Try reducing input size or increasing memory")
        elif "permission" in error.lower():
            print("Check file permissions")
    
    # Execute custom tool
    bash_function = Bash(
        cmd=f"{tool_path} --input {input_file} --output results/",
        then=tool_success,
        catch=tool_error
    )
    
    project.call(bash_function)
```

## Related Functions

- [Job](job) - Monitor long-running bash commands
- [Interactive](interactive) - Execute commands through web interface
- [Info](info) - System information gathering
- All analysis functions - May use bash for subprocess execution

## Troubleshooting

### Common Issues

1. **Command Not Found**: Executable not in PATH
2. **Permission Denied**: Insufficient permissions
3. **Timeout**: Command takes too long to execute
4. **Shell Differences**: Commands behave differently across shells
5. **Environment Variables**: Missing required environment variables

### Debug Strategies

```python
def debug_bash_execution(command: str):
    """Debug bash command execution."""
    
    def debug_success(output: str):
        print(f"Command: {command}")
        print(f"Output length: {len(output)}")
        print(f"Output preview: {output[:200]}...")
    
    def debug_error(error: str):
        print(f"Command: {command}")
        print(f"Error: {error}")
        print(f"Shell: {os.environ.get('SHELL', 'unknown')}")
        print(f"PATH: {os.environ.get('PATH', 'unknown')}")
    
    bash_function = Bash(
        cmd=command,
        then=debug_success,
        catch=debug_error
    )
    
    return bash_function
```

### Error Recovery

```python
def resilient_execution(commands: list):
    """Execute commands with automatic retry."""
    
    def retry_on_failure(cmd: str, attempts: int = 3):
        for attempt in range(attempts):
            try:
                promise = Promise(cmd)
                result = promise.execute()
                if promise.is_done():
                    return result
            except Exception as e:
                if attempt == attempts - 1:
                    print(f"Command failed after {attempts} attempts: {cmd}")
                    raise e
                else:
                    print(f"Attempt {attempt + 1} failed, retrying...")
                    time.sleep(2 ** attempt)  # Exponential backoff
    
    for cmd in commands:
        retry_on_failure(cmd)
```