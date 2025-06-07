# Interactive Function

Launch a Vue.js web interface for visual project management and analysis workflows.

## Overview

The `Interactive` function starts a modern Vue.js web application that provides a graphical user interface for Celline. It offers visual project management, sample tracking, workflow execution, and real-time analysis monitoring through an intuitive browser-based interface.

## Class Information

- **Module**: `celline.functions.interactive`
- **Class**: `Interactive`
- **Base Class**: `CellineFunction`

## Parameters

### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `host` | `str` | No | Host to bind server (default: "localhost") |
| `port` | `int` | No | Port number (default: 3000) |
| `auto_open` | `bool` | No | Auto-open browser (default: True) |

## Usage Examples

### Python API

#### Basic Web Interface

```python
from celline import Project
from celline.functions.interactive import Interactive

# Create project
project = Project("./my-project")

# Launch web interface with defaults
interactive_function = Interactive()

# Execute function (starts web server)
result = project.call(interactive_function)
```

#### Custom Server Configuration

```python
from celline import Project
from celline.functions.interactive import Interactive

# Create project
project = Project("./my-project")

# Launch with custom configuration
interactive_function = Interactive(
    host="0.0.0.0",     # Accept external connections
    port=8080,          # Custom port
    auto_open=False     # Don't auto-open browser
)

# Execute function
result = project.call(interactive_function)
```

#### Development Mode

```python
from celline import Project
from celline.functions.interactive import Interactive

# Create project
project = Project("./my-project")

# Development configuration
interactive_function = Interactive(
    host="localhost",
    port=3000,
    auto_open=True
)

# Launch development server
result = project.call(interactive_function)
```

#### Production Deployment

```python
from celline import Project
from celline.functions.interactive import Interactive

# Create project
project = Project("./my-project")

# Production configuration
interactive_function = Interactive(
    host="0.0.0.0",     # Allow external access
    port=80,            # Standard web port
    auto_open=False     # No auto-open for servers
)

# Launch production server
result = project.call(interactive_function)
```

### CLI Usage

#### Basic Usage

```bash
# Launch web interface with defaults
celline run interactive

# Access at http://localhost:3000
```

#### Custom Configuration

```bash
# Custom port
celline run interactive --port 8080

# External access
celline run interactive --host 0.0.0.0 --port 3000

# No auto-open
celline run interactive --no-open

# Combined options
celline run interactive --host 0.0.0.0 --port 8080 --no-open
```

### CLI Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--host` | `str` | localhost | Host to bind server |
| `--port`, `-p` | `int` | 3000 | Port number |
| `--no-open` | `flag` | False | Don't auto-open browser |

## Implementation Details

### Web Technology Stack

The interface uses modern web technologies:

| Component | Technology | Purpose |
|-----------|------------|---------|
| **Frontend** | Vue.js 3 | Reactive user interface |
| **Build System** | Vue CLI | Development and production builds |
| **Styling** | CSS3/SCSS | Modern responsive design |
| **State Management** | Vuex/Pinia | Application state |
| **Routing** | Vue Router | Single-page navigation |
| **Components** | Vue Components | Modular UI elements |

### Frontend Development Pipeline

The function automatically manages frontend setup:

1. **Dependency Detection**: Checks Node.js and npm availability
2. **Project Creation**: Creates Vue.js project if needed
3. **Package Installation**: Installs required dependencies
4. **TypeScript Configuration**: Sets up TypeScript support
5. **Development Server**: Launches Vue development server
6. **Browser Integration**: Opens interface in default browser

### Automatic Setup Process

```javascript
// Vue.js project structure created automatically
frontend/
├── package.json              # Project dependencies
├── tsconfig.json            # TypeScript configuration
├── vue.config.js            # Vue CLI configuration
├── src/
│   ├── App.vue              # Main application component
│   ├── main.ts              # Application entry point
│   ├── components/          # Reusable UI components
│   │   ├── SampleManager.vue
│   │   ├── WorkflowRunner.vue
│   │   └── ProgressMonitor.vue
│   ├── views/               # Page components
│   │   ├── Dashboard.vue
│   │   ├── Samples.vue
│   │   └── Analysis.vue
│   └── router/              # Navigation configuration
│       └── index.ts
└── public/
    ├── index.html           # HTML template
    └── favicon.ico
```

### Vue.js Application Features

#### Dashboard Interface

```vue
<template>
  <div class="dashboard">
    <header class="dashboard-header">
      <h1>Celline Project Dashboard</h1>
      <project-info :project="currentProject" />
    </header>
    
    <main class="dashboard-content">
      <sample-overview :samples="samples" />
      <workflow-status :workflows="activeWorkflows" />
      <storage-monitor :usage="storageUsage" />
    </main>
  </div>
</template>

<script setup lang="ts">
import { ref, onMounted } from 'vue'
import ProjectInfo from '@/components/ProjectInfo.vue'
import SampleOverview from '@/components/SampleOverview.vue'

const currentProject = ref(null)
const samples = ref([])
const activeWorkflows = ref([])

onMounted(async () => {
  await loadProjectData()
})
</script>
```

#### Sample Management

```vue
<template>
  <div class="sample-manager">
    <sample-upload @upload="handleSampleUpload" />
    <sample-table 
      :samples="samples"
      @edit="editSample"
      @delete="deleteSample"
    />
    <sample-status-chart :data="statusData" />
  </div>
</template>
```

#### Workflow Execution

```vue
<template>
  <div class="workflow-runner">
    <workflow-selector 
      :available="availableWorkflows"
      @select="selectWorkflow"
    />
    <parameter-editor 
      :workflow="selectedWorkflow"
      :parameters="workflowParameters"
    />
    <execution-controls
      @run="executeWorkflow"
      @stop="stopWorkflow"
    />
    <progress-monitor :execution="currentExecution" />
  </div>
</template>
```

## User Interface Features

### Project Management

- **Project Overview**: Visual project status and metrics
- **Sample Registry**: Interactive sample management interface
- **File Browser**: Navigate project directory structure
- **Configuration Editor**: Visual editing of project settings

### Analysis Workflows

- **Workflow Builder**: Drag-and-drop workflow creation
- **Parameter Configuration**: Form-based parameter editing
- **Execution Monitoring**: Real-time progress tracking
- **Result Visualization**: Interactive plot display

### Data Visualization

- **Quality Control Plots**: Sample QC metrics visualization
- **Cell Type Plots**: Interactive cell type distributions
- **Integration Views**: Batch correction and integration results
- **Custom Plots**: User-defined visualization panels

### System Monitoring

- **Resource Usage**: CPU, memory, and storage monitoring
- **Job Status**: Running job tracking and management
- **Log Viewer**: Real-time log display and filtering
- **Error Reporting**: User-friendly error messages

## Server Architecture

### Development Server

The function uses Vue CLI development server:

```javascript
// vue.config.js
module.exports = {
  devServer: {
    host: process.env.HOST || 'localhost',
    port: process.env.PORT || 3000,
    open: true,
    hot: true,
    overlay: {
      warnings: false,
      errors: true
    }
  },
  configureWebpack: {
    resolve: {
      alias: {
        '@': path.resolve(__dirname, 'src')
      }
    }
  }
}
```

### API Integration

Frontend communicates with Celline backend:

```typescript
// API service for backend communication
class CellineAPI {
  private baseURL: string

  constructor(baseURL: string = 'http://localhost:8000') {
    this.baseURL = baseURL
  }

  async getProjectInfo(): Promise<ProjectInfo> {
    const response = await fetch(`${this.baseURL}/api/project`)
    return response.json()
  }

  async getSamples(): Promise<Sample[]> {
    const response = await fetch(`${this.baseURL}/api/samples`)
    return response.json()
  }

  async executeWorkflow(workflow: Workflow): Promise<ExecutionResult> {
    const response = await fetch(`${this.baseURL}/api/execute`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(workflow)
    })
    return response.json()
  }
}
```

## Methods

### `register() -> str`

Returns the function identifier for registration.

### `call(project: Project) -> Project`

Main execution method that launches the web interface.

**Parameters**:
- `project`: The Celline project instance

**Returns**: Updated project instance

**Process**:
1. Validates Node.js and npm installation
2. Sets up Vue.js project structure if needed
3. Installs frontend dependencies
4. Configures development environment
5. Starts Vue development server
6. Opens browser automatically (if enabled)
7. Monitors server until shutdown

### `add_cli_args(parser: ArgumentParser) -> None`

Adds CLI-specific arguments to the argument parser.

### `cli(project: Project, args: Namespace) -> Project`

CLI entry point that processes command-line arguments and executes the function.

### `get_description() -> str`

Returns function description for help documentation.

### `get_usage_examples() -> list[str]`

Returns usage examples for CLI help.

## Error Handling

### Common Issues

1. **Node.js Missing**: Node.js not installed or not in PATH
2. **npm Issues**: Package installation failures
3. **Port Conflicts**: Requested port already in use
4. **Permission Errors**: Cannot create files or bind ports
5. **Browser Access**: Browser cannot connect to server

### Dependency Resolution

```python
def resolve_frontend_dependencies():
    # Check Node.js
    if not _check_node_installed():
        console.print("[red]Node.js required for interactive mode[/red]")
        console.print("Install from: https://nodejs.org/")
        return False
    
    # Check npm
    if not _check_npm_installed():
        console.print("[red]npm required for package management[/red]")
        return False
    
    return True
```

### Error Recovery

```python
def handle_frontend_errors(error_type):
    if error_type == "dependency_conflict":
        # Clean npm cache and retry
        _clean_npm_cache()
        _fix_package_lock(frontend_path)
        return _try_install_with_legacy_deps(frontend_path)
    
    elif error_type == "port_conflict":
        # Try alternative ports
        for port in range(3001, 3010):
            if _port_available(port):
                return _start_server_on_port(port)
    
    elif error_type == "permission_denied":
        console.print("[yellow]Try running with different permissions[/yellow]")
        return False
```

## Performance Considerations

### Development Performance

- **Hot Module Replacement**: Fast development iteration
- **Incremental Compilation**: TypeScript incremental builds
- **Asset Optimization**: Automatic code splitting
- **Memory Management**: Efficient Vue.js reactivity

### Production Optimization

```javascript
// Production build configuration
module.exports = {
  productionSourceMap: false,
  configureWebpack: config => {
    if (process.env.NODE_ENV === 'production') {
      config.optimization.splitChunks = {
        chunks: 'all',
        cacheGroups: {
          vendor: {
            name: 'vendor',
            test: /[\\/]node_modules[\\/]/,
            priority: 10,
            chunks: 'initial'
          }
        }
      }
    }
  }
}
```

### Resource Usage

| Component | Memory Usage | CPU Usage |
|-----------|--------------|-----------|
| Vue Development Server | 100-200 MB | Low |
| Node.js Runtime | 50-100 MB | Low |
| Browser Instance | 200-500 MB | Medium |
| WebSocket Connections | 10-50 MB | Low |

## Security Considerations

### Development Security

- **CORS Configuration**: Proper cross-origin settings
- **Input Validation**: Client-side input sanitization
- **Authentication**: Session management (if implemented)
- **HTTPS Support**: SSL/TLS configuration options

### Production Deployment

```typescript
// Security middleware for production
app.use(helmet())  // Security headers
app.use(cors({
  origin: process.env.ALLOWED_ORIGINS?.split(','),
  credentials: true
}))
app.use(rateLimit({
  windowMs: 15 * 60 * 1000, // 15 minutes
  max: 100 // limit each IP to 100 requests per windowMs
}))
```

## Browser Compatibility

### Supported Browsers

| Browser | Version | Features |
|---------|---------|----------|
| **Chrome** | 90+ | Full support |
| **Firefox** | 88+ | Full support |
| **Safari** | 14+ | Full support |
| **Edge** | 90+ | Full support |

### Feature Detection

```typescript
// Progressive enhancement
if ('serviceWorker' in navigator) {
  // Enable offline functionality
}

if ('WebSocket' in window) {
  // Enable real-time updates
}

if (CSS.supports('display', 'grid')) {
  // Use modern layout
}
```

## Customization

### Theme Configuration

```scss
// Custom theme variables
$primary-color: #2c3e50;
$secondary-color: #3498db;
$success-color: #27ae60;
$warning-color: #f39c12;
$error-color: #e74c3c;

// Dark mode support
@media (prefers-color-scheme: dark) {
  :root {
    --bg-color: #1a1a1a;
    --text-color: #ffffff;
  }
}
```

### Component Extension

```vue
<!-- Custom workflow component -->
<template>
  <div class="custom-workflow">
    <workflow-base :config="workflowConfig">
      <template #parameters>
        <custom-parameter-editor />
      </template>
      <template #visualization>
        <custom-results-viewer />
      </template>
    </workflow-base>
  </div>
</template>
```

## Related Functions

- [Initialize](init) - Set up project before launching interface
- [Info](info) - Display project information (used by interface)
- [Job](job) - Monitor jobs (integrated in interface)
- All analysis functions - Accessible through web interface

## Troubleshooting

### Common Issues

1. **Server Won't Start**: Check port availability and permissions
2. **Dependencies Failed**: Clear npm cache and reinstall
3. **Browser Can't Connect**: Verify firewall and host settings
4. **TypeScript Errors**: Check tsconfig.json configuration
5. **Build Failures**: Verify Node.js version compatibility

### Debug Mode

Enable detailed debugging:

```bash
# Enable Vue.js debug mode
NODE_ENV=development celline run interactive

# Enable detailed npm logging
npm config set loglevel verbose
```

### Manual Frontend Setup

For troubleshooting, set up frontend manually:

```bash
# Create Vue project manually
vue create celline-frontend
cd celline-frontend

# Install dependencies
npm install

# Start development server
npm run serve
```

### Network Configuration

For network access issues:

```bash
# Check port availability
netstat -an | grep :3000

# Test server binding
curl http://localhost:3000

# Check firewall rules
sudo ufw status
```