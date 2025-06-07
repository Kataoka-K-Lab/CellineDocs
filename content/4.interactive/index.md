---
title: Interactive Mode
navigation: false      # ← 見出し専用にしてリンクは非表示
---

This document provides a detailed explanation of Celline's web-based interactive interface.
It's a feature that allows intuitive analysis and result visualization directly in the browser.

## 🌐 Overview

Interactive mode is built with the following technologies:

- **Backend**: FastAPI (Python)
- **Frontend**: Vue.js + Nuxt.js
- **Communication**: RESTful API + WebSocket
- **Visualization**: Plotly.js + D3.js

## 🚀 Starting Interactive Mode

### Basic Startup

```bash
# Start interactive mode
celline interactive

# or
celline run interactive
```

After startup, access http://localhost:8080 in your browser.

### Startup with Custom Settings

```bash
# Change port number
celline interactive --port 8090

# Start API only (for development)
celline api --port 8000
```

### Starting from Python API

```python
from celline.functions.interactive import Interactive

# Start with default settings
interactive = Interactive()
project.call(interactive)

# Start with custom settings
interactive = Interactive(port=8090, debug=True)
project.call(interactive)
```

## 🖥️ User Interface

### Main Screen Layout

```
┌─────────────────────────────────────────────────────┐
│ Celline Interactive Dashboard                       │
├─────────────────────────────────────────────────────┤
│ [Project Info] [Sample Management] [Analysis] [Help]│
├─────────────────────────────────────────────────────┤
│ ┌─Sidebar─────┐ ┌─Main Content─────────────────────┐ │
│ │ • Projects  │ │                                │ │
│ │ • Samples   │ │     Visualization Area         │ │
│ │ • Functions │ │                                │ │
│ │ • Results   │ │                                │ │
│ │ • Settings  │ │                                │ │
│ │ └─────────────┘ └──────────────────────────────────┘ │
├─────────────────────────────────────────────────────┤
│ Status Bar: Ready | Samples: 5 | Memory: 2.1GB     │
└─────────────────────────────────────────────────────┘
```

### Project Overview

#### Project Information Panel

- **Project Name**: Current project name
- **Sample Count**: Number of added samples
- **Processing Status**: Progress of each stage
- **Resource Usage**: Memory and storage usage

```vue
<template>
  <div class="project-overview">
    <h2>{{ projectInfo.name }}</h2>
    <div class="stats-grid">
      <div class="stat-card">
        <h3>Sample Count</h3>
        <p class="stat-value">{{ projectInfo.samples.length }}</p>
      </div>
      <div class="stat-card">
        <h3>Processed</h3>
        <p class="stat-value">{{ processedCount }}</p>
      </div>
      <div class="stat-card">
        <h3>Analyzed</h3>
        <p class="stat-value">{{ analyzedCount }}</p>
      </div>
    </div>
  </div>
</template>
```

## 📊 Sample Management

### Sample List

You can add, remove, and check the status of samples.

#### Sample States

| State | Description | Actions |
|-------|-------------|---------|
| `pending` | Added, not downloaded | Download |
| `downloading` | Download in progress | Show progress |
| `downloaded` | Download completed | Count |
| `counting` | Count processing | Show progress |
| `counted` | Count completed | Preprocess |
| `processed` | Preprocessing completed | Analysis |

#### Sample Addition

```typescript
// TypeScript (Frontend)
interface AddSampleRequest {
  sample_ids: string[];
  titles?: string[];
}

async function addSamples(sampleIds: string[]) {
  const response = await fetch('/api/samples/add', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ sample_ids: sampleIds })
  });
  
  const result = await response.json();
  return result.job_id;
}
```

### Batch Operations

Batch operations on multiple samples are possible.

```javascript
// Batch download for selected samples
async function batchDownload(selectedSamples) {
  const jobs = [];
  
  for (const sample of selectedSamples) {
    const jobId = await startDownload(sample.id);
    jobs.push(jobId);
  }
  
  // Monitor jobs
  monitorJobs(jobs);
}
```

## 🔧 Function Execution

### Available Functions

The following functions can be executed from the interface:

#### Data Management
- **Add**: Sample addition
- **Download**: Data download
- **SyncDB**: Database synchronization

#### Data Processing
- **Count**: Cell Ranger count
- **Preprocess**: Preprocessing & QC

#### Analysis
- **PredictCelltype**: Cell type prediction
- **Reduce**: Dimensionality reduction
- **Integrate**: Data integration

### Function Execution Interface

```vue
<template>
  <div class="function-panel">
    <h3>{{ functionInfo.name }}</h3>
    <p>{{ functionInfo.description }}</p>
    
    <!-- Parameter input -->
    <div v-for="param in functionInfo.parameters" :key="param.name">
      <label>{{ param.label }}</label>
      <input 
        v-model="parameters[param.name]" 
        :type="param.type"
        :placeholder="param.default"
      />
    </div>
    
    <!-- Execute button -->
    <button @click="executeFunction" :disabled="isRunning">
      {{ isRunning ? 'Running...' : 'Execute' }}
    </button>
  </div>
</template>
```

### Job Monitoring

```typescript
interface JobStatus {
  job_id: string;
  status: 'pending' | 'running' | 'completed' | 'failed';
  progress: number;
  message: string;
  created_at: string;
}

class JobMonitor {
  private jobs: Map<string, JobStatus> = new Map();
  
  async startMonitoring(jobId: string) {
    const interval = setInterval(async () => {
      const status = await this.checkJobStatus(jobId);
      this.jobs.set(jobId, status);
      
      if (status.status === 'completed' || status.status === 'failed') {
        clearInterval(interval);
        this.onJobComplete(status);
      }
    }, 2000);
  }
  
  async checkJobStatus(jobId: string): Promise<JobStatus> {
    const response = await fetch(`/api/jobs/${jobId}`);
    return response.json();
  }
}
```

## 📈 Data Visualization

### QC Metrics

Visualizes quality control metrics after preprocessing.

```javascript
// QC plots using Plotly.js
function plotQCMetrics(cellData) {
  // Gene count distribution
  const geneCountTrace = {
    x: cellData.map(cell => cell.n_genes_by_counts),
    type: 'histogram',
    name: 'Gene Count Distribution'
  };
  
  // Mitochondrial gene percentage
  const mitoTrace = {
    x: cellData.map(cell => cell.pct_counts_mt),
    type: 'histogram',
    name: 'Mitochondrial Gene %'
  };
  
  Plotly.newPlot('qc-plots', [geneCountTrace, mitoTrace]);
}
```

### UMAP/t-SNE Visualization

```javascript
function plotUMAP(umapData, colorBy = 'celltype') {
  const trace = {
    x: umapData.map(point => point.umap_1),
    y: umapData.map(point => point.umap_2),
    mode: 'markers',
    type: 'scatter',
    marker: {
      color: umapData.map(point => point[colorBy]),
      colorscale: 'Viridis'
    },
    text: umapData.map(point => `Cell: ${point.cell_id}<br>Type: ${point.celltype}`)
  };
  
  const layout = {
    title: 'UMAP Visualization',
    xaxis: { title: 'UMAP 1' },
    yaxis: { title: 'UMAP 2' }
  };
  
  Plotly.newPlot('umap-plot', [trace], layout);
}
```

### Interactive Gene Expression

```vue
<template>
  <div class="gene-expression-panel">
    <div class="gene-selector">
      <input 
        v-model="selectedGene" 
        @input="updateExpression"
        placeholder="Enter gene name (e.g., CD4)"
        list="gene-list"
      />
      <datalist id="gene-list">
        <option v-for="gene in availableGenes" :value="gene" :key="gene"/>
      </datalist>
    </div>
    
    <div id="expression-plot" class="plot-container"></div>
  </div>
</template>

<script>
export default {
  data() {
    return {
      selectedGene: '',
      availableGenes: [],
      expressionData: []
    };
  },
  
  methods: {
    async updateExpression() {
      if (this.selectedGene) {
        const data = await this.fetchGeneExpression(this.selectedGene);
        this.plotGeneExpression(data);
      }
    },
    
    async fetchGeneExpression(gene) {
      const response = await fetch(`/api/expression/${gene}`);
      return response.json();
    },
    
    plotGeneExpression(data) {
      // UMAP plot with expression levels represented by color
      const trace = {
        x: data.map(point => point.umap_1),
        y: data.map(point => point.umap_2),
        mode: 'markers',
        type: 'scatter',
        marker: {
          color: data.map(point => point.expression),
          colorscale: 'Reds',
          colorbar: { title: 'Expression Level' }
        }
      };
      
      Plotly.newPlot('expression-plot', [trace]);
    }
  }
};
</script>
```

## 🔄 Real-time Updates

### WebSocket Communication

```javascript
class CellineWebSocket {
  constructor(url) {
    this.ws = new WebSocket(url);
    this.setupEventHandlers();
  }
  
  setupEventHandlers() {
    this.ws.onmessage = (event) => {
      const data = JSON.parse(event.data);
      this.handleMessage(data);
    };
    
    this.ws.onopen = () => {
      console.log('WebSocket connected');
    };
    
    this.ws.onclose = () => {
      console.log('WebSocket disconnected');
      // Reconnection logic
      setTimeout(() => this.reconnect(), 5000);
    };
  }
  
  handleMessage(data) {
    switch (data.type) {
      case 'job_update':
        this.updateJobStatus(data.job_id, data.status);
        break;
      case 'sample_update':
        this.updateSampleStatus(data.sample_id, data.status);
        break;
      case 'system_notification':
        this.showNotification(data.message);
        break;
    }
  }
}
```

### Progress Bar

```vue
<template>
  <div class="progress-container">
    <div class="progress-header">
      <span>{{ job.name }}</span>
      <span>{{ Math.round(job.progress) }}%</span>
    </div>
    <div class="progress-bar">
      <div 
        class="progress-fill" 
        :style="{ width: job.progress + '%' }"
        :class="progressClass"
      ></div>
    </div>
    <div class="progress-message">{{ job.message }}</div>
  </div>
</template>

<style scoped>
.progress-container {
  margin: 10px 0;
  padding: 10px;
  border: 1px solid #ddd;
  border-radius: 4px;
}

.progress-bar {
  width: 100%;
  height: 20px;
  background-color: #f0f0f0;
  border-radius: 10px;
  overflow: hidden;
}

.progress-fill {
  height: 100%;
  transition: width 0.3s ease;
}

.progress-fill.running {
  background-color: #007bff;
}

.progress-fill.completed {
  background-color: #28a745;
}

.progress-fill.failed {
  background-color: #dc3545;
}
</style>
```

## 📊 Result Downloads

### Export Features

```typescript
interface ExportOptions {
  format: 'csv' | 'excel' | 'pdf' | 'png';
  data_type: 'qc_metrics' | 'cell_info' | 'expression' | 'plots';
  samples?: string[];
}

async function exportResults(options: ExportOptions) {
  const response = await fetch('/api/export', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(options)
  });
  
  if (response.ok) {
    const blob = await response.blob();
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `celline_results.${options.format}`;
    a.click();
  }
}
```

### Report Generation

```vue
<template>
  <div class="report-generator">
    <h3>Analysis Report Generation</h3>
    
    <div class="report-options">
      <label>
        <input type="checkbox" v-model="includeQC"> QC Metrics
      </label>
      <label>
        <input type="checkbox" v-model="includePlots"> Visualization Plots
      </label>
      <label>
        <input type="checkbox" v-model="includeCellTypes"> Cell Type Information
      </label>
    </div>
    
    <button @click="generateReport" class="btn-primary">
      Generate Report
    </button>
  </div>
</template>

<script>
export default {
  data() {
    return {
      includeQC: true,
      includePlots: true,
      includeCellTypes: true
    };
  },
  
  methods: {
    async generateReport() {
      const options = {
        include_qc: this.includeQC,
        include_plots: this.includePlots,
        include_celltypes: this.includeCellTypes,
        format: 'pdf'
      };
      
      const response = await fetch('/api/generate-report', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(options)
      });
      
      if (response.ok) {
        const blob = await response.blob();
        this.downloadFile(blob, 'analysis_report.pdf');
      }
    }
  }
};
</script>
```

## ⚙️ Settings and Customization

### User Settings

```vue
<template>
  <div class="settings-panel">
    <h3>Settings</h3>
    
    <div class="setting-group">
      <h4>Execution Environment</h4>
      <select v-model="settings.execution_system">
        <option value="multithreading">Multithreading</option>
        <option value="PBS">PBS Cluster</option>
      </select>
      
      <label>Thread Count:</label>
      <input 
        type="number" 
        v-model="settings.nthread" 
        min="1" 
        max="64"
      />
    </div>
    
    <div class="setting-group">
      <h4>Display Settings</h4>
      <label>
        <input type="checkbox" v-model="settings.dark_mode"> Dark Mode
      </label>
      <label>
        <input type="checkbox" v-model="settings.auto_refresh"> Auto Refresh
      </label>
    </div>
    
    <button @click="saveSettings">Save Settings</button>
  </div>
</template>
```

## 🚨 Error Handling

### Error Display

```vue
<template>
  <div v-if="error" class="error-banner">
    <div class="error-content">
      <i class="icon-error"></i>
      <div class="error-message">
        <h4>{{ error.title }}</h4>
        <p>{{ error.message }}</p>
      </div>
      <button @click="dismissError" class="btn-close">×</button>
    </div>
    
    <div v-if="error.details" class="error-details">
      <details>
        <summary>Details</summary>
        <pre>{{ error.details }}</pre>
      </details>
    </div>
  </div>
</template>
```

### Automatic Recovery

```javascript
class ErrorHandler {
  constructor() {
    this.retryAttempts = 0;
    this.maxRetries = 3;
  }
  
  async handleAPIError(error, originalRequest) {
    if (error.status === 500 && this.retryAttempts < this.maxRetries) {
      this.retryAttempts++;
      await this.delay(1000 * this.retryAttempts);
      return this.retryRequest(originalRequest);
    }
    
    this.showErrorMessage(error);
  }
  
  delay(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
  }
}
```

---

> **Success**: For detailed usage of interactive mode, we recommend experiencing it by actually launching it in your browser. You can start right now with `celline interactive`!