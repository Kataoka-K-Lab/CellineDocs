<template>
  <div class="mermaid-container">
    <div ref="mermaidRef" class="mermaid-diagram"></div>
  </div>
</template>

<script setup lang="ts">
import { ref, onMounted, nextTick } from 'vue'

interface Props {
  code: string
  language?: string
}

const props = withDefaults(defineProps<Props>(), {
  language: 'mermaid'
})

const mermaidRef = ref<HTMLElement>()

onMounted(async () => {
  try {
    // Dynamically import mermaid to avoid SSR issues
    const mermaid = await import('mermaid')
    
    // Initialize mermaid with configuration
    mermaid.default.initialize({
      startOnLoad: false,
      theme: 'default',
      securityLevel: 'loose',
      themeVariables: {
        primaryColor: '#3b82f6',
        primaryTextColor: '#1f2937',
        primaryBorderColor: '#3b82f6',
        lineColor: '#6b7280',
        sectionBkgColor: '#f9fafb',
        altSectionBkgColor: '#ffffff',
        gridColor: '#e5e7eb',
        secondaryColor: '#e5e7eb',
        tertiaryColor: '#f3f4f6'
      }
    })

    await nextTick()

    if (mermaidRef.value) {
      try {
        // Generate unique ID for the diagram
        const id = `mermaid-${Math.random().toString(36).substr(2, 9)}`
        
        // Render the mermaid diagram
        const { svg } = await mermaid.default.render(id, props.code)
        
        // Insert the SVG into the DOM
        mermaidRef.value.innerHTML = svg
        
        // Add responsive styling
        const svgElement = mermaidRef.value.querySelector('svg')
        if (svgElement) {
          svgElement.style.maxWidth = '100%'
          svgElement.style.height = 'auto'
        }
      } catch (error) {
        console.error('Error rendering Mermaid diagram:', error)
        mermaidRef.value.innerHTML = `
          <div class="mermaid-error p-4 border border-red-300 rounded bg-red-50 text-red-700">
            <strong>Mermaid Diagram Error:</strong>
            <pre class="mt-2 text-sm">${error}</pre>
            <details class="mt-2">
              <summary class="cursor-pointer">Show diagram code</summary>
              <pre class="mt-2 text-xs bg-gray-100 p-2 rounded">${props.code}</pre>
            </details>
          </div>
        `
      }
    }
  } catch (error) {
    console.error('Error loading Mermaid:', error)
    if (mermaidRef.value) {
      mermaidRef.value.innerHTML = `
        <div class="mermaid-error p-4 border border-red-300 rounded bg-red-50 text-red-700">
          <strong>Mermaid not available</strong>
          <pre class="mt-2 text-sm">${props.code}</pre>
        </div>
      `
    }
  }
})
</script>

<style scoped>
.mermaid-container {
  @apply my-6 overflow-x-auto;
}

.mermaid-diagram {
  @apply flex justify-center items-center min-h-[200px];
}

/* Dark mode support */
.dark .mermaid-diagram :deep(svg) {
  filter: invert(1) hue-rotate(180deg);
}

.dark .mermaid-diagram :deep(.node rect),
.dark .mermaid-diagram :deep(.node circle),
.dark .mermaid-diagram :deep(.node ellipse),
.dark .mermaid-diagram :deep(.node polygon) {
  filter: invert(1) hue-rotate(180deg);
}

/* Responsive design */
@media (max-width: 768px) {
  .mermaid-diagram :deep(svg) {
    transform: scale(0.8);
    transform-origin: top left;
  }
}
</style>