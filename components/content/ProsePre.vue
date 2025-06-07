<script setup lang="ts">
import { ref, computed, onMounted } from 'vue'

/* Nuxt Content が差し込む <code class="language-js"> … */
const codeRef = ref<HTMLElement>()
const copied = ref(false)
const mermaidRef = ref<HTMLElement>()

// Props that might be passed by Nuxt Content
const props = defineProps<{
  code?: string
  language?: string
  filename?: string
  highlights?: number[]
  meta?: string
}>()

// Detect if this is a mermaid diagram
const isMermaid = computed(() => {
  const slots = useSlots()
  if (slots.default) {
    const content = slots.default()
    if (content && content[0] && content[0].type === 'code') {
      const codeElement = content[0]
      const className = codeElement.props?.class || ''
      return className.includes('language-mermaid')
    }
  }
  return props.language === 'mermaid' || props.meta === 'mermaid'
})

// Extract code content
const codeContent = computed(() => {
  if (props.code) return props.code
  const slots = useSlots()
  if (slots.default) {
    const content = slots.default()
    if (content && content[0] && typeof content[0].children === 'string') {
      return content[0].children
    }
    if (content && content[0] && content[0].children && Array.isArray(content[0].children)) {
      return content[0].children.map(child => {
        if (typeof child === 'string') return child
        if (child.children) return child.children
        return ''
      }).join('')
    }
  }
  return ''
})

async function copy() {
  if (!codeRef.value) return
  await navigator.clipboard.writeText(codeRef.value.innerText)
  copied.value = true
  setTimeout(() => { copied.value = false }, 1500)
}

async function renderMermaid() {
  if (!isMermaid.value || !mermaidRef.value) return
  
  try {
    const mermaid = await import('mermaid')
    
    // Initialize mermaid if not already done
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

    // Generate unique ID
    const id = `mermaid-${Math.random().toString(36).substr(2, 9)}`
    
    // Render the diagram
    const { svg } = await mermaid.default.render(id, codeContent.value)
    
    // Insert the SVG
    mermaidRef.value.innerHTML = svg
    
    // Make it responsive
    const svgElement = mermaidRef.value.querySelector('svg')
    if (svgElement) {
      svgElement.style.maxWidth = '100%'
      svgElement.style.height = 'auto'
    }
  } catch (error) {
    console.error('Error rendering Mermaid diagram:', error)
    if (mermaidRef.value) {
      mermaidRef.value.innerHTML = `
        <div class="p-4 border border-red-300 rounded bg-red-50 text-red-700">
          <strong>Mermaid Diagram Error:</strong>
          <pre class="mt-2 text-sm">${error}</pre>
        </div>
      `
    }
  }
}

onMounted(() => {
  if (isMermaid.value) {
    renderMermaid()
  }
})
</script>

<!-- ▼ 先頭改行が出ないよう <slot/> を 1 行目に直で置く -->
<template>
  <!-- Mermaid diagram rendering -->
  <div v-if="isMermaid" class="my-6 mermaid-container">
    <div ref="mermaidRef" class="mermaid-diagram flex justify-center items-center min-h-[200px] p-4 bg-white rounded-lg border"></div>
  </div>
  
  <!-- Regular code block -->
  <div v-else class="relative my-6 group/not-prose">
    <pre
      class="overflow-x-auto rounded-lg p-4 text-sm leading-relaxed
             bg-zinc-900/90 text-zinc-100">
      <code ref="codeRef" class="block"><slot /></code>
    </pre>

    <button
      @click="copy"
      class="absolute top-2 right-2 flex items-center gap-1 rounded-md
             bg-zinc-700/70 hover:bg-zinc-700 text-xs text-white px-2 py-1
             opacity-0 group-hover/not-prose:opacity-100 transition-opacity">
      <span v-if="copied">Copied ✓</span>
      <span v-else>Copy</span>
    </button>
  </div>
</template>

<style scoped>
/* Safari スクロールバーを細く */
pre::-webkit-scrollbar { height: 8px }
pre::-webkit-scrollbar-thumb { background: rgba(255,255,255,.25); border-radius: 4px }

/* Mermaid diagram styling */
.mermaid-container {
  @apply overflow-x-auto;
}

.mermaid-diagram {
  @apply w-full;
}

.mermaid-diagram :deep(svg) {
  @apply max-w-full h-auto;
}

/* Dark mode support for Mermaid */
.dark .mermaid-diagram {
  @apply bg-gray-800 border-gray-700;
}

.dark .mermaid-diagram :deep(svg) {
  filter: invert(0.9) hue-rotate(180deg);
}

/* Responsive design */
@media (max-width: 768px) {
  .mermaid-diagram :deep(svg) {
    transform: scale(0.9);
    transform-origin: center top;
  }
}
</style>
