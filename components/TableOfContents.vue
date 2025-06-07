<template>
  <!-- Desktop TOC (fixed position) -->
  <div 
    v-if="tocItems.length > 0" 
    class="hidden xl:block fixed right-4 top-1/2 transform -translate-y-1/2 max-w-xs w-64 z-50"
  >
    <div class="bg-white/80 dark:bg-gray-800/80 backdrop-blur-sm border border-gray-200 dark:border-gray-700 rounded-lg p-4 shadow-lg max-h-96 overflow-y-auto">
      <!-- Toggle Button -->
      <div class="flex items-center justify-between mb-3">
        <h3 class="text-sm font-semibold text-gray-900 dark:text-gray-100">
          Table of Contents
        </h3>
        <button 
          @click="isCollapsed = !isCollapsed"
          class="text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200 p-1"
          :title="isCollapsed ? 'Expand TOC' : 'Collapse TOC'"
        >
          <svg 
            class="w-4 h-4 transition-transform duration-200" 
            :class="{ 'rotate-180': isCollapsed }"
            fill="none" 
            stroke="currentColor" 
            viewBox="0 0 24 24"
          >
            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19 9l-7 7-7-7" />
          </svg>
        </button>
      </div>

      <!-- TOC Items -->
      <nav v-if="!isCollapsed" class="space-y-1">
        <a
          v-for="item in tocItems"
          :key="item.id"
          :href="`#${item.id}`"
          @click="scrollToHeading(item.id)"
          class="block py-1 px-2 text-sm rounded transition-colors duration-200"
          :class="[
            item.level === 2 ? 'text-gray-900 dark:text-gray-100 font-medium' : '',
            item.level === 3 ? 'text-gray-700 dark:text-gray-300 ml-3' : '',
            item.level === 4 ? 'text-gray-600 dark:text-gray-400 ml-6 text-xs' : '',
            item.level === 5 ? 'text-gray-500 dark:text-gray-500 ml-9 text-xs' : '',
            item.level === 6 ? 'text-gray-400 dark:text-gray-600 ml-12 text-xs' : '',
            'hover:bg-gray-100 dark:hover:bg-gray-700',
            activeId === item.id ? 'bg-blue-100 dark:bg-blue-900 text-blue-700 dark:text-blue-300' : ''
          ]"
        >
          {{ item.text }}
        </a>
      </nav>
    </div>
  </div>

  <!-- Mobile TOC Button (visible on smaller screens) -->
  <div v-if="tocItems.length > 0" class="xl:hidden fixed bottom-4 right-4 z-50">
    <button
      @click="showMobileToc = !showMobileToc"
      class="bg-blue-600 hover:bg-blue-700 text-white p-3 rounded-full shadow-lg transition-colors duration-200"
      :title="showMobileToc ? 'Hide TOC' : 'Show TOC'"
    >
      <svg class="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
        <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 6h16M4 10h16M4 14h16M4 18h16" />
      </svg>
    </button>

    <!-- Mobile TOC Modal -->
    <div
      v-if="showMobileToc"
      class="fixed inset-0 bg-black/50 flex items-center justify-center p-4 z-60"
      @click="showMobileToc = false"
    >
      <div
        class="bg-white dark:bg-gray-800 rounded-lg p-6 max-w-sm w-full max-h-96 overflow-y-auto"
        @click.stop
      >
        <div class="flex items-center justify-between mb-4">
          <h3 class="text-lg font-semibold text-gray-900 dark:text-gray-100">
            Table of Contents
          </h3>
          <button
            @click="showMobileToc = false"
            class="text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200 p-1"
          >
            <svg class="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M6 18L18 6M6 6l12 12" />
            </svg>
          </button>
        </div>

        <nav class="space-y-2">
          <a
            v-for="item in tocItems"
            :key="item.id"
            :href="`#${item.id}`"
            @click="scrollToHeading(item.id); showMobileToc = false"
            class="block py-2 px-3 text-sm rounded transition-colors duration-200"
            :class="[
              item.level === 2 ? 'text-gray-900 dark:text-gray-100 font-medium' : '',
              item.level === 3 ? 'text-gray-700 dark:text-gray-300 ml-4' : '',
              item.level === 4 ? 'text-gray-600 dark:text-gray-400 ml-8 text-xs' : '',
              item.level === 5 ? 'text-gray-500 dark:text-gray-500 ml-12 text-xs' : '',
              item.level === 6 ? 'text-gray-400 dark:text-gray-600 ml-16 text-xs' : '',
              'hover:bg-gray-100 dark:hover:bg-gray-700',
              activeId === item.id ? 'bg-blue-100 dark:bg-blue-900 text-blue-700 dark:text-blue-300' : ''
            ]"
          >
            {{ item.text }}
          </a>
        </nav>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, onMounted, onUnmounted, nextTick, watch } from 'vue'

interface TocItem {
  id: string
  text: string
  level: number
}

const tocItems = ref<TocItem[]>([])
const activeId = ref<string>('')
const isCollapsed = ref(false)
const showMobileToc = ref(false)

// Function to generate TOC from page headings
const generateTOC = () => {
  const headings = document.querySelectorAll('h2, h3, h4, h5, h6')
  const items: TocItem[] = []

  headings.forEach((heading, index) => {
    const level = parseInt(heading.tagName.charAt(1))
    const text = heading.textContent?.trim() || ''
    
    // Generate ID if not present
    let id = heading.id
    if (!id) {
      id = text
        .toLowerCase()
        .replace(/[^\w\s-]/g, '') // Remove special characters
        .replace(/\s+/g, '-') // Replace spaces with hyphens
        .replace(/--+/g, '-') // Replace multiple hyphens with single
        .trim()
      
      // Ensure uniqueness
      let finalId = id
      let counter = 1
      while (document.getElementById(finalId)) {
        finalId = `${id}-${counter}`
        counter++
      }
      
      heading.id = finalId
      id = finalId
    }

    items.push({
      id,
      text,
      level
    })
  })

  tocItems.value = items
}

// Function to scroll to heading smoothly
const scrollToHeading = (id: string) => {
  const element = document.getElementById(id)
  if (element) {
    element.scrollIntoView({ 
      behavior: 'smooth',
      block: 'start'
    })
  }
}

// Function to update active heading based on scroll position
const updateActiveHeading = () => {
  const headings = document.querySelectorAll('h2, h3, h4, h5, h6')
  const scrollTop = window.scrollY + 100 // Offset for better UX

  let currentActiveId = ''

  headings.forEach((heading) => {
    const rect = heading.getBoundingClientRect()
    const offsetTop = rect.top + window.scrollY

    if (offsetTop <= scrollTop) {
      currentActiveId = heading.id
    }
  })

  activeId.value = currentActiveId
}

// Throttle function for better performance
const throttle = (func: Function, limit: number) => {
  let inThrottle: boolean
  return function(this: any) {
    const args = arguments
    const context = this
    if (!inThrottle) {
      func.apply(context, args)
      inThrottle = true
      setTimeout(() => inThrottle = false, limit)
    }
  }
}

const throttledUpdateActiveHeading = throttle(updateActiveHeading, 100)

onMounted(async () => {
  // Wait for content to be rendered
  await nextTick()
  
  // Small delay to ensure all content is fully rendered
  setTimeout(() => {
    generateTOC()
    updateActiveHeading()
  }, 500)

  // Add scroll listener
  window.addEventListener('scroll', throttledUpdateActiveHeading)

  // Re-generate TOC when content changes (for dynamic content)
  const observer = new MutationObserver(() => {
    setTimeout(generateTOC, 100)
  })

  observer.observe(document.body, {
    childList: true,
    subtree: true
  })

  // Store observer for cleanup
  ;(window as any).__tocObserver = observer
})

onUnmounted(() => {
  window.removeEventListener('scroll', throttledUpdateActiveHeading)
  
  // Cleanup mutation observer
  if ((window as any).__tocObserver) {
    ;(window as any).__tocObserver.disconnect()
    delete (window as any).__tocObserver
  }
})

// Re-generate TOC when route changes
watch(() => useRoute().path, () => {
  nextTick(() => {
    setTimeout(() => {
      generateTOC()
      updateActiveHeading()
    }, 500)
  })
})
</script>

<style scoped>
/* Custom scrollbar for TOC */
.overflow-y-auto::-webkit-scrollbar {
  width: 4px;
}

.overflow-y-auto::-webkit-scrollbar-track {
  background: transparent;
}

.overflow-y-auto::-webkit-scrollbar-thumb {
  background: rgba(156, 163, 175, 0.5);
  border-radius: 2px;
}

.overflow-y-auto::-webkit-scrollbar-thumb:hover {
  background: rgba(156, 163, 175, 0.8);
}

.dark .overflow-y-auto::-webkit-scrollbar-thumb {
  background: rgba(75, 85, 99, 0.5);
}

.dark .overflow-y-auto::-webkit-scrollbar-thumb:hover {
  background: rgba(75, 85, 99, 0.8);
}
</style>