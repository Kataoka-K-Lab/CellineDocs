export default defineNuxtPlugin(async () => {
  // Only run on client side
  if (process.server) return

  try {
    // Import mermaid dynamically
    const mermaid = await import('mermaid')
    
    // Initialize mermaid with default configuration
    mermaid.default.initialize({
      startOnLoad: false,
      theme: 'default',
      securityLevel: 'loose',
      flowchart: {
        useMaxWidth: true,
        htmlLabels: true
      },
      sequence: {
        useMaxWidth: true,
        wrap: true
      },
      gantt: {
        useMaxWidth: true
      },
      journey: {
        useMaxWidth: true
      },
      timeline: {
        useMaxWidth: true
      },
      themeVariables: {
        primaryColor: '#3b82f6',
        primaryTextColor: '#1f2937',
        primaryBorderColor: '#3b82f6',
        lineColor: '#6b7280',
        sectionBkgColor: '#f9fafb',
        altSectionBkgColor: '#ffffff',
        gridColor: '#e5e7eb',
        secondaryColor: '#e5e7eb',
        tertiaryColor: '#f3f4f6',
        background: '#ffffff',
        mainBkg: '#ffffff',
        secondBkg: '#f9fafb',
        tertiaryBkg: '#f3f4f6'
      }
    })

    // Make mermaid available globally
    return {
      provide: {
        mermaid: mermaid.default
      }
    }
  } catch (error) {
    console.warn('Failed to initialize Mermaid:', error)
  }
})