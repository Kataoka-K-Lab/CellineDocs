<script setup lang="ts">
import { computed } from 'vue'

const props = defineProps<{ href?: string }>()
const isExternal = computed(() => props.href?.startsWith('http'))
</script>

<template>
  <a
    :href="props.href"
    :target="isExternal ? '_blank' : undefined"
    :rel="isExternal ? 'noopener noreferrer' : undefined"
    class="text-primary-600 underline decoration-primary-600/30 hover:text-primary-500"
    :class="[
      'rounded focus-visible:outline-none focus-visible:ring-2',
      'focus-visible:ring-primary-500 focus-visible:ring-offset-2',
      'focus-visible:ring-offset-white dark:focus-visible:ring-offset-zinc-900',
      'transition-colors'
    ]"
  >
    <slot />
  </a>
</template>