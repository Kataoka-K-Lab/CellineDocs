<script setup lang="ts">
const props = defineProps<{ id?: string }>()
</script>

<template>
  <h2
    :id="props.id"
    class="relative group text-3xl font-bold tracking-tight text-zinc-900 dark:text-zinc-100 mt-12 mb-4 border-b border-zinc-200 dark:border-zinc-700 pb-1 scroll-mt-24">
    <a
      v-if="props.id"
      :href="`#${props.id}`"
      class="absolute -left-6 top-1/2 -translate-y-1/2 opacity-0 group-hover:opacity-100 transition-opacity"
      aria-label="Anchor"
    >
      <span aria-hidden="true">#</span>
    </a>
    <slot />
  </h2>
</template>
