<script setup lang="ts">
const props = defineProps<{ id?: string }>()
</script>

<template>
  <h1
    :id="props.id"
    class="relative group text-4xl font-extrabold tracking-tight text-zinc-900 dark:text-zinc-100 mb-6 border-b border-zinc-200 dark:border-zinc-700 pb-2 scroll-mt-28">
    <a
      v-if="props.id"
      :href="`#${props.id}`"
      class="absolute -left-6 top-1/2 -translate-y-1/2 opacity-0 group-hover:opacity-100 transition-opacity"
      aria-label="Anchor"
    >
      <span aria-hidden="true">#</span>
    </a>
    <slot />
  </h1>
</template>
