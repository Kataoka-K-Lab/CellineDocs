<script setup lang="ts">
import { useColorMode } from '#imports'

const colorMode = useColorMode()

function toggle () {
  colorMode.preference = colorMode.value === 'dark' ? 'light' : 'dark'
}
</script>

<template>
  <button
    @click="toggle"
    class="fixed bottom-4 right-4 z-50 flex items-center justify-center
           h-10 w-10 rounded-full shadow-lg
           bg-zinc-800 text-white dark:bg-zinc-200 dark:text-zinc-800
           ring-1 ring-zinc-700/10 dark:ring-zinc-900/20
           hover:scale-105 active:scale-95 transition-all"
    aria-label="Toggle color scheme"
  >
    <!-- 🌙 / ☀ 切替 -->
    <span v-if="colorMode.value === 'dark'">☀︎</span>
    <span v-else>☾</span>
  </button>
</template>
