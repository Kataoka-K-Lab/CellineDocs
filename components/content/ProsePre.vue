<script setup lang="ts">
import { ref } from 'vue'

/* Nuxt Content が差し込む <code class="language-js"> … */
const codeRef = ref<HTMLElement>()
const copied  = ref(false)

async function copy () {
  if (!codeRef.value) return
  await navigator.clipboard.writeText(codeRef.value.innerText)
  copied.value = true
  setTimeout(() => { copied.value = false }, 1500)
}
</script>

<!-- ▼ 先頭改行が出ないよう <slot/> を 1 行目に直で置く -->
<template>
  <div class="relative my-6 group/not-prose">
    <pre
      class="overflow-x-auto rounded-xl p-4 text-sm leading-relaxed font-mono border border-zinc-800 bg-zinc-900/75 text-zinc-100">
      <code ref="codeRef" class="block"><slot /></code>
    </pre>

    <button
      @click="copy"
      class="absolute top-2 right-2 flex items-center gap-1 rounded-md bg-primary-600/80 hover:bg-primary-600 text-xs text-white px-2 py-1 opacity-0 group-hover/not-prose:opacity-100 transition-opacity focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary-400 focus-visible:ring-offset-2 focus-visible:ring-offset-zinc-900">
      <span v-if="copied">Copied ✓</span>
      <span v-else>Copy</span>
    </button>
  </div>
</template>

<style scoped>
/* Safari スクロールバーを細く */
pre::-webkit-scrollbar { height: 8px }
pre::-webkit-scrollbar-thumb { background: rgba(255,255,255,.25); border-radius: 4px }
</style>
