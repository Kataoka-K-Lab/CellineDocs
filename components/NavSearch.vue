<script setup lang="ts">
import { ref, onMounted } from 'vue'
import { watchDebounced, useMagicKeys } from '@vueuse/core'
import { useAsyncData } from '#app'
import { queryContent } from '#imports'

/* ───── ステート ───── */
const q = ref('')
const inputEl = ref<HTMLInputElement>()

/* ───── ⌘K / Ctrl+K でフォーカス ───── */
const { meta_k, ctrl_k } = useMagicKeys()
onMounted(() => {
  watchDebounced([meta_k, ctrl_k], ([isMetaK, isCtrlK]) => {
    if (isMetaK || isCtrlK) inputEl.value?.focus()
  })
})

/* ───── 検索 ───── */
const { data: hits, refresh } = useAsyncData(
  () => 'search-' + q.value,                 // 🔑 可変キーでキャッシュ無効化
  () =>
    q.value.trim()
      ? queryContent()
          .where({
            $or: [
              { title:  { $contains: q.value } },  // ← ここを $contains に
              { _path:  { $contains: q.value } }
            ]
          })
          .only(['_path', 'title'])
          .limit(15)
          .find()
      : []                                       // 空検索は空配列
)

/* 250 ms デバウンスでリフレッシュ */
watchDebounced(q, () => refresh(), { debounce: 250 })
</script>

<template>
  <div class="relative">
    <!-- 入力欄 -->
    <input
      ref="inputEl"
      v-model="q"
      type="search"
      placeholder="Search"
      class="w-full rounded-md pl-9 pr-16 py-2 text-sm
             border border-zinc-300 dark:border-zinc-700
             bg-white/90 dark:bg-zinc-800/70 backdrop-blur
             focus:ring-2 focus:ring-blue-500 outline-none"
    />

    <!-- 左：虫めがね -->
    <svg class="absolute left-3 top-2.5 h-4 w-4 text-zinc-400 pointer-events-none"
         viewBox="0 0 24 24" fill="none" stroke="currentColor">
      <circle cx="11" cy="11" r="8" stroke-width="2"/>
      <path d="M21 21l-5.2-5.2" stroke-width="2" stroke-linecap="round"/>
    </svg>

    <!-- 右：ショートカット表示 -->
    <span class="absolute right-3 top-1/2 -translate-y-1/2 text-[11px]
                 text-zinc-500 dark:text-zinc-400 select-none flex items-center gap-0.5">
      <kbd class="rounded border px-1 py-0.5 text-[10px]
                 border-zinc-300 dark:border-zinc-600 bg-zinc-100 dark:bg-zinc-700">⌘</kbd>
      <kbd class="rounded border px-1 py-0.5 text-[10px]
                 border-zinc-300 dark:border-zinc-600 bg-zinc-100 dark:bg-zinc-700">K</kbd>
    </span>

    <!-- サジェスト -->
    <ul v-if="q && hits?.length"
        class="absolute z-30 mt-1 w-full max-h-60 overflow-auto
               rounded-md border border-zinc-200 dark:border-zinc-700
               bg-white dark:bg-zinc-800 shadow-lg">
      <li v-for="hit in hits" :key="hit._path">
        <NuxtLink
          :to="hit._path"
          class="block px-3 py-2 text-sm truncate
                 hover:bg-blue-50 dark:hover:bg-blue-900/30"
          @click="q = ''">
          {{ hit.title }}
        </NuxtLink>
      </li>
    </ul>
  </div>
</template>
