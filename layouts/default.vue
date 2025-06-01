<script setup lang="ts">
import { ref } from 'vue'
import SidebarNavigation from '@/components/SidebarNavigation.vue'
const { page } = useContent()

/* モバイル時に左ドロワーを開閉するフラグ */
const open = ref(false)
</script>

<template>
  <ContentNavigation v-slot="{ navigation }">
    <div class="flex min-h-screen">
      <!-- ───────── モバイル・ドロワー ───────── -->
      <transition
        enter-from-class="-translate-x-full"
        enter-active-class="transform transition-transform duration-200"
        leave-to-class="-translate-x-full"
        leave-active-class="transform transition-transform duration-200">
        <aside
          v-show="open"
          class="fixed inset-y-0 left-0 z-40 w-64 bg-white dark:bg-zinc-900
                 border-r border-zinc-200 dark:border-zinc-800
                 p-6 overflow-y-auto lg:hidden">
        <div class="sticky top-0 p-6 space-y-7">
          <NuxtLink to="/" class="flex items-center gap-2">
            <img src="/celline_icon.png"
                alt="Celline"
                class="h-auto w-30 dark:invert" />
          </NuxtLink>
          <NavSearch />
          <SidebarNavigation :items="navigation" @click="open = false" />
        </div>
        </aside>
      </transition>

      <!-- ───────── デスクトップ用サイドバー ───────── -->
      <aside class="hidden lg:block w-72 shrink-0 border-r
                    border-zinc-200 dark:border-zinc-800
                    bg-white dark:bg-zinc-900">
        <div class="sticky top-0 p-6 space-y-7">
          <NuxtLink to="/" class="flex items-center gap-2">
            <img src="/celline_icon.png"
                alt="Celline"
                class="h-auto w-30 dark:invert" />
          </NuxtLink>
          <NavSearch />
          <SidebarNavigation :items="navigation" />
        </div>
      </aside>

      <!-- ───────── メインエリア ───────── -->
      <main class="flex-1 min-w-0 flex flex-col">
        <header
          class="sticky top-0 z-30 flex items-center gap-3
                 bg-white/90 dark:bg-zinc-900/90 backdrop-blur
                 border-b border-zinc-200 dark:border-zinc-800 px-4 py-3">
          <!-- ハンバーガー：lg 以上で非表示 -->
          <button
            class="lg:hidden h-10 w-10 flex items-center justify-center rounded-md
                   hover:bg-zinc-100 dark:hover:bg-zinc-800"
            @click="open = !open"
            aria-label="Open navigation">
            <svg class="h-5 w-5" viewBox="0 0 24 24" fill="none" stroke="currentColor">
              <path d="M4 6h16M4 12h16M4 18h16" stroke-width="2"
                    stroke-linecap="round" stroke-linejoin="round" />
            </svg>
          </button>

          <span class="text-zinc-400">/</span>
          <span class="font-medium truncate">{{ page?.title || 'Docs' }}</span>

          <DownloadPdfButton />

        </header>

        <article class="prose lg:prose-lg dark:prose-invert max-w-none px-4 py-8">
          <slot />
        </article>
      </main>
    </div>
    <ColorModeToggle />
  </ContentNavigation>
</template>
