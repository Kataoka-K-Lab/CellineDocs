<script setup lang="ts">
import SidebarNavigation from '@/components/SidebarNavigation.vue'
const { page } = useContent()
</script>

<template>
  <div class="flex min-h-screen">
    <!-- ❶ サイドバー：ここだけで ContentNavigation を呼ぶ -->
    <aside class="w-72 shrink-0 border-r border-zinc-200 dark:border-zinc-800
                   bg-white dark:bg-zinc-900">
      <div class="sticky top-0 p-6">
        <ContentNavigation v-slot="{ navigation }">
          <SidebarNavigation :items="navigation" />
        </ContentNavigation>
      </div>
    </aside>

    <!-- ❷ メインエリア -->
    <main class="flex-1 min-w-0 flex flex-col">
      <header class="sticky top-0 z-20 flex items-center gap-2 bg-white/90 dark:bg-zinc-900/90 backdrop-blur
                     border-b border-zinc-200 dark:border-zinc-800 px-6 py-3">
        <NuxtLink to="/" class="font-semibold text-primary-600">Docs</NuxtLink>
        <span class="text-zinc-400">/</span>
        <span class="font-medium">{{ page.title }}</span>
      </header>

      <article class="prose lg:prose-lg dark:prose-invert max-w-none px-6 py-8">
        <slot />
      </article>
    </main>
  </div>
</template>
