<script setup lang="ts">
const props = defineProps<{ items: any[] }>()
const route = useRoute()
</script>

<template>
  <ul class="space-y-1">
    <li v-for="item in props.items" :key="item._path">
      <!-- ── 親ノード ───────────────────────── -->
      <details
        v-if="item.children?.length"
        :open="route.path.startsWith(item._path)"
        class="group"
      >
        <!-- summary 部をクリックすると rotate -->
        <summary
          class="flex items-center cursor-pointer select-none font-semibold
                 transition-colors hover:text-blue-600"
        >
        <svg
          class="h-6 w-6 mr-1 transition-transform duration-200
                -rotate-90   /* ← 初期状態で右向きになる */
                group-open:rotate-0"
          viewBox="0 0 20 20" fill="currentColor"
        >
          <!-- パスは下向き矢印のまま -->
          <path fill-rule="evenodd" d="M6 6l4 4 4-4" clip-rule="evenodd" />
        </svg>

          <span v-if="!item._path">{{ item.title }}</span>
          <NuxtLink v-else :to="item._path" class="flex-1 hover:underline">
            {{ item.title }}
          </NuxtLink>
        </summary>

        <!-- 子ノード：高さでスライド -->
        <transition
          enter-active-class="overflow-hidden animate-acc-down"
          leave-active-class="overflow-hidden animate-acc-up"
        >
          <div v-show="true" class="ml-3 pl-1 border-l border-zinc-200 dark:border-zinc-700">
            <SidebarNavigation :items="item.children" />
          </div>
        </transition>
      </details>

      <!-- ── 子ノード（葉）────────────────────── -->
      <NuxtLink
        v-else
        :to="item._path"
        class="block pl-5 py-0.5 rounded-md
               transition-colors hover:bg-zinc-100 dark:hover:bg-zinc-800
               "
        :class="route.path === item._path && 'bg-blue-50 text-blue-600 dark:bg-blue-900/40'"
      >
        {{ item.title }}
      </NuxtLink>
    </li>
  </ul>
</template>

<style scoped>
details > summary::-webkit-details-marker { display: none; }
</style>
