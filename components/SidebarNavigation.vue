<script setup lang="ts">
const props = defineProps<{ items: any[] }>()
const route = useRoute()
</script>

<template>
  <ul class="space-y-1 text-[15px] leading-[22px]">
    <li v-for="item in props.items" :key="item._path">
      <!-- 親ノード -->
      <details
        v-if="item.children?.length"
        :open="route.path.startsWith(item._path)"
        class="group"
      >
        <summary
          class="flex items-center gap-2 cursor-pointer select-none font-medium
                 transition-colors hover:text-blue-600"
        >
          <!-- ▲ 矢印：24px にアップ／開閉で回転 -->
          <svg
            class="h-3 w-3 mr-1 text-zinc-500                     <!-- ← ここだけ変える -->
                  -rotate-90 group-open:rotate-0 transition-transform duration-200"
            viewBox="0 0 20 20" fill="currentColor"
          >
            <path fill-rule="evenodd" d="M6 6l4 4 4-4" clip-rule="evenodd" />
          </svg>

          <span v-if="!item._path">{{ item.title }}</span>
          <NuxtLink v-else :to="item._path" class="flex-1 truncate hover:underline">
            {{ item.title }}
          </NuxtLink>
        </summary>

        <!-- 子階層 24px インデント -->
        <transition
          enter-active-class="overflow-hidden animate-acc-down"
          leave-active-class="overflow-hidden animate-acc-up"
        >
          <div v-show="true" class="pl-6 mt-1 border-l border-zinc-200 dark:border-zinc-700">
            <SidebarNavigation :items="item.children" />
          </div>
        </transition>
      </details>

      <!-- リーフノード -->
      <NuxtLink
        v-else
        :to="item._path"
        class="flex items-center gap-2 pl-6 py-0.5 rounded-md
               transition-colors hover:bg-zinc-100 dark:hover:bg-zinc-800"
        :class="route.path === item._path &&
          'bg-blue-50 text-blue-600 dark:bg-blue-900/40 font-semibold'"
      >
        {{ item.title }}
      </NuxtLink>
    </li>
  </ul>
</template>

<style scoped>
details > summary::-webkit-details-marker { display: none; }
</style>
