<script setup lang="ts">
const props = defineProps<{ items: any[] }>()
const route = useRoute()

function isActive (item: any) {
  return route.path === item._path
}
function isInTrail (item: any) {
  return route.path.startsWith(item._path + '/') || isActive(item)
}
</script>

<template>
  <ul class="space-y-0.5 text-[15px] leading-[22px]">
    <li v-for="item in props.items" :key="item._path">
      <!-- 親ノード -->
      <details v-if="item.children?.length"
               :open="isInTrail(item)"
               class="group">
        <summary
          class="flex items-center gap-2 cursor-pointer select-none
                 px-3 py-1.5 rounded-md
                 transition-colors hover:bg-zinc-100 dark:hover:bg-zinc-800"
          :class="isInTrail(item) && 'text-blue-600 font-semibold'"
        >
          <!-- 矢印 -->
          <svg class="h-7 w-7 text-zinc-500 -rotate-90 group-open:rotate-0
                      transition-transform duration-200"
               viewBox="0 0 20 20" fill="currentColor">
            <path d="M6 6l4 4 4-4" fill-rule="evenodd" clip-rule="evenodd"/>
          </svg>

          <span v-if="!item._path">{{ item.title }}</span>
          <NuxtLink v-else :to="item._path" class="flex-1 truncate">
            {{ item.title }}
          </NuxtLink>
        </summary>

        <div v-show="true" class="pl-5 mt-0.5 border-l
                       border-zinc-200 dark:border-zinc-700">
          <SidebarNavigation :items="item.children" />
        </div>
      </details>

      <!-- リーフ -->
      <NuxtLink
        v-else
        :to="item._path"
        class="flex items-center gap-2 px-3 py-1.5 rounded-md
               transition-colors hover:bg-zinc-100 dark:hover:bg-zinc-800
               relative"
        :class="isActive(item)
          ? 'bg-blue-50 dark:bg-blue-900/40 text-blue-600 font-semibold'
          : ''"
      >
        <!-- 縦バー（active の時）-->
        <span v-if="isActive(item)"
              class="absolute left-0 top-0 h-full w-[3px] rounded-r
                     bg-blue-600 dark:bg-blue-500"/>
        {{ item.title }}
      </NuxtLink>
    </li>
  </ul>
</template>

<style scoped>
details > summary::-webkit-details-marker { display: none; }
</style>
