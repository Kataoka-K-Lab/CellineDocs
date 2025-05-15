<script setup lang="ts">
const props = defineProps<{ items: any[] }>()
const route = useRoute()
</script>

<template>
  <ul class="space-y-1">
    <li v-for="item in props.items" :key="item._path">
      <details :open="route.path.startsWith(item._path)" v-if="item.children?.length">
        <summary class="cursor-pointer font-semibold">
          <span v-if="!item._path">{{ item.title }}</span>
          <NuxtLink v-else :to="item._path" class="hover:underline">{{ item.title }}</NuxtLink>
        </summary>
        <div class="ml-4 mt-1">
          <SidebarNavigation :items="item.children" />
        </div>
      </details>

      <NuxtLink
        v-else
        :to="item._path"
        class="block pl-4"
        :class="route.path === item._path && 'text-primary-600 font-medium'"
      >
        {{ item.title }}
      </NuxtLink>
    </li>
  </ul>
</template>

<style scoped>
details > summary::-webkit-details-marker { display: none; }
</style>
