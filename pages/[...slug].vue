<script setup lang="ts">
/*   👇 どのレイアウトを当てるかは自由。default でも content でも OK */
definePageMeta({ layout: 'default' })
</script>

<template>
  <!-- navigation と toc をオフにするのがポイント -->
  <ContentDoc :navigation="false" :toc="false" />
</template>
