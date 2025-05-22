<template>
  <a
    class="text-primary-600 underline decoration-primary-600/30 hover:text-primary-500"
    :class="[
      'rounded focus-visible:outline-none focus-visible:ring-2',
      'focus-visible:ring-primary-500 focus-visible:ring-offset-2',
      'focus-visible:ring-offset-white dark:focus-visible:ring-offset-zinc-900',
      'transition-colors'
    ]"
  >
    <slot />
  </a>
</template>