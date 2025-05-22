// // https://nuxt.com/docs/guide/directory-structure/nuxt.config#nuxt-config-file
// export default defineNuxtConfig({
//   modules: [
//     '@nuxt/content',
//     '@nuxt/eslint',
//     '@nuxt/fonts',
//     '@nuxt/icon',
//     '@nuxt/image',
//     '@nuxt/scripts',
//     '@nuxt/test-utils',
//     '@nuxt/ui',
//     '@nuxtjs/tailwindcss',
//   ],
//   css: [
//     '~/assets/css/tailwind.css',
//     '~/assets/scss/main.scss'
//   ],
//   content: {
//     documentDriven: true,
//     highlight: {
//       theme: {
//         default: 'github-light',
//         dark: 'github-dark'
//       }
//     }
//   },
//   tailwindcss: {
//     viewer: false
//   },

//   app: {
//     head: {
//       meta: [
//         { name: 'viewport', content: 'width=device-width, initial-scale=1' },
//         { name: 'theme-color', content: '#2563eb' },
//         { property: 'og:site_name', content: 'Docs' }
//       ],
//       link: [
//         { rel: 'icon', type: 'image/png', href: '/favicon.png' }
//       ]
//     }
//   },
//   ui: {
//     primary: '#2563eb'
//   }
// })
import { existsSync } from 'node:fs'

export default defineNuxtConfig({
  modules: ['@nuxt/content', '@nuxtjs/tailwindcss', '@nuxtjs/color-mode'],
  css: [
    '~/assets/css/tailwind.css',
    ...(existsSync('node_modules/katex/dist/katex.min.css') ? ['katex/dist/katex.min.css'] : [])
  ],
  content: {
     documentDriven: true,
     highlight: {
      theme: {              // ダーク／ライト別テーマ
        default: 'github-dark',
        dark: 'github-dark'
      },
      preload: ['js','ts','vue','bash', 'python']
    },
    markdown: {
      remarkPlugins: ['remark-math'],
      rehypePlugins: [
        ['rehype-katex', { output: 'html', strict: false }]
      ]
    },
  },
  colorMode: {
    classSuffix: ''               // ⚠️ 'dark' / 'light' だけを付与
  },
  /* 静的サイト生成モード */
  nitro: { preset: 'static' },           // nuxi generate → .output/public

  app: {
    /* リポジトリ名に合わせて必ず書き換えてください */
    baseURL: process.env.NODE_ENV === 'production'
      ? '/CellineDocs/'                  // 例: '/celline-docs/'
      : '/'
  }
})
