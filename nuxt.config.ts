export default defineNuxtConfig({
  modules: ['@nuxt/content', '@nuxtjs/tailwindcss', '@nuxtjs/color-mode'],
  css: ['katex/dist/katex.min.css'],
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
