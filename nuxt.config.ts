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
      preload: ['js','ts','vue','bash', 'python', 'mermaid']
    },
    markdown: {
      remarkPlugins: ['remark-math'],
      rehypePlugins: [
        ['rehype-katex', { output: 'html', strict: false }],
        ['rehype-mermaid', { strategy: 'img-svg' }]
      ],
      anchorLinks: false  // 見出しの自動アンカーリンク生成を無効化
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
