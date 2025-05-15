import type { Config } from 'tailwindcss'
import defaultTheme from 'tailwindcss/defaultTheme'
import plugin from 'tailwindcss/plugin'

export default <Partial<Config>>{
  content: [
    './components/**/*.{vue,js,ts}',
    './layouts/**/*.vue',
    './pages/**/*.vue',
    './content/**/*.{md,yml,json}',
    './app.vue'
  ],
  darkMode: 'class',

  theme: {
    extend: {
      /* ---------- カラー & フォント ---------- */
      colors: {
        primary: {
          50:  '#eff6ff',
          100: '#dbeafe',
          200: '#bfdbfe',
          300: '#93c5fd',
          400: '#60a5fa',
          500: '#3b82f6',
          600: '#2563eb',
          700: '#1d4ed8',
          800: '#1e40af',
          900: '#1e3a8a'
        }
      },
      fontFamily: { sans: ['Inter', ...defaultTheme.fontFamily.sans] },

      /* ---------- Tailwind Typography ---------- */
      typography: ({ theme }) => ({
        /* --- ライト --- */
        DEFAULT: {
          css: {
            '--tw-prose-body'          : theme('colors.zinc.800'),
            '--tw-prose-headings'      : theme('colors.zinc.900'),
            '--tw-prose-links'         : theme('colors.primary.600'),
            '--tw-prose-bold'          : theme('colors.zinc.900'),
            '--tw-prose-code'          : theme('colors.zinc.800'),
            '--tw-prose-hr'            : theme('colors.zinc.200'),
            '--tw-prose-quotes'        : theme('colors.zinc.900'),
            '--tw-prose-quote-borders' : theme('colors.zinc.200'),
            '--tw-prose-captions'      : theme('colors.zinc.500'),
            '--tw-prose-kbd'           : theme('colors.zinc.900'),
            '--tw-prose-kbd-shadows'   : theme('colors.zinc.400'),

            /* テーブルのセル余白／罫線を統一 */
            'thead th, tbody td': {
              padding: '0.5rem 0.75rem',
              borderBottomWidth: '1px',
              borderColor: theme('colors.zinc.200')
            }
          }
        },

        /* --- ダーク --- */
        invert: {
          css: {
            '--tw-prose-body'          : theme('colors.zinc.100'),
            '--tw-prose-headings'      : theme('colors.zinc.100'),
            '--tw-prose-links'         : theme('colors.primary.400'),
            '--tw-prose-bold'          : theme('colors.zinc.100'),
            '--tw-prose-code'          : theme('colors.zinc.100'),
            '--tw-prose-hr'            : theme('colors.zinc.700'),
            '--tw-prose-quotes'        : theme('colors.zinc.100'),
            '--tw-prose-quote-borders' : theme('colors.zinc.700'),
            '--tw-prose-captions'      : theme('colors.zinc.400'),
            '--tw-prose-kbd'           : theme('colors.zinc.100'),
            '--tw-prose-kbd-shadows'   : theme('colors.zinc.700'),

            /* ✅ ダークで壊れるテーブルを修正 */
            'thead tr, tbody tr': { display: 'table-row' },
            'thead th': {
              backgroundColor: theme('colors.zinc.800'),
              color: theme('colors.zinc.100'),
              borderColor: theme('colors.zinc.700'),
              padding: '0.5rem 0.75rem'
            },
            'tbody td': {
              color: theme('colors.zinc.100'),
              borderColor: theme('colors.zinc.700'),
              padding: '0.5rem 0.75rem'
            },
            'tbody tr:last-child td': { borderBottomWidth: '0px' }
          }
        }
      }),

      /* ---------- サイドバー用アニメ ---------- */
      keyframes: {
        'acc-down': { from: { height: '0' }, to: { height: 'var(--radix-accordion-content-height)' } },
        'acc-up'  : { from: { height: 'var(--radix-accordion-content-height)' }, to: { height: '0' } }
      },
      animation: {
        'acc-down': 'acc-down 200ms ease-out',
        'acc-up'  : 'acc-up   200ms ease-out'
      }
    }
  },

  plugins: [
    require('@tailwindcss/typography'),
    plugin(({ addVariant }) => addVariant('dark', '.dark &'))
  ]
}
