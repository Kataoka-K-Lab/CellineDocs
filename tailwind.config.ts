import type { Config } from 'tailwindcss'
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
            colors: {
                primary: {
                    50: '#eff6ff',
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
            typography: ({ theme }) => ({
                DEFAULT: {
                    css: {
                        '--tw-prose-body': theme('colors.zinc.800'),
                        '--tw-prose-headings': theme('colors.zinc.900'),
                        '--tw-prose-links': theme('colors.primary.600'),
                        '--tw-prose-bold': theme('colors.zinc.900'),
                        '--tw-prose-code': theme('colors.zinc.800'),
                        '--tw-prose-hr': theme('colors.zinc.200'),
                        '--tw-prose-quotes': theme('colors.zinc.900'),
                        '--tw-prose-quote-borders': theme('colors.zinc.200'),
                        '--tw-prose-captions': theme('colors.zinc.500'),
                        '--tw-prose-kbd': theme('colors.zinc.900'),
                        '--tw-prose-kbd-shadows': theme('colors.zinc.400')
                    }
                },
                dark: {
                    css: {
                        '--tw-prose-body': theme('colors.zinc.100'),
                        '--tw-prose-headings': theme('colors.zinc.100'),
                        '--tw-prose-links': theme('colors.primary.400'),
                        '--tw-prose-bold': theme('colors.zinc.100'),
                        '--tw-prose-code': theme('colors.zinc.100'),
                        '--tw-prose-hr': theme('colors.zinc.700'),
                        '--tw-prose-quotes': theme('colors.zinc.100'),
                        '--tw-prose-quote-borders': theme('colors.zinc.700'),
                        '--tw-prose-captions': theme('colors.zinc.400'),
                        '--tw-prose-kbd': theme('colors.zinc.100'),
                        '--tw-prose-kbd-shadows': theme('colors.zinc.700')
                    }
                }
            })
        }
    },
    plugins: [
        require('@tailwindcss/typography'),
        plugin(({ addVariant }) => {
            addVariant('dark', '.dark &')
        })
    ]
}