/** @type {import('tailwindcss').Config} */
export default {
  darkMode: 'class',
  content: ['./index.html', './src/**/*.{js,ts,jsx,tsx}'],
  theme: {
    extend: {
      fontFamily: {
        sans: [
          '-apple-system',
          'BlinkMacSystemFont',
          '"SF Pro Display"',
          '"SF Pro Text"',
          '"Inter"',
          'system-ui',
          'sans-serif',
        ],
        mono: [
          '"SF Mono"',
          'ui-monospace',
          'SFMono-Regular',
          'Menlo',
          'Monaco',
          'Consolas',
          'monospace',
        ],
      },
      colors: {
        apple: {
          blue: '#007AFF',
          blueDark: '#0A84FF',
          green: '#34C759',
          greenDark: '#30D158',
          orange: '#FF9500',
          orangeDark: '#FF9F0A',
          red: '#FF3B30',
          redDark: '#FF453A',
          purple: '#AF52DE',
          purpleDark: '#BF5AF2',
          teal: '#5AC8FA',
          indigo: '#5856D6',
        },
      },
      boxShadow: {
        'apple-sm': '0 1px 2px rgba(0, 0, 0, 0.04), 0 0 0 0.5px rgba(0, 0, 0, 0.06)',
        'apple-md': '0 4px 20px rgba(0, 0, 0, 0.06), 0 1px 3px rgba(0, 0, 0, 0.04)',
        'apple-lg': '0 12px 40px rgba(0, 0, 0, 0.12), 0 2px 8px rgba(0, 0, 0, 0.06)',
        'apple-glass': 'inset 0 0.5px 0 rgba(255, 255, 255, 0.6), 0 8px 32px rgba(0, 0, 0, 0.08)',
      },
      backdropBlur: {
        apple: '24px',
      },
    },
  },
  plugins: [],
};
