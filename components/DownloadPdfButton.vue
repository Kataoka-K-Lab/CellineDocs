<script setup lang="ts">
/* クリックで PDF を生成 */
async function download () {
  const html2pdf = (await import('html2pdf.js')).default

  /* ───────── 元記事をクローン ───────── */
  const origin = document.querySelector('article')
  if (!origin) return

  const clone = origin.cloneNode(true) as HTMLElement

  // 1️⃣ dark クラスを除外してライト配色に
  clone.classList.remove('dark:prose-invert')

  // 2️⃣ 白背景 & 黒文字を保証（ダーク時の text-white 等を上書き）
  clone.style.background = '#ffffff'
  clone.style.color      = '#000000'

  /* ───────── PDF 生成 ───────── */
  await html2pdf()
    .set({
      margin: [12, 12],                                  // mm
      filename: (document.title || 'doc') + '.pdf',
      html2canvas: { scale: 2, useCORS: true },
      jsPDF: { unit: 'mm', format: 'a4', orientation: 'portrait' }
    })
    .from(clone)
    .toPdf()                      // ← PDF 化完了後に…
    .get('pdf')                   //    jsPDF インスタンス取得
    .then((pdf) => {
      /* ───── ページ番号を挿入 ───── */
      const total = pdf.internal.getNumberOfPages()
      const width = pdf.internal.pageSize.getWidth()
      const height = pdf.internal.pageSize.getHeight()
      pdf.setFontSize(10)

      for (let i = 1; i <= total; i++) {
        pdf.setPage(i)
        pdf.text(`${i} / ${total}`, width - 20, height - 10, { align: 'right' })
      }
    })
    .save()
}
</script>

<template>
  <button
    @click="download"
    class="flex items-center gap-1 rounded-md px-2.5 py-1.5 text-xs
           bg-zinc-100 hover:bg-zinc-200 dark:bg-zinc-800 dark:hover:bg-zinc-700
           border border-zinc-300 dark:border-zinc-600">
    <!-- ダウンロードアイコン -->
    <svg class="h-4 w-4" viewBox="0 0 24 24" fill="none" stroke="currentColor">
      <path d="M12 16V8" stroke-width="2" stroke-linecap="round"/>
      <path d="M9 13l3 3 3-3" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>
      <rect x="4" y="4" width="16" height="18" rx="2" ry="2"
            stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>
    </svg>
    PDF
  </button>
</template>
