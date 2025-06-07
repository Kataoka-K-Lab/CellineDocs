# Quick Start Tutorial

このチュートリアルでは、Cellineを使用してシングルセルRNA-seq解析を行う完全なワークフローを実際のデータを使って学習します。

## 🎯 このチュートリアルで学ぶこと

- Cellineプロジェクトの作成と設定
- 公的データベースからのサンプル取得
- データの前処理と品質管理
- インタラクティブな解析とビジュアライゼーション

## 📊 使用データ

このチュートリアルでは、マウス脳細胞のシングルセルRNA-seqデータ（GSE115189）を使用します。

- **データセット**: GSE115189
- **サンプル数**: 2サンプル
- **細胞数**: 約3,000細胞
- **種**: Mus musculus

## 🚀 ステップ1: プロジェクトセットアップ

### プロジェクトディレクトリの作成

```bash
# 作業ディレクトリを作成
mkdir celline-tutorial
cd celline-tutorial

# Cellineプロジェクトを初期化
celline init tutorial-project
```

### 設定の確認

```bash
# インストールが正しく完了しているか確認
celline list

# システム情報を表示
celline info

# 基本設定を行う
celline config
```

設定画面では以下を選択：
- **Execution system**: `multithreading` (ローカル実行)
- **Number of threads**: `2` (お使いのシステムに応じて調整)

## 🔍 ステップ2: データの取得

### サンプルの追加

```bash
# GEOからサンプルセットを追加
celline run add GSE115189

# 個別サンプルを追加する場合
# celline run add GSM3169075 GSM3169076
```

### サンプル情報の確認

```bash
# 追加されたサンプルを確認
cat samples.toml
```

出力例：
```toml
GSM3169075 = "Control sample 1"
GSM3169076 = "Control sample 2"
```

## ⬇️ ステップ3: データダウンロード

### シーケンシングデータのダウンロード

```bash
# 並列ダウンロード（2スレッド）
celline run download --nthread 2

# 進行状況を確認
ls -la resources/
```

### データ構造の確認

```bash
# ダウンロードされたファイルを確認
find resources/ -name "*.fastq.gz" | head -5
```

## 🧮 ステップ4: カウント処理

### Cell Rangerによるカウント

```bash
# 転写情報データベースの設定
celline run set_transcriptome

# Cell Rangerでカウント処理を実行
celline run count
```

この処理には時間がかかります（数時間〜）。完了すると、各サンプルに対して以下が生成されます：
- Feature-barcode matrix
- Quality metrics
- Cell/gene summary

## 🔬 ステップ5: 品質管理と前処理

### Seuratオブジェクトの作成

```bash
# Seuratオブジェクトを作成
celline run create_seurat
```

### 前処理の実行

```bash
# 品質管理と前処理を実行
celline run preprocess

# 結果の確認
ls -la data/*/cell_info.tsv
```

### QCメトリクスの確認

```bash
# 各サンプルのQC統計を確認
head data/GSM3169075/cell_info.tsv
```

出力例：
```tsv
barcode	project	sample	cell	cell_type	include	n_genes_by_counts	pct_counts_mt	doublet_score	predicted_doublets
AAACCTGAGAAGGCCT-1	GSE115189	GSM3169075	GSM3169075_1	Neuron	true	2847	1.2	0.05	false
```

## 🌐 ステップ6: インタラクティブ解析

### Webインターフェースの起動

```bash
# インタラクティブモードを開始
celline interactive
```

ブラウザで http://localhost:8080 にアクセスします。

### インターフェースの機能

1. **Sample Overview**: サンプル一覧と処理状況
2. **Quality Control**: QCメトリクスの可視化
3. **Cell Type Analysis**: 細胞型分布の確認
4. **Gene Expression**: 遺伝子発現パターンの探索

## 📈 ステップ7: 解析結果の確認

### 基本統計の確認

```python
# Pythonスクリプトでの確認例
import polars as pl
import scanpy as sc

# cell_info.tsvの読み込み
cell_info = pl.read_csv("data/GSM3169075/cell_info.tsv", separator="\t")

# 基本統計を表示
print("Cell count:", len(cell_info))
print("Included cells:", cell_info.filter(pl.col("include")).height)
print("Cell types:", cell_info["cell_type"].unique().to_list())
```

### 品質メトリクスの可視化

```python
# QCメトリクスの可視化
import matplotlib.pyplot as plt

# 遺伝子数の分布
plt.figure(figsize=(10, 6))
plt.subplot(1, 2, 1)
plt.hist(cell_info["n_genes_by_counts"], bins=50)
plt.xlabel("Number of genes")
plt.ylabel("Number of cells")

# ミトコンドリア遺伝子比率の分布
plt.subplot(1, 2, 2)
plt.hist(cell_info["pct_counts_mt"], bins=50)
plt.xlabel("Mitochondrial gene percentage")
plt.ylabel("Number of cells")

plt.tight_layout()
plt.show()
```

## 🔬 ステップ8: 高度な解析（オプション）

### 細胞型予測

```bash
# 事前学習モデルを使用した細胞型予測
celline run predict_celltype
```

### 次元削減とクラスタリング

```bash
# 次元削減（PCA、UMAP）
celline run reduce

# バッチエフェクト補正
celline run integrate
```

### カスタム解析

```python
# Pythonでのカスタム解析例
from celline import Project
from celline.data import Seurat

# プロジェクトを読み込み
project = Project("./")

# Seuratオブジェクトを取得
seurat = project.seurat("GSE115189", "GSM3169075")

# カスタム解析を実行
# ...
```

## 📊 結果の解釈

### 品質管理指標

- **n_genes_by_counts**: 200-5000の範囲が一般的
- **pct_counts_mt**: 5%以下が推奨
- **doublet_score**: 0.3以下が一般的

### 細胞フィルタリング

フィルタリング後の統計：
```bash
# フィルタリング統計を確認
celline run info
```

## 🔄 ワークフローの自動化

### バッチ処理スクリプト

```bash
#!/bin/bash
# complete_workflow.sh

# プロジェクト初期化
celline init

# サンプル追加
celline run add GSE115189

# データ処理パイプライン
celline run download --nthread 4
celline run count
celline run create_seurat
celline run preprocess

# 解析
celline run predict_celltype
celline run reduce

echo "解析完了！celline interactiveで結果を確認してください。"
```

### Python API使用例

```python
from celline import Project
from celline.functions import Add, Download, Count, Preprocess

# プロジェクト作成
project = Project("./tutorial-project")

# 解析パイプライン
project.call(Add([Add.SampleInfo(id="GSE115189")]))
project.call(Download())
project.call(Count())
project.call(Preprocess())

print("解析完了！")
```

## 🎉 おめでとうございます！

このチュートリアルでは以下を完了しました：

- ✅ Cellineプロジェクトの作成
- ✅ 公的データベースからのデータ取得
- ✅ データの前処理と品質管理
- ✅ インタラクティブな解析環境の使用

## 📚 次のステップ

- [CLI Reference](/celline/cli) - より詳細なコマンドオプション
- [Interactive Mode](/celline/interactive) - Webインターフェースの詳細機能
- [API Reference](/celline/api) - Pythonプログラミングでの使用
- [Developer Guide](/celline/development) - カスタム機能の開発

## 💡 ヒントとベストプラクティス

### メモリ管理

```bash
# 大きなデータセットの場合
celline config --nthread 1  # スレッド数を削減
```

### 進行状況の監視

```bash
# ログファイルで進行状況を確認
tail -f resources/*/log/*.log
```

### データのバックアップ

```bash
# 重要な結果をバックアップ
tar -czf results_backup.tar.gz data/ results/
```

---

::alert{type="success"}
チュートリアル完了です！更なる解析のために[Advanced Usage](/celline/cli/advanced)をご覧ください。
::