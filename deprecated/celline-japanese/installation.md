# Installation

Cellineのインストール方法について詳しく説明します。
環境に応じて最適なインストール方法を選択してください。

## 🐍 Python環境

### システム要件

| 要件 | 最小バージョン | 推奨バージョン |
|------|----------------|----------------|
| Python | 3.10 | 3.11+ |
| pip | 21.0 | 最新 |
| メモリ | 4GB | 8GB+ |
| ストレージ | 1GB | 10GB+ |

### Python仮想環境の作成

```bash
# venvを使用
python -m venv celline-env
source celline-env/bin/activate  # Linux/Mac
# または
celline-env\Scripts\activate     # Windows

# condaを使用
conda create -n celline python=3.11
conda activate celline
```

## 📦 インストール方法

### 方法1: pip経由（推奨）

```bash
# 最新安定版をインストール
pip install celline

# アップグレード
pip install --upgrade celline

# 特定バージョンを指定
pip install celline==0.1.10
```

### 方法2: UV経由

```bash
# UVがインストールされていない場合
curl -LsSf https://astral.sh/uv/install.sh | sh

# Cellineをプロジェクトに追加
uv add celline

# または、一時的に使用
uv run --with celline celline --help
```

### 方法3: 開発版インストール

```bash
# GitHubから最新の開発版をインストール
pip install git+https://github.com/YUYA556223/celline.git

# 特定のブランチを指定
pip install git+https://github.com/YUYA556223/celline.git@develop
```

## 🔧 依存関係

Cellineは以下の主要パッケージに依存しています：

### Python依存関係

```text
argparse>=1.4.0
continuousvi>=0.1.5
fastapi>=0.104.0
inquirer>=3.4.0
pandas>=2.2.3
polars>=1.26.0
pyarrow>=19.0.1
pydantic>=2.5.0
requests>=2.31.0
rich>=14.0.0
scanpy>=1.11.1
scrublet>=0.2.3
tqdm>=4.67.1
uvicorn>=0.24.0
```

### R依存関係

Cellineは一部の機能でRを使用します。以下のRパッケージが必要です：

```r
# 必須パッケージ
install.packages(c(
  "Seurat",
  "SeuratDisk", 
  "tidyverse",
  "hdf5r"
))

# オプションパッケージ
install.packages(c(
  "scater",
  "scran",
  "SingleCellExperiment",
  "BiocManager"
))
```

## 🚀 インストール後の確認

### 基本的な動作確認

```bash
# Cellineがインストールされているか確認
celline --help

# バージョン確認
python -c "import celline; print('Celline version:', celline.__version__)"

# 利用可能な関数を表示
celline list
```

### システム情報の確認

```bash
# システム情報を表示
celline info

# 設定状態を確認
celline config
```

## 🔧 外部ツールのインストール

### Cell Ranger（オプション）

10x Genomicsデータを処理する場合に必要：

```bash
# Cell Rangerのダウンロード（要アカウント登録）
# https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

# インストール後、パスを通す
export PATH=/path/to/cellranger:$PATH
```

### R環境の設定

```bash
# Rがインストールされているか確認
R --version

# Rのパスを確認
which R

# CellineでRパスを設定
celline config
# または
export R_HOME=/usr/lib/R
```

## 🐳 Docker環境での使用

### Docker image（将来実装予定）

```bash
# Docker imageの使用例（開発中）
docker run -it celline/celline:latest

# マウントしてローカルデータを使用
docker run -v $(pwd):/workspace celline/celline:latest
```

## 🌐 開発環境のセットアップ

### ソースからのインストール

```bash
# リポジトリをクローン
git clone https://github.com/YUYA556223/celline.git
cd celline

# 開発依存関係も含めてインストール
pip install -e ".[dev]"

# または、uvを使用
uv sync --all-extras
```

### 開発ツールの設定

```bash
# コードフォーマッタ
pip install ruff black

# テストツール
pip install pytest pytest-cov

# ドキュメント生成
pip install sphinx sphinx-rtd-theme
```

## 🔍 トラブルシューティング

### 一般的な問題と解決策

#### 1. Pythonバージョン関連

```bash
# Pythonバージョンが古い場合
python --version  # 3.10+であることを確認

# pyenvを使用してPythonバージョンを管理
pyenv install 3.11.0
pyenv global 3.11.0
```

#### 2. 依存関係の競合

```bash
# 依存関係の競合を解決
pip install --force-reinstall celline

# または、クリーンな環境で再インストール
pip uninstall celline
pip install celline
```

#### 3. R関連の問題

```bash
# Rが見つからない場合
sudo apt-get install r-base r-base-dev  # Ubuntu/Debian
brew install r                          # macOS

# Rパッケージのインストールに失敗する場合
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

#### 4. メモリ関連の問題

```bash
# メモリ使用量を削減
celline config --nthread 1

# スワップファイルを設定（Linux）
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

#### 5. 権限関連の問題

```bash
# ユーザーレベルでインストール
pip install --user celline

# 権限の問題を回避
pip install --no-deps celline
```

## ✅ インストール確認チェックリスト

- [ ] Python 3.10+ がインストールされている
- [ ] pip が最新バージョンである
- [ ] `celline --help` が正常に実行される
- [ ] `celline list` で関数一覧が表示される
- [ ] R環境が設定されている（Rを使用する場合）
- [ ] 必要なRパッケージがインストールされている
- [ ] 十分なディスク容量がある

## 🔄 アップグレード

### 定期的なアップグレード

```bash
# 最新バージョンにアップグレード
pip install --upgrade celline

# 特定バージョンにダウングレード
pip install celline==0.1.9

# 開発版にアップグレード
pip install --upgrade git+https://github.com/YUYA556223/celline.git
```

### 設定ファイルの移行

新しいバージョンでは設定ファイルの形式が変更される場合があります：

```bash
# 設定をバックアップ
cp setting.toml setting.toml.backup

# 新しい設定形式に移行（必要な場合）
celline config --migrate
```

---

> **注意**: インストールで問題が発生した場合は、[Troubleshooting](/celline/troubleshooting) セクションを参照してください。