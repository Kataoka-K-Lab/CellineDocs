# Configuration

Cellineの設定システムについて詳しく説明します。
プロジェクト固有の設定から、グローバル実行環境まで、効率的な解析のための設定方法を学びます。

## 🎯 設定システムの概要

Cellineは以下の階層で設定を管理します：

1. **プロジェクト設定** (`setting.toml`) - プロジェクト固有の設定
2. **サンプル設定** (`samples.toml`) - 解析対象サンプルの管理
3. **実行時設定** - コマンドライン引数による一時的な設定

## 📁 プロジェクト設定 (setting.toml)

プロジェクトのルートディレクトリに作成される設定ファイルです。

### 基本構造

```toml
[project]
name = "my-scrna-project"
version = "1.0.0"
description = "Single cell RNA-seq analysis project"

[execution]
system = "multithreading"
nthread = 4
pbs_server = ""

[R]
r_path = "/usr/bin/R"

[fetch]
wait_time = 4

[analysis]
target_species = "homo_sapiens"
reference_genome = "GRCh38"
```

### 設定セクション詳細

#### `[project]` - プロジェクト情報

| パラメータ | 型 | デフォルト | 説明 |
|-----------|---|-----------|------|
| `name` | string | directory名 | プロジェクト名 |
| `version` | string | "0.01" | プロジェクトバージョン |
| `description` | string | "" | プロジェクトの説明 |

#### `[execution]` - 実行環境

| パラメータ | 型 | デフォルト | 説明 |
|-----------|---|-----------|------|
| `system` | string | "multithreading" | 実行システム (`multithreading` / `PBS`) |
| `nthread` | integer | 1 | 並列実行スレッド数 |
| `pbs_server` | string | "" | PBSクラスターサーバー名 |

#### `[R]` - R環境設定

| パラメータ | 型 | デフォルト | 説明 |
|-----------|---|-----------|------|
| `r_path` | string | auto-detect | Rの実行パス |

#### `[fetch]` - データ取得設定

| パラメータ | 型 | デフォルト | 説明 |
|-----------|---|-----------|------|
| `wait_time` | integer | 4 | API呼び出し間のウェイト時間（秒） |

#### `[analysis]` - 解析設定

| パラメータ | 型 | デフォルト | 説明 |
|-----------|---|-----------|------|
| `target_species` | string | "" | 対象種（homo_sapiens, mus_musculus等） |
| `reference_genome` | string | "" | 参照ゲノム（GRCh38, GRCm39等） |

## 🧬 サンプル設定 (samples.toml)

解析対象のサンプル情報を管理するファイルです。通常は `celline run add` コマンドで自動生成されます。

### 基本構造

```toml
# シンプルな形式
GSM1234567 = "Control sample 1"
GSM1234568 = "Treatment sample 1"

# 詳細情報付きの形式
[GSM1234569]
title = "Control sample 2"
condition = "control"
replicate = 2
tissue = "brain"
cell_type = "mixed"

[GSM1234570]
title = "Treatment sample 2"
condition = "treatment"
replicate = 2
tissue = "brain"
cell_type = "mixed"
```

### 手動でのサンプル追加

```toml
# samples.tomlに直接記述
GSM5555555 = "Custom sample"

# または詳細情報付き
[GSM6666666]
title = "Custom detailed sample"
condition = "experimental"
batch = "batch1"
notes = "Special processing required"
```

## ⚙️ CLI設定コマンド

### 基本的な設定

```bash
# インタラクティブ設定
celline config

# 現在の設定を表示
celline config --show
```

### 実行システムの設定

```bash
# マルチスレッド実行を設定
celline config --system multithreading --nthread 8

# PBSクラスター実行を設定
celline config --system PBS --pbs-server my-cluster --nthread 16
```

### R環境の設定

```bash
# Rパスを手動設定
celline config --r-path /opt/R/4.3.0/bin/R

# 自動検出を使用
celline config --r-path auto
```

## 🔧 実行環境の詳細設定

### マルチスレッド実行

ローカルマシンでの並列実行に適しています。

```toml
[execution]
system = "multithreading"
nthread = 4  # CPUコア数に応じて調整

# 推奨設定例
# CPU 4コア: nthread = 2-3
# CPU 8コア: nthread = 4-6
# CPU 16コア: nthread = 8-12
```

**適用される処理:**
- データダウンロード
- ファイル処理
- 並列化可能な解析処理

### PBSクラスター実行

HPC環境での大規模解析に適しています。

```toml
[execution]
system = "PBS"
nthread = 16
pbs_server = "your-cluster-name"

[pbs]
queue = "normal"
walltime = "24:00:00"
memory = "64GB"
```

**PBSジョブテンプレート例:**

```bash
#!/bin/bash
#PBS -l select=1:ncpus=16:mem=64GB
#PBS -l walltime=24:00:00
#PBS -q normal
#PBS -N celline-job

cd $PBS_O_WORKDIR
module load R/4.3.0
celline run count
```

## 🌍 環境変数設定

Cellineは以下の環境変数を認識します：

```bash
# R環境
export R_HOME=/opt/R/4.3.0
export R_LIBS_USER=/home/user/R/library

# Cell Ranger
export CELLRANGER_PATH=/opt/cellranger-7.0.0

# 一時ディレクトリ
export TMPDIR=/scratch/tmp

# メモリ制限
export CELLINE_MAX_MEMORY=32G
```

### .bashrc/.zshrcでの設定

```bash
# ~/.bashrc または ~/.zshrc に追加
export CELLINE_CONFIG_DIR=$HOME/.config/celline
export CELLINE_CACHE_DIR=$HOME/.cache/celline
export CELLINE_R_PATH=/opt/R/bin/R
```

## 📊 プロファイル管理

複数の設定プロファイルを管理できます。

### プロファイルの作成

```bash
# 開発用プロファイル
mkdir -p ~/.config/celline/profiles/development
cat > ~/.config/celline/profiles/development/setting.toml << EOF
[execution]
system = "multithreading"
nthread = 2

[analysis]
debug_mode = true
verbose = true
EOF

# 本番用プロファイル  
mkdir -p ~/.config/celline/profiles/production
cat > ~/.config/celline/profiles/production/setting.toml << EOF
[execution]
system = "PBS"
nthread = 32
pbs_server = "production-cluster"

[analysis]
debug_mode = false
verbose = false
EOF
```

### プロファイルの使用

```bash
# プロファイルを指定して実行
celline --profile development run preprocess
celline --profile production run count
```

## 🔒 セキュリティ設定

### アクセス制御

```toml
[security]
# APIアクセス制限
api_allowed_hosts = ["localhost", "127.0.0.1"]
api_port = 8000

# データディレクトリの権限
data_permissions = "750"
result_permissions = "755"
```

### 認証設定

```toml
[auth]
# 将来的な機能
enable_auth = false
auth_provider = "none"  # "ldap", "oauth", etc.
```

## 📈 パフォーマンス設定

### メモリ管理

```toml
[performance]
# メモリ使用量の制限
max_memory_gb = 32
temp_dir = "/tmp/celline"

# キャッシュ設定
enable_cache = true
cache_size_gb = 10
cache_dir = "~/.cache/celline"
```

### 並列処理の調整

```toml
[parallel]
# I/O集約的タスク
io_workers = 4

# CPU集約的タスク  
cpu_workers = 8

# メモリ集約的タスク
memory_workers = 2
```

## 🔍 デバッグ設定

### ログレベル

```toml
[logging]
level = "INFO"  # DEBUG, INFO, WARNING, ERROR
log_file = "celline.log"
max_log_size_mb = 100
backup_count = 5

# モジュール別ログレベル
[logging.modules]
"celline.functions" = "DEBUG"
"celline.database" = "INFO"
"celline.api" = "WARNING"
```

### デバッグモード

```bash
# デバッグ情報付きで実行
celline --debug run preprocess

# 詳細ログ出力
celline --verbose run download
```

## 🔄 設定の継承と優先順位

設定は以下の優先順位で適用されます（高い順）：

1. **コマンドライン引数**
2. **環境変数**
3. **プロジェクト設定** (`setting.toml`)
4. **ユーザー設定** (`~/.config/celline/config.toml`)
5. **システム設定** (`/etc/celline/config.toml`)
6. **デフォルト値**

### 設定の確認

```bash
# 現在の有効な設定を表示
celline config --show-effective

# 設定ファイルの場所を表示
celline config --show-files

# 設定の詳細な解析
celline config --debug
```

## 🚨 設定のバリデーション

### 設定チェック

```bash
# 設定ファイルの構文チェック
celline config --validate

# 実行環境の検証
celline config --test-environment

# 依存関係の確認
celline config --check-dependencies
```

### 自動修復

```bash
# 設定の自動修復
celline config --repair

# デフォルト設定の復元
celline config --reset
```

## 📦 設定のエクスポート・インポート

### 設定のバックアップ

```bash
# 現在の設定をエクスポート
celline config --export > my-config.toml

# 設定を別のプロジェクトにインポート
celline config --import my-config.toml
```

### 設定テンプレート

```bash
# 設定テンプレートを生成
celline config --generate-template > template.toml

# テンプレートから新しいプロジェクトを作成
celline init --template template.toml new-project
```

---

::alert{type="info"}
より高度な設定については、[Advanced Usage](/celline/cli/advanced) を参照してください。
::