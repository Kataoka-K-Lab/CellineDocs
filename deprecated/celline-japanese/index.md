# Functions Reference

Cellineの各解析関数の詳細なリファレンスです。
すべての関数のパラメータ、使用例、実装詳細について説明します。

## 📚 関数分類

### 🗄️ データ管理
- [Add](#add) - サンプル追加とメタデータ取得
- [Download](#download) - シーケンシングデータのダウンロード
- [SyncDB](#syncdb) - データベース同期

### 🧮 データ処理
- [Count](#count) - Cell Rangerによるカウント処理
- [CreateSeuratObject](#createseuratobject) - Seuratオブジェクト作成
- [SetTranscriptome](#settranscriptome) - 転写情報データベース設定

### 🔬 品質管理・前処理
- [Preprocess](#preprocess) - 品質管理とフィルタリング
- [VCount](#vcount) - カウントデータの検証

### 📊 解析・可視化
- [PredictCelltype](#predictcelltype) - 細胞型予測
- [Reduce](#reduce) - 次元削減
- [Integrate](#integrate) - データ統合とバッチエフェクト補正
- [BatchCor](#batchcor) - バッチ相関解析

### 🔧 ユーティリティ
- [Info](#info) - システム情報表示
- [Initialize](#initialize) - 初期化とセットアップ
- [Interactive](#interactive) - インタラクティブモード
- [Job](#job) - ジョブ管理
- [Bash](#bash) - シェルコマンド実行

---

## 📊 データ管理関数

### Add

サンプルIDを公的データベースから取得し、プロジェクトに追加します。

#### クラス定義

```python
class Add(CellineFunction):
    def __init__(self, sample_id: Union[List[SampleInfo], pl.DataFrame]) -> None
```

#### パラメータ

**SampleInfo**
```python
class SampleInfo(NamedTuple):
    id: str                    # サンプルID (GSE, GSM, SRR等)
    title: Optional[str] = ""  # サンプルタイトル
```

#### 使用例

```python
from celline.functions.add import Add
import polars as pl

# SampleInfoリストを使用
samples = [
    Add.SampleInfo(id="GSE123456", title="Brain tissue study"),
    Add.SampleInfo(id="GSM789012", title="Control sample"),
    Add.SampleInfo(id="GSM789013", title="Treatment sample")
]
add_func = Add(samples)
project.call(add_func)

# DataFrameを使用
df = pl.DataFrame({
    "id": ["GSE123456", "GSM789012"],
    "title": ["Study 1", "Sample 1"]
})
add_func = Add(df)
project.call(add_func)

# CLIからの使用
# celline run add GSE123456
# celline run add GSM789012 GSM789013 --title "My samples"
# celline run add --from-file samples.csv
```

#### メソッド

##### `get_samples() -> Dict[str, str]`

現在のプロジェクトのサンプル情報を取得します。

```python
add_func = Add([])
samples = add_func.get_samples()
print(samples)  # {'GSM123456': 'Sample 1', 'GSM123457': 'Sample 2'}
```

#### サポートするデータベース

- **GEO (Gene Expression Omnibus)**
  - GSE (Studies)
  - GSM (Samples)
- **SRA (Sequence Read Archive)**
  - SRR (Runs)
- **CNCB (China National Center for Bioinformation)**
  - PRJCA (Projects)
  - CRA (Study)
  - CRR (Runs)

#### エラーハンドリング

```python
try:
    add_func = Add([Add.SampleInfo(id="INVALID_ID")])
    project.call(add_func)
except ValueError as e:
    print(f"Invalid sample ID format: {e}")
except ConnectionError as e:
    print(f"Database connection failed: {e}")
```

---

### Download

追加されたサンプルのシーケンシングデータをダウンロードします。

#### クラス定義

```python
class Download(CellineFunction):
    def __init__(
        self,
        then: Optional[Callable[[str], None]] = None,
        catch: Optional[Callable[[subprocess.CalledProcessError], None]] = None,
    ) -> None
```

#### パラメータ

- `then`: ダウンロード完了時のコールバック関数
- `catch`: エラー発生時のコールバック関数

#### 使用例

```python
from celline.functions.download import Download

# 基本的なダウンロード
download = Download()
project.call(download)

# コールバック付きダウンロード
def on_complete(sample_id):
    print(f"✓ Downloaded: {sample_id}")

def on_error(error):
    print(f"✗ Error downloading: {error}")

download = Download(then=on_complete, catch=on_error)
project.call(download)

# CLIからの使用
# celline run download
# celline run download --nthread 4
# celline run download --force
```

#### 内部データ構造

##### `JobContainer`

```python
class JobContainer(NamedTuple):
    filetype: str           # データ形式 (FASTQ, SRA, etc.)
    nthread: str           # スレッド数
    cluster_server: str    # クラスターサーバー名
    jobname: str          # ジョブ名
    logpath: str          # ログファイルパス
    sample_id: str        # サンプルID
    download_target: str  # ダウンロード先ディレクトリ
    download_source: str  # ダウンロード元URL
    run_ids_str: str      # ランIDの文字列
```

#### ダウンロード対象

- **FASTQ files** - Illumina/10x Genomicsデータ
- **SRA files** - SRAアーカイブからのダウンロード
- **BAM files** - アライメント済みデータ

#### 進行状況の監視

```python
# ログファイルでの進行状況確認
import os
from celline.utils.path import Path

def monitor_download_progress(sample_id):
    path = Path("GSE123456", sample_id)
    log_files = os.listdir(path.resources_sample_log)
    
    for log_file in log_files:
        if "download" in log_file:
            with open(f"{path.resources_sample_log}/{log_file}") as f:
                print(f.read())
```

---

### SyncDB

ローカルデータベースを最新の公的データベースと同期します。

#### クラス定義

```python
class SyncDB(CellineFunction):
    def __init__(self) -> None
```

#### 使用例

```python
from celline.functions.sync_DB import SyncDB

# データベース同期
sync = SyncDB()
project.call(sync)

# CLIからの使用
# celline run sync_DB
```

#### 同期対象データベース

- **GEO metadata** - GSE/GSM情報
- **SRA metadata** - SRR情報  
- **Transcriptome references** - 参照転写情報
- **Gene annotations** - 遺伝子アノテーション

---

## 🧮 データ処理関数

### Count

Cell Rangerを使用してFASTQファイルから発現カウントマトリクスを生成します。

#### クラス定義

```python
class Count(CellineFunction):
    def __init__(self) -> None
```

#### 使用例

```python
from celline.functions.count import Count

# カウント処理
count = Count()
project.call(count)

# CLIからの使用
# celline run count
```

#### 出力ファイル

- `filtered_feature_bc_matrix.h5` - フィルタ済み発現マトリクス
- `raw_feature_bc_matrix.h5` - 生発現マトリクス
- `metrics_summary.csv` - 品質メトリクス
- `web_summary.html` - HTMLレポート

#### Cell Ranger設定

```python
# Cell Rangerパラメータの確認
import subprocess

# Cell Rangerのバージョン確認
result = subprocess.run(['cellranger', '--version'], 
                       capture_output=True, text=True)
print(f"Cell Ranger version: {result.stdout}")

# 利用可能なコマンド
commands = ['count', 'aggr', 'reanalyze', 'mkfastq']
```

---

### CreateSeuratObject

カウントデータからSeuratオブジェクトを作成し、基本的な解析を実行します。

#### クラス定義

```python
class CreateSeuratObject(CellineFunction):
    def __init__(self) -> None
```

#### 使用例

```python
from celline.functions.create_seurat import CreateSeuratObject

# Seuratオブジェクト作成
create_seurat = CreateSeuratObject()
project.call(create_seurat)

# CLIからの使用
# celline run create_seurat
```

#### Rコード例（内部実行）

```r
# Celline内部で実行されるRコード
library(Seurat)
library(SeuratDisk)

# データ読み込み
raw <- Read10X_h5(h5_path)

# Seuratオブジェクト作成
seurat_obj <- CreateSeuratObject(raw, project = proj) %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(features = VariableFeatures(.)) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(dims = 1:20) %>%
    RunUMAP(dims = 1:20)

# 保存
saveRDS(seurat_obj, h5seurat_path)
```

#### 生成される解析結果

- **正規化**: LogNormalize
- **高可変遺伝子**: FindVariableFeatures (top 2000)
- **スケーリング**: ScaleData
- **主成分分析**: PCA (50 components)
- **近傍グラフ**: FindNeighbors
- **クラスタリング**: Leiden algorithm
- **UMAP**: 2次元埋め込み

---

### SetTranscriptome

解析で使用する転写情報データベースを設定します。

#### クラス定義

```python
class SetTranscriptome(CellineFunction):
    def __init__(self, species: Optional[str] = None, version: Optional[str] = None) -> None
```

#### パラメータ

- `species`: 対象種（homo_sapiens, mus_musculus等）
- `version`: ゲノムバージョン（GRCh38, GRCm39等）

#### 使用例

```python
from celline.functions.set_transcriptome import SetTranscriptome

# 自動検出
set_ref = SetTranscriptome()
project.call(set_ref)

# 手動指定
set_ref = SetTranscriptome(species="homo_sapiens", version="GRCh38")
project.call(set_ref)

# CLIからの使用
# celline run set_transcriptome
# celline run set_transcriptome --species homo_sapiens --version GRCh38
```

#### サポート種・バージョン

| 種 | 利用可能バージョン |
|----|------------------|
| homo_sapiens | GRCh38, GRCh37 |
| mus_musculus | GRCm39, GRCm38 |
| rattus_norvegicus | Rnor6.0 |
| danio_rerio | GRCz11 |

---

## 🔬 品質管理・前処理関数

### Preprocess

品質管理、フィルタリング、二重細胞検出を実行します。

#### クラス定義

```python
class Preprocess(CellineFunction):
    def __init__(self, target_celltype: Optional[List[str]] = None) -> None
```

#### パラメータ

- `target_celltype`: フィルタリング対象の細胞型リスト

#### 使用例

```python
from celline.functions.preprocess import Preprocess

# 基本的な前処理
preprocess = Preprocess()
project.call(preprocess)

# 特定細胞型のみを対象
preprocess = Preprocess(target_celltype=["Neuron", "Astrocyte", "Oligodendrocyte"])
project.call(preprocess)

# CLIからの使用
# celline run preprocess
# celline run preprocess --target-celltype Neuron Astrocyte
```

#### QCメトリクス

```python
# 生成されるQCメトリクス
qc_metrics = {
    "n_genes_by_counts": "検出遺伝子数",
    "total_counts": "総リード数", 
    "pct_counts_mt": "ミトコンドリア遺伝子比率",
    "doublet_score": "二重細胞スコア",
    "predicted_doublets": "二重細胞予測"
}
```

#### フィルタリング基準

| メトリクス | 下限 | 上限 |
|-----------|------|------|
| n_genes_by_counts | 200 | 5000 |
| pct_counts_mt | - | 20% |
| doublet_score | - | 0.3 |

#### 出力ファイル

- `cell_info.tsv` - 細胞レベルのQC情報
- `qc_plots.pdf` - QCプロット
- `filtering_report.txt` - フィルタリングレポート

---

### VCount

カウントデータの妥当性を検証します。

#### クラス定義

```python
class VCount(CellineFunction):
    def __init__(self) -> None
```

#### 使用例

```python
from celline.functions.vcount import VCount

# カウントデータ検証
vcount = VCount()
project.call(vcount)
```

#### 検証項目

- ファイルの存在確認
- マトリクスの次元チェック
- カウント分布の確認
- 品質メトリクスの妥当性

---

## 📊 解析・可視化関数

### PredictCelltype

機械学習モデルを使用して細胞型を予測します。

#### クラス定義

```python
class PredictCelltype(CellineFunction):
    def __init__(self, model_path: Optional[str] = None) -> None
```

#### パラメータ

- `model_path`: カスタム学習済みモデルのパス

#### 使用例

```python
from celline.functions.predict_celltype import PredictCelltype

# デフォルトモデルを使用
predict = PredictCelltype()
project.call(predict)

# カスタムモデルを使用
predict = PredictCelltype(model_path="/path/to/custom_model.pkl")
project.call(predict)

# CLIからの使用
# celline run predict_celltype
# celline run predict_celltype --model-path /path/to/model.pkl
```

#### モデル情報

##### 事前学習済みモデル

- **Human Brain**: 人脳細胞型分類（Allen Brain Atlas）
- **Mouse Brain**: マウス脳細胞型分類（Allen Brain Atlas）
- **PBMC**: 末梢血単核球細胞分類
- **Universal**: 汎用細胞型分類

##### カスタムモデル作成

```python
from celline.functions.predict_celltype import BuildCellTypeModel

# カスタムモデルの学習
model_builder = BuildCellTypeModel(
    training_data="/path/to/training.h5ad",
    features=["CD4", "CD8A", "CD19", "CD14"],
    model_type="random_forest"
)
project.call(model_builder)
```

---

### Reduce

PCA、UMAP、t-SNEによる次元削減を実行します。

#### クラス定義

```python
class Reduce(CellineFunction):
    def __init__(
        self, 
        methods: List[str] = ["pca", "umap"],
        n_components: int = 50,
        n_neighbors: int = 15
    ) -> None
```

#### パラメータ

- `methods`: 使用する次元削減手法
- `n_components`: 主成分数
- `n_neighbors`: 近傍点数（UMAP）

#### 使用例

```python
from celline.functions.reduce import Reduce

# デフォルト設定
reduce = Reduce()
project.call(reduce)

# カスタム設定
reduce = Reduce(
    methods=["pca", "umap", "tsne"],
    n_components=30,
    n_neighbors=20
)
project.call(reduce)

# CLIからの使用
# celline run reduce
# celline run reduce --methods pca umap tsne --n-components 30
```

#### 利用可能な手法

| 手法 | 説明 | パラメータ |
|-----|------|-----------|
| PCA | 主成分分析 | n_components |
| UMAP | Uniform Manifold Approximation | n_neighbors, min_dist |
| t-SNE | t-distributed Stochastic Neighbor Embedding | perplexity |

---

### Integrate

複数サンプルの統合とバッチエフェクト補正を実行します。

#### クラス定義

```python
class Integrate(CellineFunction):
    def __init__(
        self,
        method: str = "harmony",
        batch_key: str = "sample",
        n_components: int = 50
    ) -> None
```

#### パラメータ

- `method`: 統合手法（harmony, combat, mnn, cca）
- `batch_key`: バッチ情報のカラム名
- `n_components`: 統合に使用する主成分数

#### 使用例

```python
from celline.functions.integrate import Integrate

# Harmonyを使用
integrate = Integrate(method="harmony")
project.call(integrate)

# MNNを使用
integrate = Integrate(
    method="mnn",
    batch_key="dataset",
    n_components=30
)
project.call(integrate)

# CLIからの使用
# celline run integrate
# celline run integrate --method mnn --batch-key dataset
```

#### 統合手法の比較

| 手法 | 特徴 | 適用場面 |
|-----|------|----------|
| Harmony | 高速、効果的 | 大規模データセット |
| Combat | 従来的手法 | 小〜中規模データセット |
| MNN | 生物学的変動保持 | 異なるプロトコル |
| CCA | Seurat標準 | 10x Genomicsデータ |

---

### BatchCor

バッチ間の相関を解析し、バッチエフェクトを評価します。

#### クラス定義

```python
class BatchCor(CellineFunction):
    def __init__(self, batch_keys: List[str] = ["sample", "dataset"]) -> None
```

#### 使用例

```python
from celline.functions.batch_cor import BatchCor

# バッチ相関解析
batch_cor = BatchCor(batch_keys=["sample", "condition", "batch"])
project.call(batch_cor)
```

---

## 🔧 ユーティリティ関数

### Info

システム情報とプロジェクト状態を表示します。

#### クラス定義

```python
class Info(CellineFunction):
    def __init__(self) -> None
```

#### 使用例

```python
from celline.functions.info import Info

# システム情報表示
info = Info()
project.call(info)

# CLIからの使用
# celline info
```

### Initialize

プロジェクトの初期化とセットアップを実行します。

#### クラス定義

```python
class Initialize(CellineFunction):
    def __init__(self) -> None
```

#### 使用例

```python
from celline.functions.initialize import Initialize

# 初期化
init = Initialize()
project.call(init)

# CLIからの使用
# celline init
```

### Interactive

インタラクティブWebインターフェースを起動します。

#### クラス定義

```python
class Interactive(CellineFunction):
    def __init__(self, port: int = 8080) -> None
```

#### 使用例

```python
from celline.functions.interactive import Interactive

# デフォルトポートで起動
interactive = Interactive()
project.call(interactive)

# カスタムポートで起動
interactive = Interactive(port=8090)
project.call(interactive)

# CLIからの使用
# celline interactive
# celline run interactive --port 8090
```

### Job

ジョブの管理と監視を行います。

#### クラス定義

```python
class Job(CellineFunction):
    def __init__(self, action: str = "status") -> None
```

#### 使用例

```python
from celline.functions.job import Job

# ジョブステータス確認
job = Job(action="status")
project.call(job)

# ジョブキャンセル
job = Job(action="cancel")
project.call(job)
```

### Bash

シェルコマンドを実行します。

#### クラス定義

```python
class Bash(CellineFunction):
    def __init__(self, command: str, timeout: Optional[int] = None) -> None
```

#### 使用例

```python
from celline.functions.bash import Bash

# シェルコマンド実行
bash = Bash("ls -la resources/")
project.call(bash)

# タイムアウト付き実行
bash = Bash("long_running_command.sh", timeout=3600)
project.call(bash)
```

---

::alert{type="info"}
各関数の詳細な実装例とカスタマイズ方法については、[Developer Guide](/celline/development) を参照してください。
::