# Advanced Usage

Cellineの高度な機能と使用方法について説明します。
カスタムワークフロー、大規模データ処理、統合解析など、上級者向けの機能を詳しく解説します。

## 🚀 大規模データ処理

### クラスター環境での実行

#### PBS/Torqueクラスター

```bash
# PBSシステムの設定
celline config --system PBS --pbs-server my-cluster

# ジョブリソースの設定
cat > pbs_config.toml << EOF
[pbs]
queue = "normal"
walltime = "48:00:00"
memory = "128GB"
ncpus = 32
nodes = 1
EOF

# 大規模データセットの処理
celline run count --config pbs_config.toml
```

#### Slurmクラスター（カスタム実装）

```bash
#!/bin/bash
#SBATCH --job-name=celline-analysis
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16

module load R/4.3.0 cellranger/7.0.0

# 並列処理で複数サンプルを処理
celline run download --nthread 16
celline run count --nthread 16
```

### メモリ効率的な処理

```python
# 大きなデータセットでのメモリ管理
from celline import Project
from celline.utils.memory import MemoryManager

project = Project("./large-dataset")

# メモリ使用量を制限
with MemoryManager(max_memory="32GB") as mem:
    # チャンクごとに処理
    for chunk in project.iter_samples(chunk_size=10):
        project.process_chunk(chunk)
        mem.clear_cache()
```

## 🔬 カスタムワークフロー

### パイプライン設計

```python
from celline import Project, Pipeline
from celline.functions import *

# カスタムパイプラインの定義
class CustomSingleCellPipeline(Pipeline):
    def __init__(self, project: Project):
        super().__init__(project)
        
    def run(self):
        # 段階的な処理
        self.stage1_data_acquisition()
        self.stage2_quality_control()
        self.stage3_analysis()
        self.stage4_integration()
    
    def stage1_data_acquisition(self):
        """データ取得ステージ"""
        self.project.call(Add(self.sample_list))
        self.project.call(Download())
        
    def stage2_quality_control(self):
        """品質管理ステージ"""
        self.project.call(Count())
        self.project.call(Preprocess())
        
    def stage3_analysis(self):
        """解析ステージ"""
        self.project.call(PredictCelltype())
        self.project.call(Reduce())
        
    def stage4_integration(self):
        """統合ステージ"""
        self.project.call(Integrate())

# パイプラインの実行
project = Project("./my-analysis")
pipeline = CustomSingleCellPipeline(project)
pipeline.run()
```

### 条件分岐処理

```python
from celline.functions import *

project = Project("./conditional-analysis")

# 条件に基づく処理分岐
project.call_if_else(
    condition=lambda p: p.has_spatial_data(),
    true=SpatialAnalysis(),
    false=StandardAnalysis()
)

# データ型に基づく処理
if project.is_10x_data():
    project.call(Count())
elif project.is_smartseq_data():
    project.call(SmartSeqProcessing())
else:
    project.call(CustomProcessing())
```

## 📊 バッチ処理と自動化

### 複数プロジェクトの一括処理

```bash
#!/bin/bash
# batch_analysis.sh

# プロジェクトリストファイル
PROJECTS="project_list.txt"

while IFS=',' read -r project_name gse_id description; do
    echo "Processing $project_name ($gse_id)"
    
    # プロジェクトディレクトリ作成
    mkdir -p "$project_name"
    cd "$project_name"
    
    # 解析パイプライン実行
    celline init "$project_name"
    celline run add "$gse_id"
    celline run download --nthread 8
    celline run count
    celline run preprocess
    celline run predict_celltype
    
    # 結果をアーカイブ
    tar -czf "../${project_name}_results.tar.gz" results/
    
    cd ..
done < "$PROJECTS"
```

### Makefileを使用した自動化

```makefile
# Makefile for Celline analysis

# 変数定義
PROJECT_NAME := my-scrna-analysis
GSE_ID := GSE123456
NTHREAD := 8

# デフォルトターゲット
all: init download count preprocess analyze

# プロジェクト初期化
init:
	celline init $(PROJECT_NAME)

# サンプル追加
add:
	celline run add $(GSE_ID)

# データダウンロード
download: add
	celline run download --nthread $(NTHREAD)

# カウント処理
count: download
	celline run count

# 前処理
preprocess: count
	celline run preprocess

# 解析
analyze: preprocess
	celline run predict_celltype
	celline run reduce
	celline run integrate

# 結果の確認
check:
	celline info
	find results/ -name "*.csv" -o -name "*.png" | head -10

# クリーンアップ
clean:
	rm -rf resources/*/raw resources/*/tmp
	
# 完全クリーンアップ
distclean:
	rm -rf resources/ data/ results/

.PHONY: all init add download count preprocess analyze check clean distclean
```

## 🔄 統合解析

### 複数データセットの統合

```python
from celline.functions.integrate import MultiDatasetIntegration

# 複数データセットの統合解析
integration = MultiDatasetIntegration([
    "dataset1_GSE123456",
    "dataset2_GSE789012", 
    "dataset3_GSE345678"
])

project = Project("./integrated-analysis")
project.call(integration)

# バッチエフェクト補正手法の比較
methods = ["harmony", "combat", "mnn", "cca"]
for method in methods:
    integration.method = method
    project.call(integration)
    project.save_results(f"integration_{method}")
```

### メタ解析

```python
from celline.functions.meta import MetaAnalysis

# メタ解析の実行
meta = MetaAnalysis([
    {"gse": "GSE123456", "condition": "control", "tissue": "brain"},
    {"gse": "GSE789012", "condition": "disease", "tissue": "brain"},
    {"gse": "GSE345678", "condition": "treatment", "tissue": "brain"}
])

project = Project("./meta-analysis")
project.call(meta)

# 結果の統計解析
meta.perform_differential_analysis()
meta.generate_forest_plots()
```

## 🧬 カスタム解析機能

### 新しい解析関数の作成

```python
from celline.functions._base import CellineFunction
import argparse

class CustomCellTypeClassification(CellineFunction):
    """カスタム細胞型分類関数"""
    
    def __init__(self, model_path: str, confidence_threshold: float = 0.8):
        super().__init__()
        self.model_path = model_path
        self.confidence_threshold = confidence_threshold
    
    def call(self, project):
        """メイン処理"""
        import pickle
        import scanpy as sc
        
        # カスタムモデルの読み込み
        with open(self.model_path, 'rb') as f:
            model = pickle.load(f)
        
        # 各サンプルの処理
        for sample_id in project.get_samples():
            adata = project.load_sample_data(sample_id)
            
            # 特徴量の抽出
            features = self.extract_features(adata)
            
            # 予測の実行
            predictions = model.predict(features)
            confidence = model.predict_proba(features)
            
            # 信頼度フィルタリング
            high_conf_mask = confidence.max(axis=1) >= self.confidence_threshold
            
            # 結果の保存
            adata.obs['custom_celltype'] = predictions
            adata.obs['celltype_confidence'] = confidence.max(axis=1)
            adata.obs['high_confidence'] = high_conf_mask
            
            project.save_sample_data(sample_id, adata)
        
        return project
    
    def extract_features(self, adata):
        """特徴量抽出"""
        # カスタム特徴量の計算
        import numpy as np
        
        # 高可変遺伝子の発現
        hvg_expr = adata[:, adata.var.highly_variable].X.toarray()
        
        # パスウェイスコア
        pathway_scores = self.calculate_pathway_scores(adata)
        
        # 結合
        features = np.concatenate([hvg_expr, pathway_scores], axis=1)
        return features
    
    def add_cli_args(self, parser: argparse.ArgumentParser):
        """CLI引数の追加"""
        parser.add_argument('--model-path', required=True,
                          help='Path to the trained classification model')
        parser.add_argument('--confidence-threshold', type=float, default=0.8,
                          help='Confidence threshold for predictions')
    
    def get_description(self):
        return "Custom cell type classification using trained models"
```

### 解析関数の登録

```python
# カスタム関数をCellineに登録
from celline.cli.registry import get_registry

registry = get_registry()
registry.register_function(
    name="custom_classify",
    class_ref=CustomCellTypeClassification,
    module_path="my_custom_functions.classification"
)

# CLI経由で使用
# celline run custom_classify --model-path my_model.pkl
```

## 📊 データ可視化の高度な機能

### カスタムプロット作成

```python
from celline.visualization import CellinePlotter

class AdvancedVisualization(CellineFunction):
    def call(self, project):
        plotter = CellinePlotter(project)
        
        # 複数サンプルの統合UMAP
        plotter.integrated_umap(
            samples=project.get_samples(),
            color_by="celltype",
            split_by="condition"
        )
        
        # 発現ヒートマップ
        plotter.expression_heatmap(
            genes=["CD4", "CD8A", "IL2", "IFNG"],
            group_by="celltype",
            save_path="results/expression_heatmap.pdf"
        )
        
        # 軌跡解析
        plotter.trajectory_plot(
            root_cell="stem_cell",
            save_path="results/trajectory.pdf"
        )
        
        return project
```

### インタラクティブ可視化

```python
import plotly.graph_objects as go
from celline.interactive import InteractivePlotter

def create_interactive_plots(project):
    plotter = InteractivePlotter(project)
    
    # 3D UMAP
    fig = plotter.plot_3d_umap(
        color_by="celltype",
        hover_data=["sample", "condition", "batch"]
    )
    
    # インタラクティブ遺伝子発現
    expr_fig = plotter.plot_gene_expression(
        genes=["CD4", "CD8A"],
        plot_type="violin"
    )
    
    # ダッシュボードの作成
    dashboard = plotter.create_dashboard([fig, expr_fig])
    dashboard.serve(port=8050)
```

## 🔌 API統合

### 外部データベースとの連携

```python
from celline.database import ExternalDBConnector

class CustomDBIntegration(CellineFunction):
    def __init__(self, db_config: dict):
        self.db_config = db_config
    
    def call(self, project):
        # 外部データベースから追加情報を取得
        connector = ExternalDBConnector(self.db_config)
        
        for sample_id in project.get_samples():
            # メタデータの拡張
            extended_meta = connector.fetch_extended_metadata(sample_id)
            
            # 遺伝子アノテーションの更新
            gene_annotations = connector.fetch_gene_annotations(
                species=project.get_species()
            )
            
            # プロジェクトに統合
            project.update_sample_metadata(sample_id, extended_meta)
            project.update_gene_annotations(gene_annotations)
        
        return project
```

### 機械学習パイプライン統合

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier

class MLPipeline(CellineFunction):
    def __init__(self, model_config: dict):
        self.model_config = model_config
    
    def call(self, project):
        # 機械学習パイプラインの構築
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', RandomForestClassifier(**self.model_config))
        ])
        
        # 訓練データの準備
        X_train, y_train = project.prepare_training_data()
        
        # モデルの訓練
        pipeline.fit(X_train, y_train)
        
        # 予測の実行
        for sample_id in project.get_samples():
            X_test = project.prepare_test_data(sample_id)
            predictions = pipeline.predict(X_test)
            project.save_predictions(sample_id, predictions)
        
        # モデルの保存
        project.save_model(pipeline, "trained_classifier.pkl")
        
        return project
```

## 🔧 パフォーマンス最適化

### 並列処理の最適化

```python
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from celline.utils import ProcessManager

class OptimizedProcessing(CellineFunction):
    def call(self, project):
        samples = project.get_samples()
        
        # CPU集約的タスクはプロセス並列
        with ProcessPoolExecutor(max_workers=8) as executor:
            cpu_futures = [
                executor.submit(self.cpu_intensive_task, sample)
                for sample in samples
            ]
        
        # I/O集約的タスクはスレッド並列
        with ThreadPoolExecutor(max_workers=16) as executor:
            io_futures = [
                executor.submit(self.io_intensive_task, sample)
                for sample in samples
            ]
        
        # 結果の収集
        cpu_results = [f.result() for f in cpu_futures]
        io_results = [f.result() for f in io_futures]
        
        return project
```

### メモリ使用量の最適化

```python
from celline.utils.memory import MemoryProfiler, LazyLoader

class MemoryOptimizedAnalysis(CellineFunction):
    def call(self, project):
        with MemoryProfiler() as profiler:
            # 遅延読み込みを使用
            lazy_data = LazyLoader(project.data_path)
            
            # チャンクごとの処理
            for chunk in lazy_data.iter_chunks(chunk_size=1000):
                self.process_chunk(chunk)
                
                # メモリ使用量の監視
                if profiler.get_memory_usage() > profiler.memory_limit:
                    profiler.clear_cache()
                    gc.collect()
        
        return project
```

## 🚨 エラーハンドリングと復旧

### 堅牢なエラーハンドリング

```python
from celline.utils.exceptions import CellineException
from celline.utils.recovery import RecoveryManager

class RobustAnalysis(CellineFunction):
    def call(self, project):
        recovery = RecoveryManager(project)
        
        try:
            # 分析ステップの実行
            for step in self.analysis_steps:
                checkpoint = recovery.create_checkpoint(step.name)
                
                try:
                    step.execute(project)
                    recovery.mark_success(checkpoint)
                except Exception as e:
                    recovery.mark_failure(checkpoint, e)
                    
                    # 自動復旧の試行
                    if recovery.can_recover(step):
                        recovery.recover_from_checkpoint(checkpoint)
                    else:
                        raise CellineException(f"Step {step.name} failed: {e}")
        
        except CellineException as e:
            # 詳細なエラーレポートの生成
            recovery.generate_error_report(e)
            raise
        
        return project
```

---

::alert{type="warning"}
高度な機能を使用する際は、十分なテストを行い、本番環境での実行前に小規模なデータセットで動作確認を行ってください。
::