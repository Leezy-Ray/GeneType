# GeneType 基因型分析

基于 Sanger 测序色谱（.ab1）与参考序列（.dna）的基因型判断工具，用于区分 **WT（野生型）**、**HOM（纯合突变）**、**HET（杂合突变）**。

## 功能

- **两步法比对**：先定位参考序列中的最佳匹配片段，再在该片段内做比对并检测不匹配位点。
- **峰数判断**：在差异位点根据色谱峰数（单峰→纯合，双峰→杂合）判定基因型。
- **已知多态性过滤**：支持通过 `known_polymorphisms.txt` 忽略已知多态性位点，减少假阳性。
- **图形界面**：选择参考序列与测序文件（多选或整文件夹），一键分析并导出 CSV。

## 环境要求

- **Python** 3.8+
- **依赖**：Biopython、NumPy；推荐安装 **parasail** 以加速 Smith-Waterman 比对。

## 安装

```bash
# 克隆或下载项目后，进入项目目录
cd GeneType

# 安装依赖（推荐使用 bun 时仍用 pip 安装 Python 包）
pip install -r requirements.txt

# 可选：安装 parasail 以提升比对速度
pip install parasail
```

## 使用

### 图形界面（推荐）

```bash
python genotype_gui.py
```

1. 选择 **参考序列**：`.dna` 文件。
2. 选择 **测序文件**：可多选 `.ab1` 文件，或选择包含 `.ab1` 的文件夹（自动收集并按文件名排序）。
3. 点击 **分析**，结果会显示在列表中，并可导出为 CSV。

### 命令行

```bash
python genotype_analyzer.py <参考序列.dna> <测序文件1.ab1> [测序文件2.ab1 ...]
```

## 数据文件说明

| 文件/格式 | 说明 |
|-----------|------|
| `.dna` | 参考序列（如标准野生型序列）。 |
| `.ab1` | Sanger 测序色谱文件（Applied Biosystems 格式）。 |
| `known_polymorphisms.txt` | 已知多态性位点（参考序列 0-based 位置），每行一个数字；这些位点的差异在分析时会被忽略。可为空或根据 WT 样本结果补充。 |

## 项目结构

| 文件/目录 | 说明 |
|-----------|------|
| `genotype_gui.py` | 图形界面入口。 |
| `genotype_analyzer.py` | 基因型分析核心逻辑（两步法比对、峰数判断、结果输出）。 |
| `alignment_utils.py` | 序列比对工具（Smith-Waterman，支持 parasail / Biopython）。 |
| `pysanger.py` | Sanger 色谱解析（基于 PySanger）。 |
| `requirements.txt` | Python 依赖列表。 |
| `known_polymorphisms.txt` | 已知多态性位点配置。 |
| `GeneType.spec` | PyInstaller 打包配置。 |
| `build_exe.bat` | Windows 下一键打包脚本。 |
| `data/` | 示例数据与验证相关文件（见下方）。 |
| `test_new_example.py` | 对 `data/newExample` 做基因型分析并与 `data/ans2.txt` 比对。 |

### data/ 目录（示例与验证）

| 内容 | 说明 |
|------|------|
| `standard.dna` | 参考序列（示例）。 |
| `40个案例/` | 40 个案例的 .ab1 测序文件。 |
| `newExample/` | 新示例 .ab1 与 .seq 文件。 |
| `79_*.ab1`, `87_*.ab1` | 单文件示例测序。 |
| `results.csv` | 分析结果导出示例。 |
| `ans.txt` | 40个案例的标准答案（与 `genotype_analyzer.py` 默认比对）。 |
| `ans2.txt` | newExample 的标准答案（与 `test_new_example.py` 比对）。 |

## 打包为 Windows exe

在项目目录下执行：

```bat
build_exe.bat
```

完成后可执行文件位于 `dist\GeneType.exe`，可在未安装 Python 的 Windows 电脑上运行。  
详细步骤、环境建议及常见问题见 [打包说明.md](打包说明.md)。

## 输出结果

- **基因型**：WT / HOM / HET / ERR。
- **突变数**：与参考序列的差异位点数（不含已知多态性）。
- **置信度**：基于比对与峰形的置信度。
- 可导出为 CSV，便于后续统计与作图。

## 许可证与引用

- 色谱解析模块 `pysanger.py` 来源于 [PySanger](https://github.com/ponnhide/PySanger)。
- 序列比对推荐使用 [parasail](https://github.com/jeffdaily/parasail)。
