# Co-selection

本仓库是论文 **《Understanding pollutant-driven shifts of antibiotic resistome in activated
sludge: A lab-scale study》**（Li P., Yin Z., Ye L.，*Journal of Hazardous Materials*
**487** (2025) 137199，[doi:10.1016/j.jhazmat.2025.137199](https://doi.org/10.1016/j.jhazmat.2025.137199)）
的配套代码。

仓库覆盖了从原始测序数据到论文图件的全部脚本：16S rRNA 基因扩增子数据的 usearch 流程、
宏基因组的 GTDB-Tk / diamond / ARGs-OAP 流程，以及产出统计结果和图件的三个 R 分析脚本。

> English version: **[README.md](README.md)**

---

## 研究概要

实验使用来自工业园区污水处理厂的活性污泥，在 21 个实验室规模反应器（7 个处理 × 3 个平行）
中运行：3 种重金属（Cr(VI)、Cd(II)、Cu(II)）和 3 种芳香族有机物（苯酚、4-氨基苯酚、
4-硝基苯酚），外加空白对照，第 1 阶段（0–30 d）投加污染物，第 2 阶段（30–60 d）停止投加。
共采集 15 个样品（14 个污泥样品 + 种子污泥），分别进行 16S rRNA 基因扩增子测序和
鸟枪法宏基因组测序。

* 重金属和有机物处理都使 ARG 总丰度升高（相对种子污泥为 1.7–2.9 倍），其中 `bacA` 和
  `sul1` 始终是丰度最高的 ARG。
* ARG 组成在重金属组和有机物组之间明显分离（NMDS，stress = 0.11）；ARG 丰度的升高与
  污染物引起的细菌群落变化同步，而不是污染物直接作用的结果。
* 共现网络与 Mantel 检验表明：ARG 通过金属抗性基因（MRG）和芳香族降解基因（ADG）被
  **间接共选择**，且主要集中在 γ-变形菌纲。
* ARG 总丰度升高，但高风险 ARG 的比例并未上升——丰度高并不等于健康风险高。

## 仓库结构

```
Co-selection/
├── pipelines/          # 上游生物信息流程（shell，含一个 python 辅助脚本）
│   ├── 01_16S_amplicon/                usearch：拼接、去引物、UNOISE3、ASV 表、SINTAX
│   ├── 02_metagenome_binning/          GTDB-Tk 对 257 个 MAG 的物种注释与建树（Fig. 3d）
│   ├── 03_resistance_gene_annotation/  diamond 比对 SARG/MRG/ADG/MGE，生成 MAG 基因表（Fig. 3d, 3e）
│   └── 04_gene_abundance/              ARGs-OAP v3.2.2 定量（Fig. 1, 2, 5）
├── analyses/           # 下游统计与绘图（R）
│   ├── 01_16S_community/               Fig. 6、Fig. S7、Fig. S8
│   ├── 02_metagenome_resistome/        Fig. 1、2、3c、5，Fig. S2–S6
│   └── 03_MAG_gene_annotation/         Fig. 3d、Fig. 4、Fig. S4
├── data/               # 小型输入表格（见 data/README.md）
└── docs/               # 流程、目录、文件对应、输入输出、软件环境
```

两层脚本都按执行顺序编号；每个文件头部都写明了用途、对应的论文图件、主要输出以及它在
原始工作目录中的路径。两个较长的 R 脚本原本是单一的连续 notebook，因此额外提供了
"Section index"（章节目录及行号）；[`docs/04`](docs/04_script_inputs_outputs.md) 把它整理成了表格。

## 文档

| 文档 | 内容 |
| --- | --- |
| [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md) | 端到端流程，逐阶段列出输入与输出 |
| [docs/02_directory_layout.md](docs/02_directory_layout.md) | 目录说明、命名约定、环境变量 |
| [docs/03_file_mapping.md](docs/03_file_mapping.md) | 每个脚本：原路径 → 新路径，以及做了哪些改动 |
| [docs/04_script_inputs_outputs.md](docs/04_script_inputs_outputs.md) | 各脚本的输入、输出，以及"章节 → 图件"对照 |
| [docs/05_software_and_environment.md](docs/05_software_and_environment.md) | 软件版本、conda 环境、R 包、数据库 |

## 复现

**1. 上游流程。** `pipelines/` 下的脚本就是实际使用的作业脚本（保留了 `#SBATCH` 头）。
除路径外与生产运行完全一致：所有站点相关的绝对路径都改成了变量，每个脚本开头都有
`# Site-specific settings` 段（如 `ORF_DIR`、`SARG_DB`、`DIAMOND_DIR`、`IN_DIR`）。
设置这些变量后，按 [docs/01](docs/01_pipeline_overview.md) 给出的顺序提交即可。

**2. R 分析。** `analyses/` 下的脚本需要 R >= 4.2 以及
[docs/05](docs/05_software_and_environment.md) 中列出的 R 包。输入表通过 `CS_DATA` 定位，
图件写入 `CS_WORK`；两者默认指向仓库内的目录，也可以在不修改脚本的情况下覆盖：

```r
Sys.setenv(CS_DATA = "/path/to/data", CS_WORK = "/path/to/results")
source("analyses/02_metagenome_resistome/01_resistome_analysis.R")
```

请从脚本开头顺序运行：同一个 notebook 内的后续章节会复用前面创建的对象（例如第 3 节读入的
ARG type 表，在 NMDS、热图和风险等级分析中还会继续使用）。

## 数据可用性

复现分析所需的小型表格（ASV 表、SINTAX 注释、MAG 丰度、GTDB-Tk 注释、质粒统计、样品信息、
LEfSe 输入表）已随仓库提供，见 [`data/`](data/README.md)。较大的中间表（ARGs-OAP 输出、
diamond 比对结果、MAG × 基因矩阵）未随仓库分发，可用 `pipelines/` 中的脚本从原始 reads
重新生成。原始测序数据按论文中的说明可向通讯作者索取。所用数据库（SARG、BacMet、
AromaDeg、MGE、MiDAS 4.8.1、GTDB R214）均为公开数据库。

## 引用

```bibtex
@article{Li2025coselection,
  title   = {Understanding pollutant-driven shifts of antibiotic resistome in activated
             sludge: A lab-scale study},
  author  = {Li, Pengwei and Yin, Zirui and Ye, Lin},
  journal = {Journal of Hazardous Materials},
  volume  = {487},
  pages   = {137199},
  year    = {2025},
  doi     = {10.1016/j.jhazmat.2025.137199}
}
```

使用本仓库脚本时请引用上述论文。论文方法部分（Methods 2.3）也给出了本仓库地址。
