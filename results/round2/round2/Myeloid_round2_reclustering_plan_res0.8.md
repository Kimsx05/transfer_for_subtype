# Myeloid Round2 重新聚类方案

## 1. Round2 目标

Round2 仅针对 **Myeloid lineage** 进行重新聚类。

本轮目的不是提前冻结全部 Myeloid subtype，而是在剔除 Round1 中明确的非 Myeloid contamination 后，将剩余 Myeloid cells **全部放在一起重新计算低维空间和 clustering**。

本轮重点：

- 去除明确 contamination；
- 将目前偏碎的 Myeloid clustering 适度收敛；
- 保留 Mono / TAM / DC / Neutrophil 等主要结构；
- **clustering resolution 从 1.0 调整为 0.8**；
- 不新增 marker，不修改 frozen Myeloid marker framework；
- 输出 Round1 → Round2 cluster mapping，判断不同 Round1 states 是否稳定、合并或重组。

完成后停止，等待人工注释。

---

## 2. Round1 中需要排除的 cluster

从 Myeloid Round2 输入中排除以下 Round1 clusters：

- `5`
- `11`
- `19`
- `23`

当前处理理由：

- `5`：T-cell contamination
- `11`：Epithelial contamination
- `19`：Epithelial contamination
- `23`：ILC / γδ-T-like non-Myeloid contamination

重要：

**这里的“排除”仅指不进入 Myeloid Round2 reclustering，不代表从完整 Seurat 对象中物理删除这些细胞。**

这些细胞仍保留在完整对象中，后续可由对应 lineage 或人工审阅进一步处理。

---

## 3. Round2 Myeloid 输入范围

除上述 4 个明确 contamination clusters 外，其余 Round1 Myeloid clusters **全部一起进入 Round2**。

不提前冻结最终 subtype。

当前重点关注的 Round1 groups：

### Classical Monocyte family
- `1`
- `2`
- `4`
- `6`
- `17`
- `21`

其中：

- `2`
- `17`
- `21`

目前统一视为 **Classical Monocyte family**，暂不进一步细分。

---

### C1QC-family TAM
- `7`
- `15`

当前不强行合并。

Round1 中：

- `7`：偏 activated / APC-like C1QC TAM
- `15`：偏 C1QC / TREM2 TAM

由于两者在 UMAP 上并不直接相连，且中间存在 SPP1 TAM cluster，因此 Round2 需要观察：

> 两者在去除 contamination、降低 resolution 后是否仍保持独立。

---

### SPP1 TAM
- `8`

当前为较稳定 SPP1 TAM candidate。

---

### FOLR2 / LYVE1-like TAM
- `0`

当前属于 conditional subtype，不属于必须保留的 core taxonomy。

Round2 重点判断：

- 是否仍形成稳定独立 cluster；
- 是否跨 specimen / dataset 稳定；
- 是否最终需要保留为独立 `FOLR2_TAM`；
- 或是否自然并入 C1QC TAM family。

---

### Transitional TAM states
- `13`
- `22`

当前理解：

- `13`：ISG-high TAM transitional state
- `22`：MERTK / SPP1 / macrophage transitional state

本轮不要求对这两个 cluster 做精细最终命名。

Round2 重点观察：

- 是否继续形成稳定独立群；
- 是否自然并回 C1QC / SPP1 / other TAM states。

---

### Cluster 10
- `10`

当前为 IL3RA-high / unusual monocyte-like population。

主要问题：

- marker phenotype 特殊；
- Round1 显示明显 dataset bias；
- 当前不做最终 subtype 定义。

Round2 重点判断：

- 是否仍然形成独立 cluster；
- 是否仍主要由同一 dataset 驱动；
- 是否在重新聚类后自然并入 Mono / DC / other Myeloid group。

---

### DC / Neutrophil / pDC / Cycling groups

其余 Round1 中较清晰的结构继续全部保留进入 Round2，包括：

- cDC2
- cDC1
- LAMP3 DC
- pDC
- Neutrophil
- IFN-Neutrophil
- Cycling TAM / macrophage

本轮不提前冻结。

---

## 4. Round2 clustering 参数

Round2 使用：

- RNA assay
- HVG = 3000
- PCA = 50 PCs
- Harmony batch variable = `orig.ident`
- downstream dims = 1:30
- UMAP based on Harmony
- **FindClusters resolution = 0.8**

---

## 5. 为什么使用 resolution = 0.8

Round1 resolution=1.0 已成功暴露：

- Mono internal states
- TAM heterogeneity
- DC populations
- Neutrophil populations
- contamination

但当前 clustering 有一定过度碎裂。

Myeloid 本身内部异质性较高，因此本轮不降到 0.6，而选择：

**resolution = 0.8**

目的：

- 适度减少过度切分；
- 保留 Mono / TAM / DC / Neutrophil 的主要结构；
- 观察若干 transitional states 是否自然并回主群；
- 避免把真实 Myeloid heterogeneity 过度压缩。

---

## 6. DotPlot marker framework

Round2 **不新增 marker**。

继续严格使用 frozen Myeloid marker framework。

### Pan-myeloid
- `LYZ`
- `LST1`
- `TYROBP`

### Epithelial contamination / doublet screen
- `EPCAM`
- `KRT8`
- `KRT19`

### Cycling
- `MKI67`
- `TOP2A`
- `UBE2C`

### Neutrophil
- `FCGR3B`
- `CSF3R`
- `CXCR2`

### Eosinophil
- `CLC`
- `CCR3`
- `ALOX15`

### Basophil
- `CD200R1`
- `CLEC12A`
- `KLK10`

### Classical Monocyte
- `FCN1`
- `S100A8`
- `VCAN`

### Non-classical Monocyte
- `FCGR3A`
- `CDKN1C`
- `LST1`

### C1QC TAM
- `C1QC`
- `C1QA`
- `APOE`

### SPP1 TAM
- `SPP1`
- `GPNMB`
- `FN1`

### ISG TAM
- `ISG15`
- `CXCL10`
- `GBP1`

### FOLR2 / LYVE1-like TAM
- `FOLR2`
- `LYVE1`
- `SELENOP`

### TREM2 auxiliary
- `TREM2`
- `APOC1`
- `GPNMB`

### cDC1
- `CLEC9A`
- `XCR1`
- `CADM1`

### cDC2
- `CD1C`
- `FCER1A`
- `CLEC10A`

### LAMP3 DC
- `LAMP3`
- `CCR7`
- `FSCN1`

### pDC
- `LILRA4`
- `IL3RA`
- `GZMB`

本轮不修改 marker 顺序，不增加额外 marker panel。

---

## 7. Round2 重新聚类流程

基于 Round1 已完成的 Myeloid working object。

先根据 `cluster_v2_round1` 去除：

`5, 11, 19, 23`

然后对剩余 Myeloid cells 重新执行：

1. `FindVariableFeatures`
   - nfeatures = 3000

2. `ScaleData`

3. `RunPCA`
   - npcs = 50

4. `RunHarmony`
   - batch variable = `orig.ident`

5. `FindNeighbors`
   - reduction = Harmony
   - dims = 1:30

6. `FindClusters`
   - **resolution = 0.8**

7. `RunUMAP`
   - reduction = Harmony
   - dims = 1:30

8. Frozen-marker DotPlot

9. `FindAllMarkers`

10. specimen / dataset composition audit

11. Round1 → Round2 cluster mapping audit

---

## 8. Round2 cluster 字段

新增 Round2 cluster metadata，例如：

`cluster_v2_round2`

必须同时保留：

`cluster_v2_round1`

不得覆盖 Round1 cluster ID。

---

## 9. Round1 → Round2 cluster mapping audit

需要明确回答：

> 每一个 Round2 cluster 是由哪些 Round1 clusters 重新组成的？

至少输出以下文件。

### 9.1 Counts matrix

`round1_to_round2_cluster_counts.csv`

行：

`cluster_v2_round1`

列：

`cluster_v2_round2`

值：

对应 cell number。

---

### 9.2 Row-normalized proportions

`round1_to_round2_cluster_row_proportions.csv`

用于回答：

> 一个 Round1 cluster 的细胞，在 Round2 被重新分配到了哪些 clusters？

---

### 9.3 Column-normalized proportions

`round1_to_round2_cluster_column_proportions.csv`

用于回答：

> 一个新的 Round2 cluster，主要由哪些 Round1 clusters 组成？

---

### 9.4 本轮重点追踪

特别关注以下关系：

- `2 / 17 / 21` 是否与其他 Classical Mono 自然收敛；
- `7 / 15` 是否仍维持两个独立 C1QC-family TAM clusters；
- `0` 是否保持独立 FOLR2 / LYVE1-like TAM；
- `13 / 22` 是否自然并入主 TAM populations；
- `10` 是否仍保持独立，且是否仍表现明显 dataset bias；
- `8` 是否保持稳定 SPP1 TAM；
- DC / Neutrophil populations 是否保持独立稳定结构。

---

### 9.5 可视化

输出：

- `round1_to_round2_cluster_heatmap.pdf/png`

如果方便，再输出：

- `round1_to_round2_alluvial.pdf/png`

---

## 10. UMAP 输出

至少输出：

- `01_umap_cluster_round2.pdf/png`
- `02_umap_dataset_round2.pdf/png`
- `03_umap_specimen_round2.pdf/png`
- `04_umap_round1_cluster_on_round2_embedding.pdf/png`

第 4 张图使用：

- Round2 UMAP 坐标；
- Round1 cluster 上色。

用于观察：

- Round1 clusters 是否被重新打散；
- Mono clusters 是否收敛；
- C1QC TAM 7 / 15 是否仍分离；
- transitional TAM 是否自然并入主群；
- cluster 10 是否仍形成 dataset-specific island。

UMAP 沿用既定可视化风格：

- 点大小足够
- alpha 足够高
- cluster 云团清楚可见
- 白底
- cluster label 清晰

---

## 11. DotPlot 输出

继续使用 frozen Myeloid marker panel。

要求：

- 纵轴 = Round2 cluster
- 横轴 = frozen marker
- 顶部保留 marker group annotation
- dot size = Percent Expressed
- color = Average Expression
- 白底
- marker group 清楚分隔

输出：

`05_dotplot_frozen_markers_round2.pdf/png`

不得新增 marker，不得修改 frozen marker 顺序。

---

## 12. FindAllMarkers

基于 Round2 clusters 输出：

- `markers_all_round2.csv`
- `markers_top10_per_cluster_round2.csv`
- `markers_top20_per_cluster_round2.csv`

保留完整统计字段。

不根据 marker 自动完成 subtype annotation。

---

## 13. specimen / dataset audit

继续输出：

- `cluster_by_specimen_counts_round2.csv`
- `cluster_by_specimen_proportions_round2.csv`
- `cluster_by_dataset_counts_round2.csv`
- `cluster_by_dataset_proportions_round2.csv`

重点关注：

- 单一 specimen 驱动
- 单一 dataset 驱动
- 极小 cluster
- 仅少数 specimen 出现的 cluster
- cluster 10 是否仍明显 dataset-specific
- conditional TAM states 是否跨 dataset 稳定

这些均作为 warning，不自动删除。

---

## 14. annotation decision template

生成：

`annotation_decision_template_round2.csv`

字段至少：

- `cluster_round2`
- `n_cells`
- `major_round1_sources`
- `proposed_subtype`
- `confidence`
- `merge_to`
- `notes`

其中：

- cluster / n_cells 自动填写
- `major_round1_sources` 自动填写主要 Round1 来源
- 其他人工注释字段保持空白

---

## 15. 本轮禁止事项

Round2 禁止：

1. 提前冻结所有 Myeloid subtype
2. 自动把 cluster 7 和 15 合并
3. 自动把 cluster 0 定为最终 FOLR2 TAM
4. 自动给 13 / 22 过度细化命名
5. 自动删除 cluster 10
6. 新增 marker panel
7. 自动 merge clusters
8. 修改完整对象中的正式 `subtype_v2`
9. 启动 cell2location
10. 重跑其他 lineage
11. 覆盖 Round1 结果
12. 将 resolution 改为其他值

---

## 16. Round2 输出目录

在 v2 下创建：

`round2/Myeloid/`

并保留：

- code
- logs
- figures
- markers
- audit
- clustered qs

不要预建后续 round。

---

## 17. 强制停止点

完成以下内容后立即停止：

Myeloid contaminant exclusion  
→ Myeloid reclustering  
→ Round2 UMAP  
→ frozen-marker DotPlot  
→ FindAllMarkers  
→ specimen/dataset audit  
→ Round1→Round2 mapping audit  
→ annotation decision template

然后等待人工注释。

最终状态写：

`MYELOID_ROUND2_STATUS = READY_FOR_MANUAL_ANNOTATION`

如果失败：

`MYELOID_ROUND2_STATUS = PARTIAL_REVIEW_REQUIRED`
