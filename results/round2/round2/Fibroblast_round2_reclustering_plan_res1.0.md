# Fibroblast Round2 重新聚类方案

## 1. Round2 目标

Round2 仅针对 **Fibroblast / stromal lineage** 进行重新聚类。

本轮目的不是提前冻结任何最终 CAF subtype，而是在剔除 Round1 中明确的非 Fibroblast contamination 后，将剩余 stromal cells **全部放在一起重新计算低维空间和 clustering**。

本轮重点：

- 去除明确的 T / epithelial / myeloid / Schwann / mixed contamination；
- 保留剩余 fibroblast / CAF / mural-like populations 全部一起重跑；
- **clustering resolution 继续保持 1.0，不降低**；
- DotPlot 保留 frozen fibroblast marker framework；
- 仅新增一个 `MyoFibro / SMC` marker block；
- 继续输出 Round1 → Round2 cluster mapping，辅助判断不同 stromal states 是否稳定。

完成后停止，等待人工注释。

---

## 2. Round1 中需要排除的 cluster

从 Fibroblast Round2 输入中排除以下 Round1 clusters：

- `10`
- `11`
- `12`
- `19`
- `21`
- `23`

当前处理理由：

- `10`：T-cell contamination
- `11`：Epithelial contamination
- `12`：Myeloid contamination
- `19`：Epithelial / mixed contamination
- `21`：Schwann / peripheral glial contamination
- `23`：T-cell + mural mixed / doublet-like contamination

重要：

**这里的“排除”仅指不进入 Fibroblast Round2 reclustering，不代表从完整 Seurat 对象中物理删除这些细胞。**

这些细胞仍保留在完整对象中，后续可由对应 lineage 或人工审阅进一步处理。

---

## 3. Round2 Fibroblast 输入范围

除上述 6 个明确 contamination clusters 外，其余 Round1 Fibroblast clusters **全部一起进入 Round2**。

不提前冻结：

- iCAF
- Matrix_CAF
- rCAF
- vCAF
- Pericyte
- MyoFibro / SMC
- Cycling stromal
- stress-like / unresolved fibroblast
- transitional fibroblast states

Round2 不预先指定任何 cluster 的最终 subtype。

特别是：

- Round1 cluster 3 / 18 目前均表现出明显 contractile / smooth-muscle-like program；
- Round1 cluster 5 偏 Pericyte；
- Round1 cluster 22 仍保留进入 Round2，不作为 contamination 删除。

这些均留到 Round2 后结合 all markers 与 cluster mapping 再判断。

---

## 4. Round2 clustering 参数

Round2 继续使用：

- RNA assay
- HVG = 3000
- PCA = 50 PCs
- Harmony batch variable = `orig.ident`
- downstream dims = 1:30
- UMAP based on Harmony
- **FindClusters resolution = 1.0**

本轮不降低 resolution。

原因：

- Fibroblast / CAF 内部异质性本身较高；
- 当前希望在去除 contamination 后继续保留较细 stromal structure；
- MyoFibro / SMC / Pericyte / Matrix / iCAF 等可能存在连续但可区分的状态；
- 因此 Round2 继续用 resolution=1.0，避免过早将真实 stromal heterogeneity 合并。

---

## 5. DotPlot marker framework

Round2 继续使用 frozen Fibroblast marker framework。

### Pan-fibroblast
- `COL1A1`
- `DCN`
- `LUM`

### Epithelial contamination / doublet screen
- `EPCAM`
- `KRT8`
- `KRT19`

### Cycling
- `MKI67`
- `TOP2A`
- `UBE2C`

### Matrix_CAF
- `COL1A1`
- `MMP11`
- `POSTN`

### iCAF
- `PLA2G2A`
- `CFD`
- `CXCL12`

### vCAF
- `MCAM`
- `NOTCH3`
- `COL18A1`

### Pericyte
- `RGS5`
- `CSPG4`
- `PDGFRB`

### apCAF
- `CD74`
- `HLA-DRA`
- `CIITA`

### tCAF
- `MME`
- `NDRG1`
- `ENO1`

### ifnCAF
- `IDO1`
- `CXCL10`
- `IFIT1`

### rCAF
- `CCL19`
- `CCL21`
- `COL14A1`

### dCAF / proliferating fibroblast
- `MKI67`
- `TOP2A`
- `UBE2C`

---

## 6. 新增 DotPlot marker block

在原 frozen fibroblast DotPlot 基础上，**仅新增以下一个 marker group**：

### MyoFibro / SMC
- `ACTA2`
- `TAGLN`
- `MYH11`

要求：

- 这三个 marker 必须作为一个独立顶部注释 block；
- block 名称固定为：

`MyoFibro / SMC`

- 不再拆成两个 marker groups；
- 本轮不新增其他 MyoFibro 或 SMC marker；
- 不修改其他 frozen marker 的顺序与内容。

说明：

这个 block 仅用于显示 contractile stromal axis，**不代表 Round2 最终一定将 MyoFibro 与 SMC 合并为同一个 biological subtype**。

最终是否需要区分：

- MyoFibro
- SMC

留到 Round2 all markers + cluster mapping 后再人工决定。

---

## 7. Round2 重新聚类流程

基于 Round1 已完成的 Fibroblast working object。

先根据 `cluster_v2_round1` 去除：

`10, 11, 12, 19, 21, 23`

然后对剩余 cells 重新执行：

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
   - **resolution = 1.0**

7. `RunUMAP`
   - reduction = Harmony
   - dims = 1:30

8. Frozen-marker DotPlot
   - 原 marker framework
   - 加入 `MyoFibro / SMC = ACTA2 / TAGLN / MYH11`

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

### 9.2 Row-normalized proportions

`round1_to_round2_cluster_row_proportions.csv`

用于回答：

> 一个 Round1 cluster 的细胞，在 Round2 被重新分配到了哪些 clusters？

### 9.3 Column-normalized proportions

`round1_to_round2_cluster_column_proportions.csv`

用于回答：

> 一个新的 Round2 cluster，主要由哪些 Round1 clusters 组成？

特别关注：

- Round1 `3 / 18` 是否继续维持 contractile / MyoFibro-SMC-like identity；
- Round1 `5` 是否仍保持独立 Pericyte identity；
- Round1 `9 / 15 / 22` 是否重新组织为 Matrix_CAF；
- Round1 `0 / 1 / 4 / 13` 是否重新组织为 iCAF；
- Round1 `2 / 7 / 8 / 14 / 16 / 17` 是否自然并入稳定 CAF states，或继续形成独立群；
- Round1 `20` 是否保持独立 Cycling stromal identity。

### 9.4 可视化

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

第 4 张图使用 Round2 UMAP 坐标、按 Round1 cluster 上色。

用于观察：

- Round1 cluster 是否被重新打散；
- contractile stromal populations 是否仍然独立；
- iCAF / Matrix_CAF 是否保持稳定；
- unresolved clusters 是否重新并入主 stromal populations。

UMAP 沿用已经修正后的可视化风格：

- 点大小足够
- alpha 足够高
- cluster 云团清楚可见
- 白底
- cluster label 清晰

---

## 11. DotPlot 输出

继续使用 frozen Fibroblast marker panel，并新增：

`MyoFibro / SMC = ACTA2 / TAGLN / MYH11`

要求：

- 纵轴 = Round2 cluster
- 横轴 = marker
- 顶部保留 marker group annotation
- dot size = Percent Expressed
- color = Average Expression
- 白底
- marker group 清楚分隔
- `MyoFibro / SMC` 作为独立 marker block 展示

输出：

`05_dotplot_frozen_markers_plus_myo_smc_round2.pdf/png`

除上述 3 个 marker 外，不新增其他 marker。

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
- stress / unresolved fibroblast 是否有明显 dataset bias

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

1. 提前冻结任何 Fibroblast / CAF subtype
2. 自动将 MyoFibro 与 SMC 合并为最终 subtype
3. 自动删除新的疑难 stromal cluster
4. 除 `ACTA2 / TAGLN / MYH11` 外新增 marker
5. 修改其他 frozen marker
6. 自动 merge clusters
7. 修改完整对象中的正式 `subtype_v2`
8. 启动 cell2location
9. 重跑其他 lineage
10. 覆盖 Round1 结果
11. 将 resolution 改为 0.6

---

## 16. Round2 输出目录

在 v2 下创建：

`round2/Fibroblast/`

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

Fibroblast contaminant exclusion  
→ Fibroblast reclustering  
→ Round2 UMAP  
→ frozen-marker DotPlot + MyoFibro/SMC block  
→ FindAllMarkers  
→ specimen/dataset audit  
→ Round1→Round2 mapping audit  
→ annotation decision template

然后等待人工注释。

最终状态写：

`FIBRO_ROUND2_STATUS = READY_FOR_MANUAL_ANNOTATION`

如果失败：

`FIBRO_ROUND2_STATUS = PARTIAL_REVIEW_REQUIRED`
