# Endothelial Round2 重新聚类方案

## 1. Round2 目标

Round2 仅针对 **Endothelial lineage** 进行重新聚类。

本轮目的不是提前冻结任何 Endothelial subtype，而是在剔除 Round1 中明确的非 Endothelial 污染 cluster 后，将剩余 Endothelial cells **全部放在一起重新计算低维空间和 clustering**。

重点希望解决：

- Capillary vs Angiogenic/Tip 的边界；
- Round1 中多个 Venous / Inflammatory clusters 是否自然收敛；
- 多个 Arterial clusters 是否自然收敛；
- 低信息量 / stress-like EC 是否会重新归入稳定 vascular subtype；
- Round1 的高分辨率切分是否主要来自 resolution=1.0 的过度细分。

本轮仍然只做到：

- reclustering
- UMAP
- frozen-marker DotPlot
- FindAllMarkers
- specimen / dataset audit
- Round1 → Round2 cluster 来源审计

随后停止，等待人工注释。

---

## 2. Round1 中需要排除的 cluster

从 Endothelial Round2 输入中排除以下 Round1 clusters：

- `7`
- `11`
- `18`
- `20`
- `22`

当前处理理由：

- `7`：T-cell contamination
- `11`：Epithelial contamination
- `18`：Neutrophil / granulocyte contamination
- `20`：mixed / non-Endothelial contamination
- `22`：Macrophage contamination

重要：

**这里的“排除”仅指不进入 Endothelial Round2 reclustering，不代表从完整 Seurat 对象中物理删除这些细胞。**

这些细胞仍保留在完整对象中，后续可由对应 lineage 或人工审阅进一步处理。

---

## 3. Round2 Endothelial 输入范围

除上述 5 个明确非 Endothelial clusters 外，其余 Round1 Endothelial clusters **全部一起进入 Round2**。

即不提前冻结：

- Capillary
- Tip / Angiogenic
- Venous / Inflammatory
- Arterial
- Lymphatic
- Cycling
- Stress-like / unresolved EC

Round2 不预先指定任何 cluster 的最终 subtype。

特别是 Round1 中：

- cluster 2：目前更偏 `Capillary-like / angiogenic capillary`
- cluster 4：目前更偏 `Angiogenic / Tip EC`

但这两个标签在 Round2 **不冻结**。

目标是让去除污染后的 Endothelial cells 重新形成内部结构，再结合 marker 与 cluster mapping 判断二者的关系。

---

## 4. Capillary vs Tip 的判读原则

本轮不依赖单 marker 判断。

### Capillary backbone
更关注稳定 microvascular / capillary program，例如：

- `PLVAP`
- `RGCC`
- `CA4`
- `CD36`
- `EMCN`
- `GPIHBP1`

### Angiogenic / Tip program
更关注显著 angiogenic activation，例如：

- `CXCR4`
- `KDR`
- `ESM1`
- `APLN`
- `DLL4`
- `ANGPT2`

判读逻辑：

> 如果 cluster 以稳定 capillary identity 为主体，同时叠加一定 KDR / ESM1，则优先考虑 `Capillary-like / angiogenic capillary`。

> 如果 angiogenic / tip program 本身就是最突出差异来源，尤其 `CXCR4 / ESM1 / APLN / DLL4` 明显增强，则更偏 `Angiogenic / Tip EC`。

Round2 后再结合新的 cluster structure 与 Round1 → Round2 mapping 判断是否存在：

**Capillary → angiogenic capillary → Tip**

这样的连续谱。

---

## 5. Round2 不新增 marker

Round2 **不新增新的 marker panel**。

继续严格使用已经冻结的 Endothelial framework：

### Pan-Endothelial
- `PECAM1`
- `VWF`
- `CDH5`

### Epithelial contamination / doublet screening
- `EPCAM`
- `KRT8`
- `KRT19`

### Cycling
- `MKI67`
- `TOP2A`
- `UBE2C`

### Angiogenic / Tip EC
- `CXCR4`
- `KDR`
- `ESM1`

### Venous / Inflammatory EC
- `ACKR1`
- `SELE`
- `VCAM1`

### Capillary EC
- `CA4`
- `CD36`
- `RGCC`

### Arterial EC
- `GJA5`
- `EFNB2`
- `HEY1`

### Lymphatic EC
- `PROX1`
- `LYVE1`
- `CCL21`

### MHC-II auxiliary module
- `HLA-DRA`
- `HLA-DPA1`
- `CD74`

### Hypoxia / metallothionein auxiliary module
- `MT1X`
- `MT1E`
- `MT2A`

本轮不新增额外 marker，不修改 frozen marker 顺序。

---

## 6. Round2 重新聚类流程

基于 Round1 已完成的 Endothelial working object。

先根据 `cluster_v2_round1` 去除：

`7, 11, 18, 20, 22`

然后对剩余 Endothelial cells 重新执行：

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
   - **resolution = 0.6**

7. `RunUMAP`
   - reduction = Harmony
   - dims = 1:30

8. Frozen-marker DotPlot

9. `FindAllMarkers`

10. specimen / dataset composition audit

11. Round1 → Round2 cluster mapping audit

---

## 7. 为什么 Round2 使用 resolution = 0.6

Round1 使用 resolution=1.0 的主要价值已经实现：

- 暴露 Endothelial 内部异质性；
- 暴露非 Endothelial contamination；
- 显示多个 Venous / Arterial / Tip / Capillary states。

在明确 contamination 被排除后，Round2 更强调：

- 减少同一 vascular subtype 的过度切分；
- 观察多个 Round1 cluster 是否自然合并；
- 获得更适合后续 cell2location 的中等分辨率 taxonomy。

因此 Round2 clustering resolution 固定为：

**0.6**

---

## 8. Round2 cluster 字段

新增 Round2 cluster metadata，例如：

`cluster_v2_round2`

必须同时保留：

`cluster_v2_round1`

不得覆盖 Round1 cluster ID。

---

## 9. Round1 → Round2 cluster mapping audit

这是 Round2 的关键审计。

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

特别关注：

- Round1 `4 / 8 / 21` 是否重新聚成统一 Tip / Angiogenic cluster；
- Round1 `2 / 5 / 6` 是否重新聚成 Capillary cluster；
- Round1 `0 / 1 / 3 / 10 / 13 / 15` 是否收敛为 Venous / Inflammatory EC；
- Round1 `9 / 12 / 19` 是否收敛为 Arterial EC；
- Round1 `17` 是否自然并入稳定 vascular subtype；
- Round1 `14` 是否仍维持独立 Lymphatic identity。

---

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

第 4 张图需要：

- 使用 Round2 UMAP 坐标；
- 按 Round1 cluster 上色。

用于观察：

- Round1 cluster 是否被重新打散；
- 多个相似 vascular states 是否重新合并；
- Capillary 与 Tip 是否仍然保持稳定边界；
- stress / low-information clusters 是否重新归入主 lineage。

UMAP 沿用已经修正后的可视化风格：

- 点大小足够
- alpha 足够高
- cluster 云团清楚可见
- 白底
- cluster label 清晰

---

## 11. DotPlot 输出

继续使用 frozen Endothelial marker panel。

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
- stress / hypoxia-like cluster 是否有明显 dataset bias

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

1. 提前冻结任何 Endothelial subtype
2. 自动把 cluster 2 判为 Capillary
3. 自动把 cluster 4 判为 Tip
4. 新增 marker panel
5. 自动 merge clusters
6. 自动删除新的疑难 EC cluster
7. 修改完整对象中的正式 `subtype_v2`
8. 启动 cell2location
9. 重跑其他 lineage
10. 覆盖 Round1 结果

---

## 16. Round2 输出目录

在 v2 下创建：

`round2/Endothelial/`

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

Endothelial contaminant exclusion  
→ Endothelial reclustering  
→ Round2 UMAP  
→ frozen-marker DotPlot  
→ FindAllMarkers  
→ specimen/dataset audit  
→ Round1→Round2 mapping audit  
→ annotation decision template

然后等待人工注释。

最终状态写：

`ENDO_ROUND2_STATUS = READY_FOR_MANUAL_ANNOTATION`

如果失败：

`ENDO_ROUND2_STATUS = PARTIAL_REVIEW_REQUIRED`
