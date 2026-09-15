# B cell Round2 重新聚类方案

## 1. Round2 目标

Round2 仅针对 **B cell lineage** 进行重新聚类。

本轮目的不是提前冻结某些 B subtype，而是在剔除 Round1 中明确的非 B 污染 cluster 后，将其余 B cells **全部放在一起重新计算低维空间和 clustering**，观察 Naïve / Memory / IFN / GC / Activated 等状态能否得到更清晰、稳定的分离。

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

从 B Round2 输入中排除以下 Round1 clusters：

- `12`
- `13`
- `14`
- `16`
- `17`

当前处理理由：

- `12`：明显非 B，T-like
- `13`：明确不是 B，也不是 epithelial；具体身份本轮不继续追究，直接从 B Round2 排除
- `14`：明显非 B，T/NK / γδ-like
- `16`：Epithelial contamination
- `17`：明显非 B，Neutrophil / granulocyte-like

重要：

**这里的“排除”仅指不进入 B Round2 reclustering，不代表从完整 Seurat 对象中物理删除这些细胞。**

这些细胞后续仍保留在完整对象中，可由对应 major lineage 或人工审阅进一步处理。

---

## 3. Round2 B 输入范围

除上述 5 个明确非 B clusters 外，Round1 中其余 B clusters **全部一起进入 Round2**。

包括：

`0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15`

其中：

- `10`
- `11`

虽然 Round1 均表现出明显 GC B program，但 **本轮不冻结、不提前写死为 GC B，也不从 reclustering 中移除**。

Round2 仍将它们与其他 B cells 一起重新聚类。

原则：

> 不提前锁定任何 B subtype，让剩余 B cells 在去除明显污染后重新形成内部结构。

---

## 4. Round2 不新增 marker

Round2 **不新增新的 marker panel**。

继续严格使用已经冻结的 B-cell framework：

### Pan-B
- `CD79A`
- `MS4A1`
- `CD19`

### Epithelial contamination / doublet screening
- `EPCAM`
- `KRT8`
- `KRT19`

### Naïve B
- `IGHD`
- `TCL1A`
- `FCER2`

### Memory B
- `CD27`
- `TNFRSF13B`
- `AIM2`

### GC B
- `BCL6`
- `AICDA`
- `RGS13`

### IFN-stimulated B
- `ISG15`
- `STAT1`
- `MX1`

### Cycling B
- `MKI67`
- `TOP2A`
- `UBE2C`

本轮不额外加入：

- Activated B marker panel
- Atypical memory B marker panel
- Stress B marker panel

如这些状态真实存在，应优先由：

- Round2 cluster structure
- `FindAllMarkers`
- frozen DotPlot
- specimen / dataset stability

自然体现出来，再由人工决定是否建立新的 subtype 名称。

---

## 5. Round2 重新聚类流程

基于 Round1 已完成对象，优先复用已有 B lineage working object。

### 输入
Round1 B clustered object。

先根据 `cluster_v2_round1` 去除：

`12, 13, 14, 16, 17`

然后对剩余 B cells 重新执行：

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
   - resolution = 0.6

Round2 将 clustering granularity 调低至 0.6，原因是 Round1 已完成污染 cluster 暴露与排除；本轮更强调在纯化后的 B lineage 内减少过度切分。

7. `RunUMAP`
   - reduction = Harmony
   - dims = 1:30

8. Frozen-marker DotPlot

9. `FindAllMarkers`

10. specimen / dataset composition audit

11. Round1 → Round2 cluster mapping audit

---

## 6. Round2 cluster 字段

新增清晰的 Round2 cluster metadata，例如：

`cluster_v2_round2`

必须同时保留：

`cluster_v2_round1`

用于追踪 cluster 重组过程。

不得覆盖 Round1 cluster ID。

---

## 7. Round1 → Round2 cluster mapping audit

这是 Round2 的新增重点审计。

需要明确回答：

> 每一个 Round2 cluster 是由哪些 Round1 clusters 重新组成的？

至少输出以下文件。

### 7.1 Counts matrix

`round1_to_round2_cluster_counts.csv`

行：

`cluster_v2_round1`

列：

`cluster_v2_round2`

值：

对应 cell number。

---

### 7.2 Row-normalized proportions

`round1_to_round2_cluster_row_proportions.csv`

用于回答：

> 一个 Round1 cluster 的细胞，在 Round2 被重新分配到了哪些 clusters？

例如：

Round1 cluster 2：
- 70% → Round2 cluster A
- 25% → Round2 cluster B
- 5% → Round2 cluster C

---

### 7.3 Column-normalized proportions

`round1_to_round2_cluster_column_proportions.csv`

用于回答：

> 一个新的 Round2 cluster，主要来源于哪些 Round1 clusters？

例如：

Round2 cluster 5：
- 45% 来自 Round1 cluster 0
- 30% 来自 Round1 cluster 6
- 20% 来自 Round1 cluster 15
- 5% 其他

这可以帮助判断：

- 新 cluster 是否是多个旧 cluster 的合理重组；
- 是否仍然主要复制 Round1 的原始边界；
- 某个新 subtype 是否有多个 Round1 cluster 共同支持；
- 是否出现一个旧 cluster 被重新拆成多个生物学上更清晰的新群。

---

### 7.4 可视化

建议输出：

`round1_to_round2_cluster_heatmap.pdf/png`

如果方便，再输出：

`round1_to_round2_alluvial.pdf/png`

用于直观看 Round1 → Round2 的 cell reassignment。

---

## 8. UMAP 输出

至少输出：

- `01_umap_cluster_round2.pdf/png`
- `02_umap_dataset_round2.pdf/png`
- `03_umap_specimen_round2.pdf/png`
- `04_umap_round1_cluster_on_round2_embedding.pdf/png`

第 4 张图非常重要：

在 **Round2 UMAP 坐标** 上用 Round1 cluster 上色。

这样可以直观看：

- Round1 cluster 是否被重新打散；
- 原来的 GC / IFN / Naïve / Memory 等区域是否在新的 embedding 中重新组织；
- 某些 Round1 cluster 是否只是技术性边界。

UMAP 继续使用已经修正后的可视化风格：

- 点大小足够
- alpha 足够高
- cluster 云团清楚可见
- 白底
- cluster label 清晰

---

## 9. DotPlot 输出

继续使用 frozen B marker panel。

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

## 10. FindAllMarkers

基于 Round2 clusters 输出：

- `markers_all_round2.csv`
- `markers_top10_per_cluster_round2.csv`
- `markers_top20_per_cluster_round2.csv`

保留完整统计字段。

不根据 marker 自动完成 subtype annotation。

---

## 11. specimen / dataset audit

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

这些均作为 warning，不自动删除。

---

## 12. annotation decision template

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

## 13. 本轮禁止事项

Round2 禁止：

1. 冻结 10/11 为 GC B 后跳过 reclustering
2. 提前给任何 Round2 cluster 写最终 subtype
3. 新增 marker panel
4. 自动创建 Activated B / Atypical B / Stress B 标签
5. 自动 merge clusters
6. 自动删除新的疑难 cluster
7. 修改完整对象的正式 `subtype_v2`
8. 启动 cell2location
9. 重跑其他 lineage
10. 覆盖 Round1 结果

---

## 14. Round2 输出目录

在 v2 下创建：

`round2/`

但本轮只处理 B，因此建议：

```text
round2/
├── code/
├── logs/
├── B/
└── ROUND2_B_SUMMARY.md
```

不要预建 round3。

---

## 15. 强制停止点

完成以下内容后立即停止：

B contaminant exclusion  
→ B reclustering  
→ Round2 UMAP  
→ frozen-marker DotPlot  
→ FindAllMarkers  
→ specimen/dataset audit  
→ Round1→Round2 mapping audit  
→ annotation decision template

然后等待人工注释。

最终状态写：

`B_ROUND2_STATUS = READY_FOR_MANUAL_ANNOTATION`

如果失败：

`B_ROUND2_STATUS = PARTIAL_REVIEW_REQUIRED`
