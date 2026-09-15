# T/NK Round2 Harmony 敏感性分析与重新聚类方案

## 1. Round2 目标

Round2 仅针对 **T/NK lineage** 进行重新聚类与 batch-correction sensitivity analysis。

本轮与前几个 lineage 不同：

- **不删除任何 Round1 cluster**
- 当前未发现需要在 Round2 前明确剔除的 contamination cluster
- 所有 Round1 T/NK cells 全部进入 Round2
- 重点解决两个问题：
  1. Round1 中部分 cluster 存在明显 dataset-specific island
  2. resolution=1.0 下 T/NK 内部结构偏碎

因此本轮同时：

- 将 clustering resolution 从 **1.0 调整为 0.8**
- 对 Harmony 做 **4 套敏感性分析**
- 其余 preprocessing / dimensionality-reduction 参数全部保持一致
- 比较 batch mixing 与 biological structure preservation
- 选出最终推荐的 Harmony setting

本轮不自动冻结最终 subtype。

---

## 2. Round2 输入范围

Round1 T/NK 中的 **全部 cells** 进入 Round2。

不删除任何 cluster。

当前重点追踪的 Round1 structures：

### CD8 GZMK family
重点关注：
- `1`
- `3`
- `4`
- `16`

目前认为这些 cluster 很可能属于同一个 `CD8_GZMK` biological family，只是受到不同 stress / IEG / MT / epithelial-RNA state 影响而被拆开。

另外：
- `2`
- `18`

也可能与 CD8 GZMK / γδ-like cytotoxic continuum 有关，但目前与主 GZMK cluster 距离较远，并存在明显 dataset specificity，因此需重点观察。

---

### CD8 Exhausted family
重点关注：
- `0`
- `9`

两者均具有明显 exhaustion program，但 UMAP 上分离较远。

目前重点判断：

- 是否主要由 dataset effect 导致；
- 是否在更强 Harmony correction 下靠近；
- 或者 `9` 是否仍保留稳定 NK / γδ-like cytotoxic program。

---

### CD4 Helper family
重点关注：
- `7`
- `8`
- `10`
- `13`
- `14`

当前可以先看作 broad `CD4_Helper` family。

其中可能包含：
- Naive / Memory-like
- Conventional helper
- Th17-like
- CXCL13 / IL21 helper / Tfh-like

本轮不提前决定是否拆分成最终多个 subtype。

---

### Treg family
重点关注：
- `5`
- `6`
- `15`

当前均属于 Treg family。

---

### NK family
重点关注：
- `11`
- `12`

其中：
- `11` 偏 XCL1 / CD56bright-like
- `12` 偏 FCGR3A / FGFBP2 cytotoxic NK

---

## 3. Round2 clustering resolution

Round1：

`resolution = 1.0`

Round2：

**`resolution = 0.8`**

原因：

- 当前 T/NK clustering 有一定过度碎裂；
- 一些 cluster 更像状态差异，而不是独立稳定 subtype；
- 但 T/NK 内部真实异质性较高，不希望降到 0.6 后过度合并；
- 因此采用 0.8 作为中等强度的收敛。

所有 4 套 Harmony sensitivity runs 均固定：

**resolution = 0.8**

不得在不同 sensitivity run 中改变 resolution。

---

## 4. 统一 preprocessing 参数

四套 sensitivity analysis 的 preprocessing 必须完全一致：

- assay = RNA
- HVG = 3000
- PCA = 50 PCs
- downstream dims = 1:30
- Harmony sigma = 保持当前默认 / frozen setting
- Harmony max_iter = 20
- UMAP based on Harmony dims 1:30
- FindNeighbors based on Harmony dims 1:30
- FindClusters resolution = 0.8

不得在四套 runs 中改变其他参数。

---

## 5. Harmony 四套敏感性分析

### Run A — Baseline

用于作为当前方法基准。

- `group.by.vars = orig.ident`
- `theta = 2`
- `lambda = 1`
- `resolution = 0.8`

命名：

`A_baseline_origident_theta2_lambda1`

---

### Run B — Moderate correction

单纯加强 Harmony diversity penalty。

- `group.by.vars = orig.ident`
- `theta = 4`
- `lambda = 1`
- `resolution = 0.8`

命名：

`B_origident_theta4_lambda1`

---

### Run C — Strong correction

作为较强 correction sensitivity upper bound。

- `group.by.vars = orig.ident`
- `theta = 6`
- `lambda = 0.5`
- `resolution = 0.8`

命名：

`C_origident_theta6_lambda0.5`

目的：

- 判断更强 batch correction 是否进一步消除 dataset-specific islands；
- 同时检查是否出现 overcorrection；
- 不因 UMAP 更混合就自动判定为最佳结果。

---

### Run D — Dataset-level correction

专门针对当前 dataset-specific island 问题。

- `group.by.vars = dataset`
- `theta = 4`
- `lambda = 1`
- `resolution = 0.8`

命名：

`D_dataset_theta4_lambda1`

目的：

- 判断当前主要 separation 是否更接近 dataset-level effect；
- 比较 `dataset` correction 与 `orig.ident` correction 的差异；
- 避免无限加强 specimen-level correction。

---

## 6. 禁止组合校正 dataset + orig.ident

本轮不得使用：

`group.by.vars = c("dataset", "orig.ident")`

原因：

- `orig.ident` 很可能嵌套于 `dataset`
- 两者高度相关
- 双重 correction 容易造成 overcorrection
- 本轮目标是比较：
  - specimen-level correction
  - dataset-level correction

而不是叠加两者。

---

## 7. 每套 sensitivity run 必须独立输出

每套 run 都需要完整输出：

### UMAP
- cluster
- dataset
- specimen / orig.ident
- Round1 cluster

建议命名：

- `01_umap_cluster_<RUN>.pdf/png`
- `02_umap_dataset_<RUN>.pdf/png`
- `03_umap_specimen_<RUN>.pdf/png`
- `04_umap_round1_cluster_<RUN>.pdf/png`

---

### Frozen-marker DotPlot

每套 run 都用同一 frozen T/NK marker panel。

输出：

`05_dotplot_frozen_markers_<RUN>.pdf/png`

不得为不同 run 使用不同 marker panel。

---

### FindAllMarkers

每套 run 输出：

- `markers_all_<RUN>.csv`
- `markers_top10_per_cluster_<RUN>.csv`
- `markers_top20_per_cluster_<RUN>.csv`

---

### Composition audit

每套 run 输出：

- `cluster_by_dataset_counts_<RUN>.csv`
- `cluster_by_dataset_proportions_<RUN>.csv`
- `cluster_by_specimen_counts_<RUN>.csv`
- `cluster_by_specimen_proportions_<RUN>.csv`

---

### Round1 → Round2 mapping

每套 run 输出：

- `round1_to_round2_cluster_counts_<RUN>.csv`
- `round1_to_round2_cluster_row_proportions_<RUN>.csv`
- `round1_to_round2_cluster_column_proportions_<RUN>.csv`
- `round1_to_round2_cluster_heatmap_<RUN>.pdf/png`

如果方便：

- `round1_to_round2_alluvial_<RUN>.pdf/png`

---

## 8. Frozen T/NK marker framework

继续使用已有 frozen framework，不新增 marker。

### Pan-T
- `CD3D`
- `CD3E`
- `TRAC`

### Pan-NK
- `NCR1`
- `KLRD1`
- `NKG7`

### Epithelial / doublet screen
- `EPCAM`
- `KRT8`
- `KRT19`

### CD4 / CD8
- `CD4`
- `CD8A`
- `CD8B`

### Cycling
- `MKI67`
- `TOP2A`
- `UBE2C`

### CD4 Naive / Memory-like
- `CCR7`
- `TCF7`
- `LEF1`

### CD4 Helper / Conventional
- `IL7R`
- `LTB`
- `CD40LG`

### Treg
- `FOXP3`
- `IL2RA`
- `CTLA4`

### CD4 Cytotoxic
- `GZMB`
- `PRF1`
- `KLRG1`

### Tfh / CXCL13-like
- `CXCL13`
- `CXCR5`
- `ICOS`

### CD8 Naive / Memory-like
- `TCF7`
- `KLF2`
- `SELL`

### CD8 GZMK / Effector-memory-like
- `GZMK`
- `CCL5`
- `EOMES`

### CD8 Cytotoxic
- `FGFBP2`
- `GZMB`
- `CX3CR1`

### CD8 Exhausted
- `TOX`
- `ENTPD1`
- `LAG3`

### γδ T
- `TRDC`
- `TRGC1`
- `TRGC2`

### Cycling T
- `MKI67`
- `TOP2A`
- `UBE2C`

### CD56bright / XCL1-like NK
- `XCL1`
- `XCL2`
- `NCAM1`

### CD56dim / Cytotoxic NK
- `FCGR3A`
- `FGFBP2`
- `CX3CR1`

### Adaptive NK
- `KLRC2`
- `GZMH`
- `ZBTB38`

### Stressed / Dysfunctional NK
- `ATF3`
- `DUSP1`
- `GADD45B`

### Cycling NK
- `MKI67`
- `TOP2A`
- `UBE2C`

### Auxiliary tissue-resident / NKG2A-like
- `KLRC1`
- `ZNF683`
- `CXCR6`

---

## 9. 重点 biological questions

四套 sensitivity analysis 必须重点比较以下问题。

### 9.1 CD8 GZMK family

观察 Round1：

`1 / 3 / 4 / 16`

是否在新的 embedding / clustering 中：

- 更靠近；
- 合并；
- 或仍形成多个相邻 subclusters。

同时观察：

`2 / 18`

是否：

- 仍形成明显 dataset-specific island；
- 或重新靠近 CD8 GZMK family；
- 或维持 γδ-like / cytotoxic distinct program。

---

### 9.2 Exhausted CD8

重点比较：

`0 / 9`

判断：

- 更强 Harmony correction 是否使两者靠近；
- dataset-level correction 是否比 orig.ident correction 更有效；
- `9` 是否仍保留明显 NK / γδ-like cytotoxic structure。

不得仅因 UMAP 靠近就自动合并。

---

### 9.3 CD4 Helper family

观察：

`7 / 8 / 10 / 13 / 14`

是否在 resolution=0.8 后：

- 仍保持合理内部结构；
- 或出现过度合并。

本轮可暂时视为 broad `CD4_Helper` family，但不自动冻结最终 subtype。

---

### 9.4 Treg

观察：

`5 / 6 / 15`

是否保持稳定 Treg structure。

如果 strong correction 导致 Treg 被明显揉入其他 CD4 populations，视为 overcorrection warning。

---

### 9.5 NK

观察：

`11 / 12`

是否继续保持：

- XCL1 / CD56bright-like
- FCGR3A / FGFBP2 cytotoxic NK

如果 strong correction 导致两者明显丢失生物学区分，视为 overcorrection warning。

---

## 10. Harmony sensitivity 的评判标准

不得以：

> “UMAP 看起来混得最均匀”

作为唯一标准。

最终比较必须同时考虑：

### Batch mixing
- dataset segregation 是否下降
- specimen segregation 是否下降
- dataset-specific islands 是否减弱

### Biological preservation
- CD4 / CD8 / NK 大类是否仍合理
- Treg 是否稳定
- CD8 GZMK / Exhausted / Cytotoxic structure 是否保留
- NK bright / dim structure 是否保留
- γδ-like population 是否被过度吞并

### Cluster stability
- Round1 → Round2 mapping 是否具有可解释性
- 是否出现明显过度重组
- 是否有大量原 biologically coherent clusters 被强行揉在一起

---

## 11. 推荐增加定量 batch-mixing audit

如果代码实现方便，建议四套 run 都输出以下定量指标：

- dataset silhouette score
- biological cluster silhouette score
- nearest-neighbor dataset mixing
- iLISI 或等价 local inverse Simpson metric

这些用于辅助比较，不作为单独自动决策标准。

如果当前环境中实现这些指标会显著增加复杂度，可标记为 optional。

---

## 12. 最终 comparison summary

完成四套 runs 后生成：

`TNK_HARMONY_SENSITIVITY_SUMMARY.md`

至少包含：

- A / B / C / D 参数
- cell number
- cluster number
- dataset mixing summary
- specimen mixing summary
- biological structure preservation
- Round1→Round2 mapping interpretation
- overcorrection warning
- dataset-specific island warning
- 推荐主方案

并生成：

`TNK_HARMONY_SENSITIVITY_COMPARISON.csv`

建议字段：

- run
- harmony_group
- theta
- lambda
- resolution
- n_clusters
- dataset_mixing_score
- specimen_mixing_score
- biological_preservation_score
- overcorrection_warning
- notes
- recommended

---

## 13. 本轮禁止事项

Round2 禁止：

1. 删除任何 Round1 T/NK cluster
2. 将 2 / 18 直接视为 contamination
3. 将 0 / 9 自动合并
4. 将 7 / 8 / 10 / 13 / 14 自动冻结为单一最终 subtype
5. 新增 marker panel
6. 修改 frozen marker
7. 不同 Harmony runs 使用不同 resolution
8. 使用 `dataset + orig.ident` 双重 Harmony correction
9. 仅根据 UMAP 视觉效果选择最佳方案
10. 自动写入正式 `subtype_v2`
11. 启动 cell2location
12. 重跑其他 lineage

---

## 14. 输出目录

在：

`v2/round2/TNK/`

下建立：

- `A_baseline_origident_theta2_lambda1/`
- `B_origident_theta4_lambda1/`
- `C_origident_theta6_lambda0.5/`
- `D_dataset_theta4_lambda1/`
- `comparison/`

每套 run 独立保存：

- qs object
- figures
- markers
- audit
- logs

---

## 15. 强制停止点

完成以下内容后停止：

All T/NK cells  
→ common preprocessing  
→ 4 Harmony sensitivity runs  
→ resolution=0.8 clustering  
→ UMAP / DotPlot / FindAllMarkers  
→ specimen / dataset audit  
→ Round1→Round2 mapping  
→ cross-run comparison  
→ recommendation summary

然后等待人工选择最终 Harmony setting。

最终状态：

`TNK_ROUND2_SENSITIVITY_STATUS = READY_FOR_MANUAL_SELECTION`

如果失败：

`TNK_ROUND2_SENSITIVITY_STATUS = PARTIAL_REVIEW_REQUIRED`
