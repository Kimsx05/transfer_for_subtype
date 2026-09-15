# TNK Harmony sensitivity summary

All four runs used the identical 124,746-cell set, 3,000-HVG logic, shared 50-PC PCA, dimensions 1:30, UMAP settings, resolution 0.8, seed, and frozen marker panel. No combined `dataset + orig.ident` correction was used.

## Quantitative auxiliary comparison

- A_baseline_origident_theta2_lambda1: clusters=16; dataset neighbor mixing=0.4486; specimen neighbor mixing=0.8628; Round1-family preservation=0.8704; Round1→Round2 NMI=0.6813; >90% dataset clusters=4; overcorrection warning=FALSE
- B_origident_theta4_lambda1: clusters=14; dataset neighbor mixing=0.6049; specimen neighbor mixing=0.9037; Round1-family preservation=0.7722; Round1→Round2 NMI=0.5768; >90% dataset clusters=1; overcorrection warning=TRUE
- C_origident_theta6_lambda0.5: clusters=14; dataset neighbor mixing=0.6161; specimen neighbor mixing=0.9066; Round1-family preservation=0.7508; Round1→Round2 NMI=0.5561; >90% dataset clusters=1; overcorrection warning=TRUE
- D_dataset_theta4_lambda1: clusters=29; dataset neighbor mixing=0.4073; specimen neighbor mixing=0.5390; Round1-family preservation=0.7954; Round1→Round2 NMI=0.4759; >90% dataset clusters=11; overcorrection warning=TRUE

## Biological preservation audit

- A_baseline_origident_theta2_lambda1 / CD8_GZMK_core: dominant Round2 cluster 2 (0.319); clusters covering 80% = 3
- A_baseline_origident_theta2_lambda1 / CD8_GZMK_island_review: dominant Round2 cluster 7 (0.654); clusters covering 80% = 3
- A_baseline_origident_theta2_lambda1 / Exhausted_CD8_review: dominant Round2 cluster 0 (0.64); clusters covering 80% = 2
- A_baseline_origident_theta2_lambda1 / CD4_Helper_review: dominant Round2 cluster 5 (0.242); clusters covering 80% = 4
- A_baseline_origident_theta2_lambda1 / Treg_review: dominant Round2 cluster 1 (0.642); clusters covering 80% = 2
- A_baseline_origident_theta2_lambda1 / NK_review: dominant Round2 cluster 8 (0.575); clusters covering 80% = 2
- B_origident_theta4_lambda1 / CD8_GZMK_core: dominant Round2 cluster 1 (0.327); clusters covering 80% = 3
- B_origident_theta4_lambda1 / CD8_GZMK_island_review: dominant Round2 cluster 1 (0.489); clusters covering 80% = 3
- B_origident_theta4_lambda1 / Exhausted_CD8_review: dominant Round2 cluster 0 (0.864); clusters covering 80% = 1
- B_origident_theta4_lambda1 / CD4_Helper_review: dominant Round2 cluster 4 (0.34); clusters covering 80% = 3
- B_origident_theta4_lambda1 / Treg_review: dominant Round2 cluster 2 (0.746); clusters covering 80% = 2
- B_origident_theta4_lambda1 / NK_review: dominant Round2 cluster 7 (0.578); clusters covering 80% = 2
- C_origident_theta6_lambda0.5 / CD8_GZMK_core: dominant Round2 cluster 1 (0.325); clusters covering 80% = 3
- C_origident_theta6_lambda0.5 / CD8_GZMK_island_review: dominant Round2 cluster 1 (0.485); clusters covering 80% = 3
- C_origident_theta6_lambda0.5 / Exhausted_CD8_review: dominant Round2 cluster 0 (0.832); clusters covering 80% = 1
- C_origident_theta6_lambda0.5 / CD4_Helper_review: dominant Round2 cluster 4 (0.314); clusters covering 80% = 4
- C_origident_theta6_lambda0.5 / Treg_review: dominant Round2 cluster 6 (0.503); clusters covering 80% = 2
- C_origident_theta6_lambda0.5 / NK_review: dominant Round2 cluster 8 (0.566); clusters covering 80% = 2
- D_dataset_theta4_lambda1 / CD8_GZMK_core: dominant Round2 cluster 1 (0.225); clusters covering 80% = 6
- D_dataset_theta4_lambda1 / CD8_GZMK_island_review: dominant Round2 cluster 2 (0.256); clusters covering 80% = 8
- D_dataset_theta4_lambda1 / Exhausted_CD8_review: dominant Round2 cluster 3 (0.338); clusters covering 80% = 5
- D_dataset_theta4_lambda1 / CD4_Helper_review: dominant Round2 cluster 7 (0.207); clusters covering 80% = 6
- D_dataset_theta4_lambda1 / Treg_review: dominant Round2 cluster 0 (0.729); clusters covering 80% = 2
- D_dataset_theta4_lambda1 / NK_review: dominant Round2 cluster 10 (0.497); clusters covering 80% = 2

## Recommendation

Recommended for manual selection: `A_baseline_origident_theta2_lambda1`. This recommendation maximizes a balanced auxiliary score combining dataset/specimen neighbor mixing with Round1-family neighborhood preservation and Round1→Round2 NMI, while excluding runs flagged for overcorrection. It is not based on UMAP appearance alone and is not a final subtype annotation.

Manual review must still prioritize CD8 GZMK 1/3/4/16 and islands 2/18, Exhausted 0/9, Helper 7/8/10/13/14, Treg 5/6/15, and NK 11/12 using the paired comparison figures, mapping tables, DotPlots, and markers.

TNK_ROUND2_SENSITIVITY_STATUS = READY_FOR_MANUAL_SELECTION
