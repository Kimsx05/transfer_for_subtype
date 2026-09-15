# Subtype v2 Round2 execution summary

Completed UTC: 2026-09-15 08:54:29 UTC

## Input, exclusion, and cluster results

- B: Round1 input 16,247; excluded 965 (clusters 12/13/14/16/17); Round2 15,282 cells; resolution 0.6; 12 clusters.
- Endothelial: Round1 input 17,039; excluded 1,843 (clusters 7/11/18/20/22); Round2 15,196 cells; resolution 0.6; 16 clusters.
- Fibroblast: Round1 input 27,615; excluded 3,166 (clusters 10/11/12/19/21/23); Round2 24,449 cells; resolution 1.0; 20 clusters.
- Myeloid: Round1 input 40,393; excluded 4,309 (clusters 5/11/19/23); Round2 36,084 cells; resolution 0.8; 21 clusters.
- TNK: Round1/Round2 input 124,746; excluded 0; resolution 0.8; A=16, B=14, C=14, D=29 clusters.

## Composition warnings

- B: dataset >70% clusters 0/4/6; cluster 0 >90% InhouseData. No specimen >50% cluster.
- Endothelial: dataset >70% clusters 0/1/2/3/4/5/6/7/10/11/12/13/14/15; specimen >50% cluster 2; dataset >90% clusters 2/4/5/10.
- Fibroblast: dataset >70% clusters 0/3/4/6/8/10/14/19; specimen >50% cluster 19; dataset >90% clusters 10/19. Cluster 19 has 97 cells, six specimens, and 100% InhouseData.
- Myeloid: dataset >70% clusters 4/6/10/14/20; dataset >90% clusters 4/6. No specimen >50% cluster.
- TNK A: dataset >90% clusters 6/7/9/13; no specimen >50% cluster.
- TNK B: cluster 13 is both specimen >50% and dataset >90%.
- TNK C: cluster 13 is both specimen >50% and dataset >90%.
- TNK D: specimen >50% clusters 15/16/18/19/20/21/22/23/24/25/26/27/28; dataset >70% clusters 15/16/17/18/19/21/22/23/24/25/26/27/28; 11 are >90% dataset.

## Main Round1 to Round2 reorganizations

- B: Round1 2+7 converge into Round2 0; Round1 3+4 into Round2 2; Round1 0+5 into Round2 3; Round1 10 and 11 remain largely separate as Round2 9 and 10; Round1 15 mainly maps to Round2 11.
- Endothelial: Round1 0+3 dominate Round2 0; 2+8 dominate Round2 1; 4+13 dominate Round2 2; 5+2+10 contribute to Round2 3; Round1 14 remains almost entirely Round2 8; Round1 12 and 16 remain highly concentrated in Round2 6 and 12; several smaller states remain split.
- Fibroblast: Round1 3 remains concentrated in Round2 4; Round1 5 concentrates in Round2 2; Round1 9 in Round2 7; Round1 20 in Round2 17; Round1 18 splits mainly across Round2 14/18; Round1 9/15/22 do not collapse into one cluster.
- Myeloid: Round1 1/2/21 primarily map to Round2 0/2/0; Round1 7 and 15 remain primarily separated (Round2 9 and 5); Round1 0 mainly maps to Round2 1; Round1 10 remains highly concentrated in Round2 6, which is 98.2% InhouseData; Round1 13/22 distribute across Round2 17/1 and other TAM-associated clusters.

## TNK Harmony sensitivity conclusion

- A: dataset mixing 0.4486; specimen mixing 0.8628; Round1-family preservation 0.8704; NMI 0.6813; 16 clusters.
- B: stronger mixing (0.6049/0.9037) but preservation and NMI fall to 0.7722/0.5768; 14 clusters; overcorrection warning.
- C: strongest mixing (0.6161/0.9066) but further preservation/NMI loss to 0.7508/0.5561; 14 clusters; overcorrection warning, including stronger Treg redistribution.
- D: 29 clusters, weak specimen mixing 0.5390, 13 specimen-driven clusters, 11 >90% dataset clusters, family preservation 0.7954 and NMI 0.4759; strong fragmentation/overcorrection warning.
- Recommended for manual selection: A_baseline_origident_theta2_lambda1. It best preserves Round1 biological-family neighborhoods and mapping coherence; B/C improve batch mixing but cross the predefined preservation-loss warning, while D markedly fragments the structure and retains strong sample/dataset islands. A still requires review of its four dataset-specific clusters and is not accepted as a final subtype annotation.

## Integrity and stopping point

All eight final objects were re-read successfully. RNA counts/data, Round1 and Round2 cluster metadata, shared/new PCA, Harmony, and UMAP are present. `subtype_v2` is absent. No automatic annotation, merge, deletion from the master object, or cell2location was performed.

B_ROUND2_STATUS = READY_FOR_MANUAL_ANNOTATION
ENDO_ROUND2_STATUS = READY_FOR_MANUAL_ANNOTATION
FIBRO_ROUND2_STATUS = READY_FOR_MANUAL_ANNOTATION
MYELOID_ROUND2_STATUS = READY_FOR_MANUAL_ANNOTATION
TNK_ROUND2_SENSITIVITY_STATUS = READY_FOR_MANUAL_SELECTION
