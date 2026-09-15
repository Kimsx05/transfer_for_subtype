# Round2 execution manifest

Generated UTC: 2026-09-15 01:25:51 UTC

## Frozen MD files

- `B_round2_reclustering_plan_res0.6.md` (md5: `c1bca0182a5807cc399b414c9b123d6b` )
- `Endothelial_round2_reclustering_plan_res0.6.md` (md5: `44dd41462d3c64411385d28e216a2d2f` )
- `Fibroblast_round2_reclustering_plan_res1.0.md` (md5: `00d134bf5f408374fac678def7c9fd64` )
- `Myeloid_round2_reclustering_plan_res0.8.md` (md5: `db75879a58682eb3cbce54e57b18f373` )
- `TNK_round2_harmony_sensitivity_res0.8.md` (md5: `51a867707e806885be64f3f7d3709c5f` )

## Execution matrix

- B: Round1 input `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round1/B/B_round1_clustered.qs`; input=16247; exclusion=[12,13,14,16,17]; excluded cells=965; Round2 input=15282; resolution=0.6; Harmony=orig.ident; theta/default; lambda/default; marker changes=none; output=`/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round2/B`
- Endothelial: Round1 input `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round1/Endothelial/Endothelial_round1_clustered.qs`; input=17039; exclusion=[7,11,18,20,22]; excluded cells=1843; Round2 input=15196; resolution=0.6; Harmony=orig.ident; theta/default; lambda/default; marker changes=none; output=`/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round2/Endothelial`
- Fibroblast: Round1 input `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round1/Fibroblast/Fibroblast_round1_clustered.qs`; input=27615; exclusion=[10,11,12,19,21,23]; excluded cells=3166; Round2 input=24449; resolution=1; Harmony=orig.ident; theta/default; lambda/default; marker changes=append only MyoFibro / SMC: ACTA2,TAGLN,MYH11; output=`/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round2/Fibroblast`
- Myeloid: Round1 input `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round1/Myeloid/Myeloid_round1_clustered.qs`; input=40393; exclusion=[5,11,19,23]; excluded cells=4309; Round2 input=36084; resolution=0.8; Harmony=orig.ident; theta/default; lambda/default; marker changes=none; output=`/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round2/Myeloid`
- TNK: Round1 input `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round1/TNK/TNK_round1_clustered.qs`; input=124746; exclusion=[]; excluded cells=0; Round2 input=124746; resolution=0.8; Harmony=A orig.ident theta2 lambda1; B orig.ident theta4 lambda1; C orig.ident theta6 lambda0.5; D dataset theta4 lambda1; max_iter_harmony=20; marker changes=none; output=`/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/round2/TNK`

## Conflict and completeness audit

- No substantive conflicts or missing required inputs detected.
- Round1 objects are read-only inputs; Round1 cluster metadata will be retained.
- Round2 exclusions apply only to lineage working copies.
- Stopping point: after reclustering, visualization, markers, composition, mapping, and TNK sensitivity comparison; before annotation/cell2location.
