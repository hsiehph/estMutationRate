# estMutationRate
This is a simple pipeline to estimate rate of substitution using pairwise sequence alignments

Usage:
  Simply clone this repo, edit the config.yaml file, and run the snakemake.

# Config.yaml:

---
<b> input dir that contains all sequences in the FASTA format </b>
dir_inFASTA: ./inputFASTA
<br> list1 must contain all possible IDs </br>
list1_ids: ["chr5", "chm1_cen5v2", "clint_PTR_cen5_h1v1","clint_PTR_cen5_h2v1","susie_pab_cen5_h2v1","ag07107_mmu_cen5_h1v1","ag07107_mmu_cen5_h2v1"]
<br> list2 contains only human sample IDs </br>
list2_ids: ["chr5", "chm1_cen5v2"]
