# estMutationRate
This is a simple pipeline to estimate rate of substitution using pairwise sequence alignments

Usage:
  Simply clone this repo, edit the config.yaml file, and run the snakemake.

# Config.yaml:

---
<b> #input dir that contains all sequences in the FASTA format </b>
<p>dir_inFASTA: ./inputFASTA</p>

<b> #list1 must contain all possible IDs </b>
<p>list1_ids: ["chr5", "chm1_cen5v2", "clint_PTR_cen5_h1v1","clint_PTR_cen5_h2v1","susie_pab_cen5_h2v1","ag07107_mmu_cen5_h1v1","ag07107_mmu_cen5_h2v1"]</p>

<b> #list2 contains only human sample IDs </b>
<p>list2_ids: ["chr5", "chm1_cen5v2"]</p>
