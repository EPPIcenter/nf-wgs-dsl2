---
layout: default
title: Running GVCF only
nav_order: 8
has_children: false
---

# Running GVCF Only

To run GVCF only, simply add the --gvcf_only parameter

Important note: if the pipeline run starts from the GVCF workflow, the `--inputdir` must point to the `.sorted.dup.pf.bam` files directory, which is one of the final outputs if the QC workflow. 

```
nextflow run main.nf \
-profile sge,apptainer \
--inputdir path/input_directory_bam \
--outputdir path/output_directory \
--gvcf_only true
```



