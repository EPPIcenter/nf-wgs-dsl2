---
layout: default
title: Running QC only
nav_order: 6
has_children: false
---

# Running QC Only

To run QC only, add the --qc_only parameter
```
nextflow run main.nf \
--qc_only \
-profile docker \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory \
--trimadapter path/adapters/TruSeq2-PE.fa 
```


