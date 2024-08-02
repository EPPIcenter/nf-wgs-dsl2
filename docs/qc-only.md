---
layout: default
title: Running QC only
nav_order: 7
has_children: false
---

# Running QC Only

To run QC only, add the --qc_only parameter
```
nextflow run main.nf \
-profile sge,apptainer \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory \
--qc_only true
```

