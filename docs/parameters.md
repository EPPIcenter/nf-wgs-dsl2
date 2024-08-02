---
layout: default
title: Pipeline Parameters
nav_order: 3
has_children: false
---

# Modifying Parameters

Parameters can be found in `nextflow.config` 


|Parameters|Description|
|---|---|
|qc_only|If enabled, only QC workflow is run (default 'false')|
|gvcf_only|If enabled, only gVCF workflow is run (default 'false')|
|inputdir|The folder that contains the input files (default 'data')|
|outputdir|The folder where you want the resulting data to be save (default 'results/results')|
|trimadapter|The adapter used for initial trimming of reads (default 'NexteraPE-custom.fa')|


The nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on.





Below is an example of how you may use the above parameters on the command line:

```bash
nextflow run main.nf --readDIR /wynton/scratch/data --outDIR /wynton/scratch/results -profile sge,apptainer --target v4 -config conf/custom.config --omega_a 1e-120 --band_size 16 --pool pseudo
```






