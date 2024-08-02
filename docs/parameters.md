---
layout: default
title: Pipeline Parameters
nav_order: 3
has_children: false
---

# Modifying Parameters

## Overview
- `main.nf`: WGS workflow 
- `nextflow.config`: config file
- `workflows` 
  - `qc.nf`: QC sub-workflow 
  - `gvcf.nf`: GVCF sub-workflow
- `config`
  - `Apptainer`: file used to build nf-wgs-dsl2.sif  
  - `Dockerfile`: file for building docker image 
  - `base.config`: base config file 
  - `envs`: conda envs (under construction :construction:)
- `refs`: reference files used by both `QC_workflow` and `gVCF_workflow`
  - `adapters`: folder containing trimmomatic adapter files
  - `genomes`: reference genome files and more
  - `run_quality_report.Rmd`: r script for quality report used in `QC_workflow`
- *`data`: suggested directory for input files*
- *`results`: suggested directory for output*

## Parameters

### nextflow.config
|Parameters|Description|
|---|---|
|qc_only|If enabled, only QC workflow is run (default 'false')|
|gvcf_only|If enabled, only gVCF workflow is run (default 'false')|
|inputdir|The folder that contains the input files (default 'data')|
|outputdir|The folder where you want the resulting data to be save (default 'results/results')|
|trimadapter|The adapter used for initial trimming of reads (default 'NexteraPE-custom.fa')|

The nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on.

#### Parameters in main.nf
|Parameters|Description|
|---|---|
|refdir|Reference directory is assumed to be located here '$projectDir/refs/genomes'|
|rscript|Rscript for run quality report is assumed to be located here "$projectDir/refs/run_quality_report.Rmd'|
|reads|If running the pipeline starting from QC, inputdir is assumed to contain the raw reads fastq.gz files|
|bams|If running the pipeline starting from gVCF, inputdir is assumed to contain the pf bam files and .csi index files|




Below is an example of how you may use the above parameters on the command line:

```bash
nextflow run main.nf --readDIR /wynton/scratch/data --outDIR /wynton/scratch/results -profile sge,apptainer --target v4 -config conf/custom.config --omega_a 1e-120 --band_size 16 --pool pseudo
```






