# Plasmodium Falciparum WGS Pipeline 
## (Nextflow DSL 2)

Adapted from: 
- https://github.com/Karaniare/Optimized_GATK4_pipeline (shell script)
- https://github.com/jhoneycuttr/nf-wgs (Nextflow DSL 1)

Documentation, please refer to: 
https://eppicenter.github.io/nf-wgs-dsl2/

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

