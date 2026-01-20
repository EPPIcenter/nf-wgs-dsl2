# WGS Nextflow Pipeline 
## (nf-wgs-dsl2)

Adapted from: 
- https://github.com/Karaniare/Optimized_GATK4_pipeline (shell script)
- https://github.com/jhoneycuttr/nf-wgs (Nextflow DSL 1)

**Documentation, please refer to: https://eppicenter.github.io/nf-wgs-dsl2/**

## Overview
- `main.nf`: WGS workflow 
- `nextflow.config`: config file
- `workflows` 
  - `qc.nf`: QC sub-workflow 
  - `gvcf.nf`: GVCF sub-workflow
  - `vqsr.nf`: Variant Quality Score Recalibration (VQSR) sub-workflow
- `config`
  - `Apptainer`: file used to build nf-wgs-dsl2.sif  
  - `Dockerfile`: file for building docker image 
  - `base.config`: base config file 
  - `envs`: conda envs (under construction :construction:)
- `refs`: reference files used by both `QC_workflow` and `gVCF_workflow`
  - `adapters`: folder containing trimmomatic adapter files
  - `genomes`: reference genome files and more
  - `malariagen_crosses`: MalariaGEN genetic cross VCF files for VQSR training
  - `run_quality_report.Rmd`: r script for quality report used in `QC_workflow`
  - `download_malariagen_crosses.sh`: script to download MalariaGEN cross data
- *`data`: suggested directory for input files*
- *`results`: suggested directory for output*

## Workflows

This pipeline includes three main workflows:

1. **QC Workflow** (`qc.nf`): Quality control, read trimming, alignment, and BAM processing
2. **gVCF Workflow** (`gvcf.nf`): Per-sample variant calling to generate gVCF files
3. **VQSR Workflow** (`vqsr.nf`): Joint genotyping and variant quality score recalibration

See [VQSR_README.md](VQSR_README.md) for detailed documentation on the VQSR workflow.

## Parameters

### nextflow.config
|Parameters|Description|
|---|---|
|qc_only|If enabled, only QC workflow is run (default 'false')|
|gvcf_only|If enabled, only gVCF workflow is run (default 'false')|
|vqsr_only|If enabled, only VQSR workflow is run (default 'false')|
|inputdir|The folder that contains the input files (default 'data')|
|outputdir|The folder where you want the resulting data to be save (default 'results/results')|
|trimadapter|The adapter used for initial trimming of reads (default 'NexteraPE-custom.fa')|
|cross_vcfs_dir|Directory containing MalariaGEN cross VCF files for VQSR (default 'refs/malariagen_crosses')|
|vqsr_snp_filter_level|Truth sensitivity filter level for SNPs in VQSR (default '99.0')|
|vqsr_indel_filter_level|Truth sensitivity filter level for INDELs in VQSR (default '99.0')|
|concat_chromosomes|Concatenate all chromosomes into single VCF (default 'true')|

