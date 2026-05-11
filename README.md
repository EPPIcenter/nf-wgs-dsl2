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
  - `filter.nf`: Joint genotyping and hard-filtering sub-workflow
  - `vqsr.nf`: VQSR sub-workflow (**not yet functional** :construction:)
- `conf`
  - `Apptainer`: file used to build nf-wgs-dsl2.sif  
  - `Dockerfile`: file for building docker image 
  - `base.config`: base config file 
  - `envs`: conda envs (under construction :construction:)
- `refs`: reference files used by both `QC_workflow` and `gVCF_workflow`
  - `adapters`: folder containing trimmomatic adapter files
  - `genomes`: reference genome files and more
  - `malariagen_crosses`: MalariaGEN genetic cross VCF files (for future VQSR use)
  - `run_quality_report.Rmd`: r script for quality report used in `QC_workflow`
  - `download_malariagen_crosses.sh`: script to download MalariaGEN cross data
- *`data`: suggested directory for input files*
- *`results`: suggested directory for output*

## Workflows

This pipeline includes three active workflows:

1. **QC Workflow** (`qc.nf`): Quality control, read trimming, alignment, and BAM processing
2. **gVCF Workflow** (`gvcf.nf`): Per-sample variant calling to generate gVCF files
3. **Filter Workflow** (`filter.nf`): Joint genotyping (GenomicsDBImport + GenotypeGVCFs) followed by hard filtering

> **Note:** A VQSR workflow (`vqsr.nf`) is included in the repository but is **not yet functional**.

## Parameters

### nextflow.config
|Parameter|Description|
|---|---|
|`qc_only`|If enabled, only the QC workflow is run (default: `false`)|
|`gvcf_only`|If enabled, only the gVCF workflow is run (default: `false`)|
|`filter_only`|If enabled, only the filter workflow is run (default: `false`)|
|`inputdir`|Folder containing input files (default: `data`)|
|`outputdir`|Folder where results are saved (default: `results`)|
|`trim_adapter`|Adapter file for Trimmomatic (default: `refs/adapters/TruSeq3-PE.fa`)|
|`genomes_dir`|Directory containing reference genome files (default: `refs/genomes`)|
|`sif_path`|Path to the Apptainer/Singularity image (default: `conf/nf-wgs-dsl2.sif`)|
|`ploidy`|Organism ploidy for HaplotypeCaller and GenotypeGVCFs (default: `1`)|
|`output_zarr`|Convert final VCF to zarr format (default: `false`)|
|`filter_QD`|QualByDepth hard filter threshold (default: `2.0`)|
|`filter_FS`|FisherStrand hard filter threshold (default: `60.0`)|
|`filter_MQ`|RMSMappingQuality hard filter threshold (default: `30.0`)|
|`filter_MQRankSum`|MappingQualityRankSumTest hard filter threshold (default: `-12.5`)|
|`filter_ReadPosRankSum`|ReadPosRankSumTest hard filter threshold (default: `-8.0`)|
|`filter_SOR`|StrandOddsRatio hard filter threshold (default: `3.0`)|
|`filter_GQ`|Genotype quality genotype filter threshold (default: `10`)|
|`filter_DP`|Depth genotype filter threshold (default: `5`)|
|`filter_AD`|Allele depth genotype filter threshold (default: `2`)|

