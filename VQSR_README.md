# Variant Quality Score Recalibration (VQSR) Workflow

This workflow performs joint genotyping and variant quality score recalibration using MalariaGEN genetic crosses as high-confidence training data, following the methodology from the Pf7 paper.

## Overview

The VQSR workflow includes the following steps:

1. **GenomicsDB Import**: Consolidate gVCF files across all samples per chromosome
2. **Joint Genotyping**: Call variants across all samples simultaneously using GenotypeGVCFs
3. **Variant Selection**: Separate SNPs and INDELs for independent recalibration
4. **Variant Recalibration**: Build recalibration models using MalariaGEN crosses as truth sets
5. **Apply VQSR**: Filter variants based on VQSR scores
6. **Merge Variants**: Combine recalibrated SNPs and INDELs
7. **Concatenate (optional)**: Combine all chromosomes into a single VCF file

## Setup

### 1. Download MalariaGEN Cross VCF Files

First, download the genetic cross VCF files that will be used as training data:

```bash
cd refs
chmod +x download_malariagen_crosses.sh
./download_malariagen_crosses.sh
```

This will download three genetic crosses:
- 7G8 × GB4
- HB3 × Dd2
- 3D7 × HB3

These crosses contain high-confidence variant calls and are used as truth sets for the recalibration model.

### 2. Verify Directory Structure

Ensure you have the following directory structure:

```
refs/
  malariagen_crosses/
    7G8_GB4.vcf.gz
    7G8_GB4.vcf.gz.tbi
    HB3_Dd2.vcf.gz
    HB3_Dd2.vcf.gz.tbi
    3D7_HB3.vcf.gz
    3D7_HB3.vcf.gz.tbi
  genomes/
    Pf3D7.fasta
    (and other reference files)
```

## Usage

### Run Full Pipeline (QC → gVCF → VQSR)

Run the complete pipeline from FASTQ files through variant calling to recalibration:

```bash
nextflow run main.nf \
  --inputdir /path/to/fastq_files \
  --outputdir /path/to/results \
  -profile apptainer
```

### Run VQSR Only

If you already have gVCF files and want to run only the VQSR workflow:

```bash
nextflow run main.nf \
  --vqsr_only \
  --inputdir /path/to/gvcf_files \
  --outputdir /path/to/results \
  -profile apptainer
```

**Note**: The input directory should contain gVCF files with naming pattern: `SAMPLE.chr*.g.vcf` and `SAMPLE.chr*.g.vcf.idx`

### Run QC and gVCF Only (Skip VQSR)

```bash
nextflow run main.nf \
  --gvcf_only \
  --inputdir /path/to/bam_files \
  --outputdir /path/to/results \
  -profile apptainer
```

## Parameters

Key parameters that can be modified in `nextflow.config` or via command line:

### Required Parameters
- `--inputdir`: Directory containing input files (FASTQ, BAM, or gVCF depending on workflow mode)
- `--outputdir`: Directory for output files

### VQSR-Specific Parameters
- `--cross_vcfs_dir`: Directory containing MalariaGEN cross VCF files (default: `$projectDir/refs/malariagen_crosses`)
- `--vqsr_snp_filter_level`: Truth sensitivity filter level for SNPs (default: 99.0)
  - Recommended range: 90.0-99.9
  - Higher values = more sensitive but less specific
- `--vqsr_indel_filter_level`: Truth sensitivity filter level for INDELs (default: 99.0)
- `--concat_chromosomes`: Concatenate all chromosomes into single VCF (default: true)

### Workflow Mode Parameters
- `--qc_only`: Run only QC workflow (default: false)
- `--gvcf_only`: Run only gVCF calling workflow (default: false)
- `--vqsr_only`: Run only VQSR workflow (default: false)

## Example Command with Custom Parameters

```bash
nextflow run main.nf \
  --inputdir /path/to/gvcf_files \
  --outputdir /path/to/results \
  --vqsr_only \
  --vqsr_snp_filter_level 99.5 \
  --vqsr_indel_filter_level 98.0 \
  --concat_chromosomes true \
  -profile apptainer \
  -resume
```

## Output Files

The VQSR workflow generates the following outputs in `outputdir`:

### GenomicsDB
- `genomicsdb/genomicsdb_chr*/`: GenomicsDB workspaces per chromosome

### Joint Genotyping
- `joint_vcf/joint_chr*.vcf.gz`: Joint-called VCF files per chromosome
- `joint_vcf/joint_chr*.vcf.gz.tbi`: Index files

### Filtered Variants
- `filtered_vcf/snps_chr*.vcf.gz`: SNPs only per chromosome
- `filtered_vcf/indels_chr*.vcf.gz`: INDELs only per chromosome

### VQSR Results
- `vqsr/snps_chr*.recal`: SNP recalibration model
- `vqsr/snps_chr*.tranches`: SNP tranches file
- `vqsr/snps_chr*.plots.R`: R script for SNP recalibration plots
- `vqsr/indels_chr*.recal`: INDEL recalibration model
- `vqsr/indels_chr*.tranches`: INDEL tranches file
- `vqsr/indels_chr*.plots.R`: R script for INDEL recalibration plots

### Recalibrated Variants
- `recalibrated_vcf/snps_recal_chr*.vcf.gz`: Recalibrated SNPs per chromosome
- `recalibrated_vcf/indels_recal_chr*.vcf.gz`: Recalibrated INDELs per chromosome

### Final Output
- `final_vcf/recalibrated_chr*.vcf.gz`: Merged SNPs and INDELs per chromosome
- `final_vcf/recalibrated_all.vcf.gz`: All chromosomes concatenated (if `concat_chromosomes=true`)

## VQSR Annotations

The following annotations are used for building the recalibration model:

### SNPs
- `QD`: Quality by Depth
- `MQ`: RMS Mapping Quality
- `MQRankSum`: Mapping Quality Rank Sum Test
- `ReadPosRankSum`: Read Position Rank Sum Test
- `FS`: Fisher Strand Bias
- `SOR`: Strand Odds Ratio

### INDELs
- `QD`: Quality by Depth
- `MQRankSum`: Mapping Quality Rank Sum Test
- `ReadPosRankSum`: Read Position Rank Sum Test
- `FS`: Fisher Strand Bias
- `SOR`: Strand Odds Ratio

## Understanding VQSR Filter Levels

The `vqsr_snp_filter_level` and `vqsr_indel_filter_level` parameters control the truth sensitivity cutoff:

- **99.0** (default): Retains 99% of true variants, good balance
- **99.5**: More sensitive, retains more variants but may include more false positives
- **99.9**: Very sensitive, use for research where recall is critical
- **90.0**: More stringent, higher confidence but may lose true variants

## Troubleshooting

### Not Enough Variants for VQSR
If you get errors about insufficient variants for building the model:
- Ensure you have enough samples (recommended: 30+)
- Try lowering `--max-gaussians` in the VariantRecalibrator steps
- Consider combining multiple batches of samples

### Memory Issues
If processes fail with out-of-memory errors:
- Check `conf/base.config` and increase memory for `big_mem` label
- Reduce the number of samples processed at once
- Use `--batch-size` parameter in GenomicsDBImport

### Missing Cross VCF Files
If the pipeline fails to find cross VCF files:
- Verify files are in `refs/malariagen_crosses/`
- Check file names match the expected pattern
- Ensure index files (.tbi) are present

## References

1. **MalariaGEN Pf7 Paper**:
   - Genomic epidemiology of artemisinin resistant malaria
   - Nature (2021) https://doi.org/10.1038/s41586-021-03819-6

2. **GATK Best Practices**:
   - https://gatk.broadinstitute.org/hc/en-us/articles/360035531112

3. **MalariaGEN Data**:
   - https://www.malariagen.net/data/pf7-release-7

## Notes

- The VQSR workflow is designed for multi-sample cohorts (minimum 30 samples recommended)
- For small sample sizes (<30), consider hard filtering instead of VQSR
- The genetic crosses provide high-quality training data specific to *P. falciparum*
- Recalibration models are built independently for SNPs and INDELs
- The `--trust-all-polymorphic` flag is used as crosses are highly polymorphic
