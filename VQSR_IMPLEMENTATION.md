# VQSR Workflow Implementation Summary

## What Was Created

### 1. Main Workflow File
**File**: `workflows/vqsr.nf`

A complete Nextflow workflow for variant quality score recalibration with the following processes:

- **genomicsdb_import**: Consolidates gVCF files across samples using GenomicsDB
- **joint_genotype**: Performs joint genotyping with GenotypeGVCFs
- **select_snps**: Extracts SNPs from joint VCF
- **select_indels**: Extracts INDELs from joint VCF
- **variant_recalibrator_snps**: Builds recalibration model for SNPs using MalariaGEN crosses
- **variant_recalibrator_indels**: Builds recalibration model for INDELs using MalariaGEN crosses
- **apply_vqsr_snps**: Applies VQSR filtering to SNPs
- **apply_vqsr_indels**: Applies VQSR filtering to INDELs
- **merge_recalibrated_variants**: Merges filtered SNPs and INDELs per chromosome
- **concat_chromosomes**: Optionally concatenates all chromosomes into single VCF

### 2. Updated Configuration Files

**File**: `nextflow.config`
Added parameters:
- `vqsr_only`: Flag to run only VQSR workflow
- `cross_vcfs_dir`: Location of MalariaGEN cross VCF files
- `vqsr_snp_filter_level`: Truth sensitivity threshold for SNPs (default: 99.0)
- `vqsr_indel_filter_level`: Truth sensitivity threshold for INDELs (default: 99.0)
- `concat_chromosomes`: Whether to merge all chromosomes (default: true)

**File**: `main.nf`
Updated to:
- Include VQSR workflow
- Support `--vqsr_only` mode
- Enable full pipeline: QC → gVCF → VQSR
- Handle channel transformations between workflows

### 3. Helper Scripts

**File**: `refs/download_malariagen_crosses.sh`
- Script to download MalariaGEN genetic cross VCF files
- Downloads 7G8×GB4, HB3×Dd2, and 3D7×HB3 crosses
- Creates proper directory structure

**File**: `run_vqsr.sh`
- Example run script for VQSR workflow
- Template with customizable parameters
- Includes reporting options

### 4. Documentation

**File**: `VQSR_README.md`
- Comprehensive documentation for VQSR workflow
- Setup instructions
- Usage examples
- Parameter descriptions
- Troubleshooting guide
- Output file descriptions

**File**: `refs/malariagen_crosses/README.md`
- Information about MalariaGEN crosses
- Download instructions
- File verification steps
- Alternative download sources

**File**: `README.md` (updated)
- Added VQSR workflow to overview
- Updated parameter table
- References to VQSR documentation

## How to Use

### Quick Start

1. **Download training data**:
   ```bash
   cd refs
   ./download_malariagen_crosses.sh
   ```

2. **Run full pipeline** (FASTQ → recalibrated VCF):
   ```bash
   nextflow run main.nf \
     --inputdir /path/to/fastq_files \
     --outputdir /path/to/results \
     -profile apptainer
   ```

3. **Run VQSR only** (if you have gVCF files):
   ```bash
   nextflow run main.nf \
     --vqsr_only \
     --inputdir /path/to/gvcf_files \
     --outputdir /path/to/results \
     -profile apptainer
   ```

### Workflow Modes

The pipeline now supports four modes:

1. **Full pipeline** (default): QC → gVCF → VQSR
2. **QC only**: `--qc_only`
3. **gVCF only**: `--gvcf_only`
4. **VQSR only**: `--vqsr_only`

## VQSR Methodology

Following the MalariaGEN Pf7 paper approach:

1. **Training Data**: Uses genetic crosses as high-confidence truth sets
   - 7G8 × GB4
   - HB3 × Dd2
   - 3D7 × HB3

2. **Variant Annotations**: Uses standard GATK annotations
   - SNPs: QD, MQ, MQRankSum, ReadPosRankSum, FS, SOR
   - INDELs: QD, MQRankSum, ReadPosRankSum, FS, SOR

3. **Recalibration Model**: Gaussian mixture model with:
   - `--max-gaussians 4`: Up to 4 clusters
   - `--trust-all-polymorphic`: For highly polymorphic P. falciparum genome
   - `prior=15.0`: High confidence in training data

4. **Filtering**: Applies truth sensitivity filter (default 99.0%)
   - Retains 99% of variants found in training data
   - Adjustable via `--vqsr_snp_filter_level` and `--vqsr_indel_filter_level`

## Output Structure

```
results/
├── genomicsdb/              # GenomicsDB workspaces per chromosome
├── joint_vcf/               # Joint-called VCFs per chromosome
├── filtered_vcf/            # Separated SNPs and INDELs
├── vqsr/                    # Recalibration models and tranches
├── recalibrated_vcf/        # VQSR-filtered variants
└── final_vcf/               # Merged variants (SNPs + INDELs)
    ├── recalibrated_chr*.vcf.gz     # Per-chromosome
    └── recalibrated_all.vcf.gz      # All chromosomes (optional)
```

## Key Features

1. **Per-chromosome processing**: Parallelized for efficiency
2. **GenomicsDB**: Efficient storage and retrieval of gVCF data
3. **Separate SNP/INDEL models**: Independent recalibration for each variant type
4. **MalariaGEN crosses**: High-quality training data specific to P. falciparum
5. **Flexible filtering**: Adjustable truth sensitivity thresholds
6. **Optional concatenation**: Combine chromosomes into single VCF

## Requirements

- **Samples**: Recommended minimum 30 samples for robust VQSR
- **Coverage**: Adequate sequencing depth (recommended 30x+)
- **Tools**: GATK 4.x, bcftools (included in container)
- **Memory**: Processes labeled with memory requirements in base.config
- **Storage**: GenomicsDB and intermediate files require substantial disk space

## Customization

### Adjust Filter Levels
```bash
nextflow run main.nf \
  --vqsr_snp_filter_level 99.5 \
  --vqsr_indel_filter_level 98.0
```

### Use Custom Training Data
Place your own VCF files in `refs/malariagen_crosses/` and they will be automatically used.

### Modify Annotations
Edit `workflows/vqsr.nf` processes `variant_recalibrator_snps` and `variant_recalibrator_indels` to change the `-an` (annotation) parameters.

### Change Max Gaussians
Edit the `--max-gaussians` parameter in recalibrator processes if you have fewer samples or want simpler models.

## Troubleshooting

### Common Issues

1. **Insufficient variants**: Need more samples or combine batches
2. **Memory errors**: Increase memory allocation in base.config
3. **Missing cross VCFs**: Run download_malariagen_crosses.sh
4. **Wrong file format**: Ensure gVCF files match naming pattern `*.chr*.g.vcf`

### Testing

For testing with small datasets:
- Use `--max-gaussians 2` for fewer samples
- Start with high filter levels (99.9) and adjust downward
- Test with single chromosome first

## Next Steps

1. Download MalariaGEN cross VCF files
2. Test workflow with a small dataset
3. Optimize memory/CPU allocation for your cluster
4. Review VQSR plots and tranches files to assess model quality
5. Adjust filter levels based on your quality requirements

## References

- **GATK Best Practices**: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112
- **MalariaGEN Pf7**: https://www.malariagen.net/data/pf7-release-7
- **Pf7 Paper**: https://doi.org/10.1038/s41586-021-03819-6
