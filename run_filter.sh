#!/bin/bash
#SBATCH --job-name=filter_variants
#SBATCH --output=logs/filter_%j.out
#SBATCH --error=logs/filter_%j.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=ncpu

# Variant Filtering Workflow for Lab Mutants
# Date: 2026-02-02

# Set Nextflow options
export NXF_OPTS='-Xms1g -Xmx4g'

# Load required modules
module purge
module load Singularity/3.11.3
module load Python/3.11.5-GCCcore-13.2.0
module load Nextflow/25.04.4
ml git 
ssh-add ~/.ssh/id_ed25519

# Set paths
GVCF_DIR="/nemo/stp/babs/working/whittog/projects/saterialea/christine.collins/cc1026/analysis/run_variant_calling/results"
OUTPUT_DIR="/nemo/stp/babs/working/whittog/projects/saterialea/christine.collins/cc1026/analysis/run_variant_calling/results_filtered"
WORK_DIR="/nemo/stp/babs/working/whittog/projects/saterialea/christine.collins/cc1026/analysis/run_variant_calling/work_filtered"

# Create output and work directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${WORK_DIR}
mkdir -p logs

# Prepare input directory with symlinks organized by sample and chromosome
INPUT_DIR="${OUTPUT_DIR}/gvcf_input"
rm -rf ${INPUT_DIR}  # Clean any previous run
mkdir -p ${INPUT_DIR}

echo "================================================================"
echo "Variant Filtering Workflow - Lab Mutants"
echo "Date: $(date)"
echo "================================================================"
echo ""
echo "Preparing gVCF input files..."
echo ""

# Get list of unique samples from chr1 directory
cd ${GVCF_DIR}/chr1
SAMPLE_COUNT=0

for GVCF_FILE in *.g.vcf; do
    SAMPLE_NAME=$(basename ${GVCF_FILE} .chr1.g.vcf)
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    
    echo "Sample ${SAMPLE_COUNT}: ${SAMPLE_NAME}"
    
    # Link all chromosomes for this sample
    for chr in {1..14}; do
        GVCF="${GVCF_DIR}/chr${chr}/${SAMPLE_NAME}.chr${chr}.g.vcf"
        IDX="${GVCF}.idx"
        
        if [ -f "${GVCF}" ]; then
            ln -sf "${GVCF}" "${INPUT_DIR}/${SAMPLE_NAME}.chr${chr}.g.vcf"
            ln -sf "${IDX}" "${INPUT_DIR}/${SAMPLE_NAME}.chr${chr}.g.vcf.idx"
        else
            echo "  WARNING: Missing ${GVCF}"
        fi
    done
done

echo ""
echo "================================================================"
echo "Total samples: ${SAMPLE_COUNT}"
echo "Total gVCF files: $(ls ${INPUT_DIR}/*.g.vcf 2>/dev/null | wc -l)"
echo "Expected: $((SAMPLE_COUNT * 14))"
echo "================================================================"
echo ""
echo "Running filtering workflow..."
echo "  - Joint genotyping across all samples"
echo "  - Concatenating all chromosomes"
echo "  - Hard filtering (QD, FS, MQ, MQRankSum, ReadPosRankSum, SOR)"
echo "  - Converting to zarr format"
echo ""

# Run filtering workflow
nextflow run /nemo/stp/babs/working/whittog/pipelines/nf-wgs-pf-dsl2/main.nf \
  -profile slurm,singularity \
  --filter_only \
  --inputdir ${INPUT_DIR} \
  --outputdir ${OUTPUT_DIR} \
  --output_zarr true \
  -work-dir ${WORK_DIR} \
  -resume \
  -with-report ${OUTPUT_DIR}/filter_report.html \
  -with-timeline ${OUTPUT_DIR}/filter_timeline.html \
  -with-trace ${OUTPUT_DIR}/filter_trace.txt

EXIT_CODE=$?

echo ""
echo "================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Filtering workflow completed successfully!"
    echo ""
    echo "Output files:"
    echo "  - Joint VCF: ${OUTPUT_DIR}/joint_vcf/joint_all_chr.vcf.gz"
    echo "  - Filtered VCF: ${OUTPUT_DIR}/filtered_vcf/filtered_variants.vcf.gz"
    echo "  - Zarr format: ${OUTPUT_DIR}/zarr/variants.zarr"
    echo ""
    echo "Next steps:"
    echo "  1. Load variants.zarr in a Python notebook"
    echo "  2. Filter by parental strain comparisons"
    echo "  3. Annotate identified mutations"
else
    echo "✗ Filtering workflow failed with exit code: ${EXIT_CODE}"
    echo ""
    echo "Check error logs:"
    echo "  - SLURM log: logs/filter_${SLURM_JOB_ID}.err"
    echo "  - Nextflow log: .nextflow.log"
    echo "  - Work directory: ${WORK_DIR}"
fi
echo "================================================================"
echo "Completed: $(date)"
echo "================================================================"
