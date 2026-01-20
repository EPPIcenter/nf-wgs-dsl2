#!/bin/bash

# Example run script for VQSR workflow
# Modify the paths and parameters according to your setup

# Set paths
INPUT_DIR="/path/to/gvcf_files"  # Directory containing *.chr*.g.vcf files
OUTPUT_DIR="/path/to/results"
WORK_DIR="/path/to/work"  # Nextflow work directory

# VQSR parameters
SNP_FILTER=99.0    # Truth sensitivity for SNPs (90.0-99.9)
INDEL_FILTER=99.0  # Truth sensitivity for INDELs (90.0-99.9)

# Run VQSR workflow
nextflow run main.nf \
  --vqsr_only \
  --inputdir ${INPUT_DIR} \
  --outputdir ${OUTPUT_DIR} \
  --vqsr_snp_filter_level ${SNP_FILTER} \
  --vqsr_indel_filter_level ${INDEL_FILTER} \
  --concat_chromosomes true \
  -profile apptainer \
  -work-dir ${WORK_DIR} \
  -resume \
  -with-report vqsr_report.html \
  -with-timeline vqsr_timeline.html \
  -with-dag vqsr_dag.html

echo "VQSR workflow complete!"
echo "Check output in: ${OUTPUT_DIR}"
