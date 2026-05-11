#!/bin/bash
#SBATCH --job-name=pf-wgs_pipeline
#SBATCH --output=logs/pf-wgs_%j.out
#SBATCH --error=logs/pf-wgs_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=ncpu

# Full pipeline execution: QC -> gVCF -> Filtering
# Results saved to: output/dir

# Set Nextflow options
export NXF_OPTS='-Xms1g -Xmx4g'

# Load required modules 
module purge
module load Singularity/3.11.3
module load Python/3.11.5-GCCcore-13.2.0
module load Nextflow/25.04.4
module load git
ssh-add ~/.ssh/id_ed25519

# Create logs directory if it doesn't exist
mkdir -p logs

# Run full pipeline (QC -> gVCF -> Filtering)
nextflow run /path/to/repo/nf-wgs-pf-dsl2/main.nf \
  -profile slurm,singularity \
  --inputdir inputs \
  --outputdir path/to/outputs/results_prod \
  -resume
