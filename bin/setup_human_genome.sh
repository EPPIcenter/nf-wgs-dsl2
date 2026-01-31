#!/bin/bash
# Script to set up combined Pf3D7 + Human genome for pipeline with human contamination

# Get the project root directory (assuming this script is in bin/)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
GENOMES_DIR="${PROJECT_DIR}/refs/genomes"
HUMAN_REF="/nemo/stp/babs/reference/Genomics/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta"

echo "Setting up combined Pf3D7 + Human genome..."
echo "Project directory: $PROJECT_DIR"
echo "Genomes directory: $GENOMES_DIR"

# Create symlink to human genome
if [ ! -e "$GENOMES_DIR/genome.fa" ]; then
    echo "Creating symlink to human genome..."
    ln -s "$HUMAN_REF/genome.fa" "$GENOMES_DIR/genome.fa"
    ln -s "$HUMAN_REF/genome.fa.fai" "$GENOMES_DIR/genome.fa.fai"
    ln -s "$HUMAN_REF/genome.dict" "$GENOMES_DIR/genome.dict"
    echo "✓ Human genome symlinks created"
else
    echo "✓ Human genome files already exist"
fi

# Create combined Pf3D7 + Human genome
if [ ! -f "$GENOMES_DIR/Pf3D7_human.fa" ]; then
    echo "Creating combined Pf3D7 + Human genome..."
    cat "$GENOMES_DIR/Pf3D7.fasta" "$HUMAN_REF/genome.fa" > "$GENOMES_DIR/Pf3D7_human.fa"
    echo "✓ Combined genome created: Pf3D7_human.fa"
    
    echo "Indexing combined genome with samtools..."
    ml SAMtools/1.18-GCC-12.3.0
    samtools faidx "$GENOMES_DIR/Pf3D7_human.fa"
    echo "✓ samtools index created"
    
    echo "Creating dictionary with GATK..."
    ml GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8
    gatk CreateSequenceDictionary -R "$GENOMES_DIR/Pf3D7_human.fa" -O "$GENOMES_DIR/Pf3D7_human.dict"
    echo "✓ GATK dictionary created"
    
    echo "Indexing combined genome with BWA..."
    ml BWA/0.7.17-GCCcore-12.2.0
    bwa index "$GENOMES_DIR/Pf3D7_human.fa"
    echo "✓ BWA index created"
else
    echo "✓ Combined genome already exists: Pf3D7_human.fa"
fi

# Create human.bed file (BED file containing human chromosome regions)
if [ ! -f "$GENOMES_DIR/human.bed" ]; then
    echo "Creating human.bed file from human genome..."
    # Extract human chromosome names and lengths from the combined fasta index
    # Human chromosomes start after the Pf chromosomes in the combined file
    awk '$1 ~ /^chr/ {print $1"\t0\t"$2}' "$GENOMES_DIR/Pf3D7_human.fa.fai" > "$GENOMES_DIR/human.bed"
    echo "✓ human.bed created with $(wc -l < $GENOMES_DIR/human.bed) human chromosome regions"
else
    echo "✓ human.bed already exists"
fi

echo ""
echo "Setup complete! Files created:"
echo "  - Pf3D7_human.fa (combined Pf + human genome)"
echo "  - Pf3D7_human.fa.fai (samtools index)"
echo "  - Pf3D7_human.dict (GATK dictionary)"
echo "  - Pf3D7_human.fa.{amb,ann,bwt,pac,sa} (BWA indices)"
echo "  - human.bed (human chromosome regions)"
echo "  - genome.fa -> symlink to human reference"
echo ""
echo "You can now run the pipeline with human contamination detection enabled."
