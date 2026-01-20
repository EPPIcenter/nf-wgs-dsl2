#!/bin/bash

# Download MalariaGEN genetic crosses VCF files for VQSR training
# These crosses are from the Pf7 release and are used as high-confidence truth sets
# Reference: https://www.malariagen.net/data/pf7-release-7

# Create directory for cross VCFs
CROSS_DIR="malariagen_crosses"
mkdir -p $CROSS_DIR
cd $CROSS_DIR

echo "Downloading MalariaGEN genetic crosses VCF files..."

# Base URL for Pf7 data
BASE_URL="https://www.malariagen.net/sites/default/files/Pf7_vcf"

# Download genetic crosses
# 7G8 x GB4 cross
echo "Downloading 7G8xGB4 cross..."
wget -c ${BASE_URL}/crosses/7G8_GB4.vcf.gz
wget -c ${BASE_URL}/crosses/7G8_GB4.vcf.gz.tbi

# HB3 x Dd2 cross
echo "Downloading HB3xDd2 cross..."
wget -c ${BASE_URL}/crosses/HB3_Dd2.vcf.gz
wget -c ${BASE_URL}/crosses/HB3_Dd2.vcf.gz.tbi

# 3D7 x HB3 cross
echo "Downloading 3D7xHB3 cross..."
wget -c ${BASE_URL}/crosses/3D7_HB3.vcf.gz
wget -c ${BASE_URL}/crosses/3D7_HB3.vcf.gz.tbi

echo "Download complete!"
echo ""
echo "Cross VCF files have been downloaded to: $(pwd)"
echo ""
echo "These files contain high-confidence variant calls from genetic crosses"
echo "and will be used as training data for Variant Quality Score Recalibration (VQSR)."
echo ""
echo "References:"
echo "  - MalariaGEN Pf7: https://www.malariagen.net/data/pf7-release-7"
echo "  - Paper: https://www.nature.com/articles/s41586-021-03819-6"

cd ..
