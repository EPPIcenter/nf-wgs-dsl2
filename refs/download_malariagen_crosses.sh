#!/bin/bash
# Create directory for cross VCFs
CROSS_DIR="malariagen_crosses"
mkdir -p $CROSS_DIR
cd $CROSS_DIR

echo "Downloading MalariaGEN genetic crosses VCF files..."

# go to https://www.malariagen.net/data_package/pf-crosses-1-0/ and click FTP, which will link you to SANGER, chose guest
# on personal laptop
cp /Volumes/1.0/3d7_hb3.combined.final.vcf.gz* /Volumes/babs/working/whittog/pipelines/nf-wgs-pf-dsl2/refs/malariagen_crosses/
cp /Volumes/1.0/7g8_gb4.combined.final.vcf.gz* /Volumes/babs/working/whittog/pipelines/nf-wgs-pf-dsl2/refs/malariagen_crosses/
cp /Volumes/1.0/hb3_dd2.combined.final.vcf.gz* /Volumes/babs/working/whittog/pipelines/nf-wgs-pf-dsl2/refs/malariagen_crosses/
