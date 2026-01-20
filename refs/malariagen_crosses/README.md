# MalariaGEN Genetic Crosses for VQSR

This directory should contain VCF files from MalariaGEN genetic crosses used as training data for Variant Quality Score Recalibration (VQSR).

## Required Files

The VQSR workflow expects the following files in this directory:

1. **7G8 × GB4 cross**:
   - `7G8_GB4.vcf.gz`
   - `7G8_GB4.vcf.gz.tbi`

2. **HB3 × Dd2 cross**:
   - `HB3_Dd2.vcf.gz`
   - `HB3_Dd2.vcf.gz.tbi`

3. **3D7 × HB3 cross**:
   - `3D7_HB3.vcf.gz`
   - `3D7_HB3.vcf.gz.tbi`

## Download Instructions

Run the download script from the `refs` directory:

```bash
cd refs
./download_malariagen_crosses.sh
```

## Alternative Sources

If the automated download script doesn't work, you can manually download the files from:

### MalariaGEN Pf7 Data Portal
- Main page: https://www.malariagen.net/data/pf7-release-7
- FTP site: ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf7/

### Manual Download from FTP
```bash
wget ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf7/crosses/7G8_GB4.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf7/crosses/7G8_GB4.vcf.gz.tbi
wget ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf7/crosses/HB3_Dd2.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf7/crosses/HB3_Dd2.vcf.gz.tbi
wget ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf7/crosses/3D7_HB3.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf7/crosses/3D7_HB3.vcf.gz.tbi
```

## File Verification

After downloading, verify that you have all required files:

```bash
ls -lh *.vcf.gz*
```

You should see 6 files total (3 VCF files + 3 index files).

## Using Custom Cross VCFs

If you have your own high-confidence cross VCF files, you can use them instead. Ensure they:
1. Are bgzip-compressed (`.vcf.gz`)
2. Have tabix indices (`.vcf.gz.tbi`)
3. Use the same reference genome (Pf3D7)
4. Contain biallelic variants with high-quality genotype calls

Place your custom VCF files in this directory and they will be automatically used by the VQSR workflow.

## About the Crosses

These genetic crosses are progeny from controlled laboratory crosses between *P. falciparum* strains:

- **7G8 × GB4**: South American × Asian strains
- **HB3 × Dd2**: Central American × Indochina strains  
- **3D7 × HB3**: African × Central American strains

The crosses provide high-confidence variant calls because:
- Controlled genetic backgrounds
- Deep sequencing coverage
- Segregation patterns validate true variants
- Multiple progeny for validation

## References

1. Miles, A., Iqbal, Z., Vauterin, P. et al. (2016). Indels, structural variation, and recombination drive genomic diversity in Plasmodium falciparum. Genome Res. 26, 1288–1299.

2. Pf3k Project (2016). Pilot data release 5. https://www.malariagen.net/data/pf3k-5

3. MalariaGEN et al. (2021). An open dataset of Plasmodium falciparum genome variation in 7,000 worldwide samples. Wellcome Open Res. 6, 42.
