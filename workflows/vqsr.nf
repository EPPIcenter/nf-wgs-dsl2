#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Variant Quality Score Recalibration (VQSR) Workflow
// This workflow performs joint genotyping and variant recalibration using MalariaGEN crosses

// Import gVCFs into GenomicsDB
process genomicsdb_import {
    
    tag "GenomicsDB import chr${chrom}"
    label 'big_mem'
    
    publishDir "${params.outputdir}/genomicsdb", mode:'copy'
       
    input:
    tuple val(chrom), path(gvcfs), path(gvcf_indices)
    path genomes_dir

    output:
    tuple val(chrom), path("genomicsdb_chr${chrom}")

    script:
    def gvcf_args = gvcfs.collect { "-V $it" }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenomicsDBImport \
    ${gvcf_args} \
    --genomicsdb-workspace-path genomicsdb_chr${chrom} \
    -L $genomes_dir/core_chr${chrom}.list \
    --reader-threads ${task.cpus} \
    --batch-size 50
    """
}

// Joint genotyping with GenotypeGVCFs
process joint_genotype {
    
    tag "Joint genotype chr${chrom}"
    label 'big_mem'
    
    publishDir "${params.outputdir}/joint_vcf", mode:'copy'
       
    input:
    tuple val(chrom), path(genomicsdb)
    path genomes_dir

    output:
    tuple val(chrom), path("joint_chr${chrom}.vcf.gz"), path("joint_chr${chrom}.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \
    -R $genomes_dir/Pf3D7.fasta \
    -V gendb://${genomicsdb} \
    -O joint_chr${chrom}.vcf.gz \
    --heterozygosity 0.0029 \
    --indel-heterozygosity 0.0017 \
    -L $genomes_dir/core_chr${chrom}.list
    """
}

// Select SNPs for recalibration
process select_snps {
    
    tag "Select SNPs chr${chrom}"
    label 'med_mem'
    
    publishDir "${params.outputdir}/filtered_vcf", mode:'copy'
       
    input:
    tuple val(chrom), path(vcf), path(vcf_index)
    path genomes_dir

    output:
    tuple val(chrom), path("snps_chr${chrom}.vcf.gz"), path("snps_chr${chrom}.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" SelectVariants \
    -R $genomes_dir/Pf3D7.fasta \
    -V ${vcf} \
    -select-type SNP \
    -O snps_chr${chrom}.vcf.gz
    """
}

// Select INDELs for recalibration
process select_indels {
    
    tag "Select INDELs chr${chrom}"
    label 'med_mem'
    
    publishDir "${params.outputdir}/filtered_vcf", mode:'copy'
       
    input:
    tuple val(chrom), path(vcf), path(vcf_index)
    path genomes_dir

    output:
    tuple val(chrom), path("indels_chr${chrom}.vcf.gz"), path("indels_chr${chrom}.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" SelectVariants \
    -R $genomes_dir/Pf3D7.fasta \
    -V ${vcf} \
    -select-type INDEL \
    -O indels_chr${chrom}.vcf.gz
    """
}

// Variant Recalibration for SNPs
process variant_recalibrator_snps {
    
    tag "VQSR SNPs chr${chrom}"
    label 'big_mem'
    
    publishDir "${params.outputdir}/vqsr", mode:'copy'
       
    input:
    tuple val(chrom), path(snp_vcf), path(snp_vcf_index)
    path genomes_dir
    path cross_vcfs

    output:
    tuple val(chrom), path("snps_chr${chrom}.recal"), path("snps_chr${chrom}.recal.idx"), path("snps_chr${chrom}.tranches"), path("snps_chr${chrom}.plots.R")

    script:
    def resource_args = cross_vcfs.collect { 
        "--resource:cross,known=false,training=true,truth=true,prior=15.0 $it" 
    }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" VariantRecalibrator \
    -R $genomes_dir/Pf3D7.fasta \
    -V ${snp_vcf} \
    ${resource_args} \
    -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
    -mode SNP \
    -O snps_chr${chrom}.recal \
    --tranches-file snps_chr${chrom}.tranches \
    --rscript-file snps_chr${chrom}.plots.R \
    --max-gaussians 8 \
    --trust-all-polymorphic
    """
}

// Variant Recalibration for INDELs
process variant_recalibrator_indels {
    
    tag "VQSR INDELs chr${chrom}"
    label 'big_mem'
    
    publishDir "${params.outputdir}/vqsr", mode:'copy'
       
    input:
    tuple val(chrom), path(indel_vcf), path(indel_vcf_index)
    path genomes_dir
    path cross_vcfs

    output:
    tuple val(chrom), path("indels_chr${chrom}.recal"), path("indels_chr${chrom}.recal.idx"), path("indels_chr${chrom}.tranches"), path("indels_chr${chrom}.plots.R")

    script:
    def resource_args = cross_vcfs.collect { 
        "--resource:cross,known=false,training=true,truth=true,prior=15.0 $it" 
    }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" VariantRecalibrator \
    -R $genomes_dir/Pf3D7.fasta \
    -V ${indel_vcf} \
    ${resource_args} \
    -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
    -mode INDEL \
    -O indels_chr${chrom}.recal \
    --tranches-file indels_chr${chrom}.tranches \
    --rscript-file indels_chr${chrom}.plots.R \
    --max-gaussians 8 \
    --trust-all-polymorphic
    """
}

// Apply VQSR to SNPs
process apply_vqsr_snps {
    
    tag "Apply VQSR SNPs chr${chrom}"
    label 'med_mem'
    
    publishDir "${params.outputdir}/recalibrated_vcf", mode:'copy'
       
    input:
    tuple val(chrom), path(snp_vcf), path(snp_vcf_index), path(recal), path(recal_idx), path(tranches), path(plots)
    path genomes_dir

    output:
    tuple val(chrom), path("snps_recal_chr${chrom}.vcf.gz"), path("snps_recal_chr${chrom}.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" ApplyVQSR \
    -R $genomes_dir/Pf3D7.fasta \
    -V ${snp_vcf} \
    --recal-file ${recal} \
    --tranches-file ${tranches} \
    -mode SNP \
    --lod-score-cutoff ${params.vqsr_snp_lod_cutoff} \
    -O snps_recal_chr${chrom}.vcf.gz
    """
}

// Apply VQSR to INDELs
process apply_vqsr_indels {
    
    tag "Apply VQSR INDELs chr${chrom}"
    label 'med_mem'
    
    publishDir "${params.outputdir}/recalibrated_vcf", mode:'copy'
       
    input:
    tuple val(chrom), path(indel_vcf), path(indel_vcf_index), path(recal), path(recal_idx), path(tranches), path(plots)
    path genomes_dir

    output:
    tuple val(chrom), path("indels_recal_chr${chrom}.vcf.gz"), path("indels_recal_chr${chrom}.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" ApplyVQSR \
    -R $genomes_dir/Pf3D7.fasta \
    -V ${indel_vcf} \
    --recal-file ${recal} \
    --tranches-file ${tranches} \
    -mode INDEL \
    --lod-score-cutoff ${params.vqsr_indel_lod_cutoff} \
    -O indels_recal_chr${chrom}.vcf.gz
    """
}

// Merge recalibrated SNPs and INDELs
process merge_recalibrated_variants {
    
    tag "Merge variants chr${chrom}"
    label 'med_mem'
    
    publishDir "${params.outputdir}/final_vcf", mode:'copy'
       
    input:
    tuple val(chrom), path(snp_vcf), path(snp_vcf_index), path(indel_vcf), path(indel_vcf_index)
    path genomes_dir

    output:
    tuple val(chrom), path("recalibrated_chr${chrom}.vcf.gz"), path("recalibrated_chr${chrom}.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" MergeVcfs \
    -I ${snp_vcf} \
    -I ${indel_vcf} \
    -O recalibrated_chr${chrom}.vcf.gz
    """
}

// Concatenate all chromosomes into a single VCF
process concat_chromosomes {
    
    tag "Concatenate all chromosomes"
    label 'big_mem'
    
    publishDir "${params.outputdir}/final_vcf", mode:'copy'
       
    input:
    path(vcfs)
    path(indices)

    output:
    tuple path("recalibrated_all.vcf.gz"), path("recalibrated_all.vcf.gz.tbi")

    script:
    def vcf_list = vcfs.collect { "$it" }.join(' ')
    """
    bcftools concat -O z -o recalibrated_all.vcf.gz ${vcf_list}
    bcftools index -t recalibrated_all.vcf.gz
    """
}

// Convert VCF to zarr format using scikit-allel (Pf7 methods)
process vcf_to_zarr {
    
    tag "Convert VCF to zarr format"
    label 'big_mem'
    
    publishDir "${params.outputdir}/final_zarr", mode:'copy'
       
    input:
    tuple path(vcf), path(vcf_index)

    output:
    path("recalibrated_all.zarr")

    script:
    """
    #!/usr/bin/env python3
    import allel
    import zarr
    import numcodecs
    
    # Read VCF and convert to zarr
    # Following Pf7 methods: zarr v2.4.0 format using scikit-allel v1.2.1
    allel.vcf_to_zarr(
        '${vcf}',
        'recalibrated_all.zarr',
        group='.',
        fields='*',
        alt_number=7,
        log=None,
        compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False)
    )
    
    print("VCF successfully converted to zarr format")
    """
}

workflow.onComplete { 
    println ( workflow.success ? "\nVQSR workflow complete!": "Oops .. something went wrong" )
}

workflow VQSR {

    take: 
        gvcf_ch // channel with gVCF files from previous workflow
    
    main:
        // Create chromosome channel
        chrom_ch = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
        
        // Group gVCFs by chromosome
        // Expected input: tuple(sample_id, chr, gvcf, gvcf_idx)
        gvcf_grouped = gvcf_ch
            .map { sample_id, chr, gvcf, gvcf_idx -> 
                tuple(chr, gvcf, gvcf_idx) 
            }
            .groupTuple(by: 0)
        
        // Load MalariaGEN cross VCFs (these should be in params.cross_vcfs_dir)
        cross_vcfs_ch = Channel.fromPath("${params.cross_vcfs_dir}/*.vcf.gz")
        
        // Import into GenomicsDB
        genomicsdb_ch = genomicsdb_import(gvcf_grouped, params.genomes_dir)
        
        // Joint genotyping
        joint_vcf_ch = joint_genotype(genomicsdb_ch, params.genomes_dir)
        
        // Select variants
        snps_ch = select_snps(joint_vcf_ch, params.genomes_dir)
        indels_ch = select_indels(joint_vcf_ch, params.genomes_dir)
        
        // Collect cross VCFs for use in recalibration
        cross_vcfs_collected = cross_vcfs_ch.collect()
        
        // Variant recalibration
        snp_recal_ch = variant_recalibrator_snps(snps_ch, params.genomes_dir, cross_vcfs_collected)
        indel_recal_ch = variant_recalibrator_indels(indels_ch, params.genomes_dir, cross_vcfs_collected)
        
        // Apply VQSR
        snps_vqsr_ch = apply_vqsr_snps(
            snps_ch.join(snp_recal_ch, by: 0), 
            params.genomes_dir
        )
        indels_vqsr_ch = apply_vqsr_indels(
            indels_ch.join(indel_recal_ch, by: 0), 
            params.genomes_dir
        )
        
        // Merge SNPs and INDELs
        merged_ch = merge_recalibrated_variants(
            snps_vqsr_ch.join(indels_vqsr_ch, by: 0),
            params.genomes_dir
        )
        
        // Concatenate all chromosomes and convert to zarr
        if (params.concat_chromosomes) {
            all_vcfs = merged_ch.map { chrom, vcf, idx -> vcf }.collect()
            all_indices = merged_ch.map { chrom, vcf, idx -> idx }.collect()
            final_vcf_ch = concat_chromosomes(all_vcfs, all_indices)
            
            // Convert final VCF to zarr format (Pf7 methods)
            if (params.output_zarr) {
                zarr_ch = vcf_to_zarr(final_vcf_ch)
            }
        }

    emit:
        merged_ch
}
