/*
 * Variant Filtering Workflow for Lab Mutants
 * 
 * This workflow is designed for lab-generated mutants where you expect
 * few variants (induced mutations + background differences).
 * 
 * Steps:
 * 1. Joint genotype gVCFs per chromosome using GenomicsDBImport + GenotypeGVCFs
 * 2. Concatenate all chromosomes into a single VCF
 * 3. Apply hard filters (quality thresholds appropriate for mutation detection)
 * 4. Convert to zarr format for downstream analysis
 */

workflow FILTER {
    take:
        gvcf_ch  // tuple(sample_id, chr, gvcf, gvcf_idx)
    
    main:
        // Get reference genome and interval files
        genomes_dir = Channel.fromPath(params.genomes_dir, checkIfExists: true).first()
        
        // Group gVCFs by chromosome for GenomicsDB import
        gvcf_by_chr = gvcf_ch
            .map { sample_id, chr, gvcf, gvcf_idx -> 
                tuple(chr, sample_id, gvcf, gvcf_idx) 
            }
            .groupTuple(by: 0)
        
        // Step 1: Import gVCFs into GenomicsDB (per chromosome)
        genomicsdb_import(gvcf_by_chr, genomes_dir)
        
        // Step 2: Joint genotype each chromosome
        joint_genotype(genomicsdb_import.out, genomes_dir)
        
        // Step 3: Concatenate all chromosomes  
        // Sort VCFs by chromosome number then split into vcfs and indices
        sorted_outputs = joint_genotype.out
            .toSortedList { a, b -> 
                // Sort by chromosome number (a[0] and b[0] are the chr values)
                a[0].toInteger() <=> b[0].toInteger()
            }
            .map { sorted_list ->
                // Extract VCFs and indices into separate lists
                def vcfs = sorted_list.collect { it[1] }  // it[1] is the VCF file
                def indices = sorted_list.collect { it[2] }  // it[2] is the index file
                tuple(vcfs, indices)
            }
        
        concat_chromosomes(sorted_outputs, genomes_dir)
        
        // Step 4: Apply hard filters
        hard_filter(concat_chromosomes.out, genomes_dir)
        
        // Step 5: Convert to zarr (optional, controlled by param)
        if (params.output_zarr) {
            vcf_to_zarr(hard_filter.out.all_variants, hard_filter.out.pass_only)
        }
    
    emit:
        filtered_all = hard_filter.out.all_variants
        filtered_pass = hard_filter.out.pass_only
        zarr_all = params.output_zarr ? vcf_to_zarr.out.zarr_all : Channel.empty()
        zarr_pass = params.output_zarr ? vcf_to_zarr.out.zarr_pass : Channel.empty()
}

// Import gVCFs into GenomicsDB workspace (per chromosome)
process genomicsdb_import {
    
    tag "GenomicsDB import chr${chrom}"
    label 'big_mem'
    
    input:
    tuple val(chrom), val(sample_ids), path(gvcfs), path(gvcf_indices)
    path genomes_dir
    
    output:
    tuple val(chrom), path("genomicsdb_chr${chrom}")
    
    script:
    def sample_map = sample_ids.withIndex().collect { id, idx -> 
        "${id}\t${gvcfs[idx]}" 
    }.join('\n')
    """
    # Create sample map file
    cat > sample_map.txt << EOF
${sample_map}
EOF
    
    # Import to GenomicsDB
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenomicsDBImport \\
        --sample-name-map sample_map.txt \\
        --genomicsdb-workspace-path genomicsdb_chr${chrom} \\
        -L ${genomes_dir}/core_chr${chrom}.list \\
        --reader-threads ${task.cpus}
    """
}

// Joint genotype from GenomicsDB
process joint_genotype {
    
    tag "Joint genotype chr${chrom}"
    label 'big_mem'
    
    input:
    tuple val(chrom), path(genomicsdb)
    path genomes_dir
    
    output:
    tuple val(chrom), path("joint_chr${chrom}.vcf.gz"), path("joint_chr${chrom}.vcf.gz.tbi")
    
    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \\
        -R ${genomes_dir}/Pf3D7.fasta \\
        -V gendb://${genomicsdb} \\
        -L ${genomes_dir}/core_chr${chrom}.list \\
        -ploidy ${params.ploidy} \
        -G StandardAnnotation \\
        -G AS_StandardAnnotation \\
        -A ExcessHet \\
        -A InbreedingCoeff \\
        -O joint_chr${chrom}.vcf.gz
    """
}

// Concatenate all chromosomes
process concat_chromosomes {
    
    tag "Concatenate all chromosomes"
    label 'med_mem'
    
    publishDir "${params.outputdir}/joint_vcf", mode:'copy'
    
    input:
    tuple path(vcfs), path(indices)
    path genomes_dir
    
    output:
    tuple path("joint_all_chr.vcf.gz"), path("joint_all_chr.vcf.gz.tbi")
    
    script:
    // VCFs are already sorted by chromosome number
    def vcf_inputs = vcfs.collect { "-I ${it}" }.join(' ')
    """
    # Concatenate VCFs using GATK GatherVcfs
    gatk GatherVcfs \\
        ${vcf_inputs} \\
        -O joint_all_chr.vcf.gz
    
    # Index the concatenated VCF
    gatk IndexFeatureFile -I joint_all_chr.vcf.gz
    """
}

// Apply hard filters (appropriate for lab mutants)
process hard_filter {
    
    tag "Hard filter variants"
    label 'med_mem'
    
    publishDir "${params.outputdir}/filtered_vcf", mode:'copy'
    
    input:
    tuple path(vcf), path(vcf_index)
    path genomes_dir
    
    output:
    tuple path("filtered_all_variants.vcf.gz"), path("filtered_all_variants.vcf.gz.tbi"), emit: all_variants
    tuple path("filtered_pass_only.vcf.gz"), path("filtered_pass_only.vcf.gz.tbi"), emit: pass_only
    path "filtering_stats.txt"
    
    script:
    """
    # Apply hard filters using GATK VariantFiltration
    # Site-level filters for technical artifacts
    # Genotype-level filters for quality and clonal samples (heterozygous = artifact)
    gatk --java-options "-Xmx${task.memory.toGiga()}g" VariantFiltration \\
        -R ${genomes_dir}/Pf3D7.fasta \\
        -V ${vcf} \\
        -O filtered_all_variants.vcf.gz \\
        --filter-name "QD_filter"           --filter-expression "QD < ${params.filter_QD}" \\
        --filter-name "FS_filter"           --filter-expression "FS > ${params.filter_FS}" \\
        --filter-name "MQ_filter"           --filter-expression "MQ < ${params.filter_MQ}" \\
        --filter-name "MQRankSum_filter"    --filter-expression "MQRankSum < ${params.filter_MQRankSum}" \\
        --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < ${params.filter_ReadPosRankSum}" \\
        --filter-name "SOR_filter"          --filter-expression "SOR > ${params.filter_SOR}" \\
        --genotype-filter-name "GQ_filter"  --genotype-filter-expression "GQ < ${params.filter_GQ}" \\
        --genotype-filter-name "DP_filter"  --genotype-filter-expression "DP < ${params.filter_DP}" \\
        --genotype-filter-name "AD_filter"  --genotype-filter-expression "AD[1] < ${params.filter_AD}"
    
    
    # Extract only PASS variants (site and genotype filters)
    # Keep all annotations including AD (Allele Depth)
    gatk SelectVariants \\
        -R ${genomes_dir}/Pf3D7.fasta \\
        -V filtered_all_variants.vcf.gz \\
        --exclude-filtered \\
        --set-filtered-gt-to-nocall \\
        --keep-original-ac \\
        -O filtered_pass_only.vcf.gz
    
    # Index the final VCF using GATK
    gatk IndexFeatureFile -I filtered_pass_only.vcf.gz
    
    # Generate filtering statistics using GATK tools
    echo "=== Filtering Statistics ===" > filtering_stats.txt
    echo "Total variants before filtering:" >> filtering_stats.txt
    gatk CountVariants -V filtered_all_variants.vcf.gz >> filtering_stats.txt
    echo "" >> filtering_stats.txt
    echo "PASS variants after all filters:" >> filtering_stats.txt
    gatk CountVariants -V filtered_pass_only.vcf.gz >> filtering_stats.txt
    """
}

// Convert VCF to zarr format for analysis
process vcf_to_zarr {
    
    tag "Convert to zarr"
    label 'med_mem'
    conda "${projectDir}/conf/envs/zarr_env.yml"
    
    publishDir "${params.outputdir}/zarr", mode:'copy'
    
    input:
    tuple path(vcf_all), path(vcf_all_index)
    tuple path(vcf_pass), path(vcf_pass_index)
    
    output:
    path "variants_all.zarr", emit: zarr_all
    path "variants_pass.zarr", emit: zarr_pass
    
    script:
    """
    #!/usr/bin/env python3
    import allel
    import zarr
    import numcodecs
    
    # Convert all variants (with FILTER annotations) to zarr
    print("\\n=== Converting all variants (with FILTER column) ===")
    allel.vcf_to_zarr(
        '${vcf_all}',
        'variants_all.zarr',
        fields='*',
        overwrite=True,
        compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=2)
    )
    
    # Convert PASS-only variants to zarr
    print("\\n=== Converting PASS-only variants ===")
    allel.vcf_to_zarr(
        '${vcf_pass}',
        'variants_pass.zarr',
        fields='*',
        overwrite=True,
        compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=2)
    )
    
    # Print summaries
    print("\\n=== All Variants Zarr Summary ===")
    callset_all = zarr.open_group('variants_all.zarr', mode='r')
    variants_all = callset_all['variants']
    print(f"Total variants: {len(variants_all['POS'])}")
    if 'FILTER_PASS' in variants_all:
        n_pass = variants_all['FILTER_PASS'][:].sum()
        print(f"PASS variants: {n_pass}")
    
    print("\\n=== PASS-Only Zarr Summary ===")
    callset_pass = zarr.open_group('variants_pass.zarr', mode='r')
    variants_pass = callset_pass['variants']
    print(f"Total variants: {len(variants_pass['POS'])}")
    print(f"Chromosomes: {sorted(set(variants_pass['CHROM'][:]))}")
    """
}
