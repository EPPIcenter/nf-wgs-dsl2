#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// QC Workflow 

//trimmomatic read trimming
process trim_reads {
    
    tag "trim ${pair_id}"    

    publishDir "${params.outputdir}/intermediate_files",
        saveAs: {filename ->
            if (filename.indexOf("_paired.fq.gz") > 0) "trimmed_pairs/$filename"
            else if (filename.indexOf("_unpaired.fq.gz") > 0) "unpaired/$filename"
            else filename
    }
        
    input:
    tuple val(pair_id), path(reads)
    path trim_adapter

    output:
    tuple val(pair_id), path("trimmed_${pair_id}_R{1,2}_paired.fq.gz"),
    path("trimmed_${pair_id}_R{1,2}_unpaired.fq.gz")

    script:
    """
    trimmomatic PE ${reads[0]} ${reads[1]} \
    "trimmed_${pair_id}_R1_paired.fq.gz" "trimmed_${pair_id}_R1_unpaired.fq.gz" \
    "trimmed_${pair_id}_R2_paired.fq.gz" "trimmed_${pair_id}_R2_unpaired.fq.gz" \
    ILLUMINACLIP:$trim_adapter:2:30:10 LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:5:20 -threads ${task.cpus}
    """
}

//fastqc on each trimmed read pair
process fastqc {
    
    tag "FASTQC on ${pair_id}"

    publishDir "${params.outputdir}/intermediate_files/fastqc"

    input:
    tuple val(pair_id), path(paired_reads), path(unpaired_reads)

    output:
    path("fastqc_${pair_id}")

    conda 'bioconda::fastqc'

    script:
    """
    mkdir fastqc_${pair_id}
    fastqc -o fastqc_${pair_id} -q ${paired_reads}
    """  
}  

//multiqc report
process multiqc {
    
    tag "multiqc on all trimmed_fastqs"

    publishDir "${params.outputdir}/final_qc_reports", mode:'copy'

    input:
    path(fastqc_results)

    output:
    path('multiqc_report.html')  

    conda 'bioconda::multiqc'

    script:
    """    
    multiqc .
    """
}

// bwa alignment
process bwa_align {
    
    tag "align ${pair_id}"

    publishDir "${params.outputdir}/intermediate_files/align"

    input:
    tuple val(pair_id), path(paired_reads), path(unpaired_reads)
    path genomes_dir

    output:
    tuple val(pair_id), path("${pair_id}.sam")
    
    script:
    """
    bwa mem -t ${task.cpus} \
    -M -R "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:illumina\\tSM:${pair_id}\\tPU:${pair_id}" \
    $genomes_dir/Pf3D7_human.fa ${paired_reads} > ${pair_id}.sam
    """
}

// sam format converter
process sam_convert {
    
    tag "sam format converter ${pair_id}"

    publishDir "${params.outputdir}/intermediate_files/align"

    input:
    tuple val(pair_id), path(sam_file)
    path genomes_dir

    output:
    tuple val(pair_id), path("${pair_id}.bam")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g" SamFormatConverter \
    -R $genomes_dir/Pf3D7_human.fa \
    -I ${sam_file} \
    -O ${pair_id}.bam

    # rm ${sam_file}
    """
}

// sam clean
process sam_clean {
    
    tag "sam cleaning ${pair_id}"

    publishDir "${params.outputdir}/intermediate_files/align"

    input:
    tuple val(pair_id), path(bam_file)
    path genomes_dir

    output:
    tuple val(pair_id), path("${pair_id}.clean.bam")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g" CleanSam \
    -R $genomes_dir/Pf3D7_human.fa \
    -I ${bam_file} \
    -O ${pair_id}.clean.bam

    # rm ${bam_file}
    """
}

// sam file sorting
process sam_sort {
    
    tag "sam sorting ${pair_id}"
    scratch true

    publishDir "${params.outputdir}/intermediate_files/align"

    input:
    tuple val(pair_id), path(clean_bam)
    path genomes_dir

    output:
    tuple val(pair_id), path("${pair_id}.sorted.bam")

    afterScript "rm -rf TMP"

    script: 
    """
    mkdir -p TMP

    # sam file sorting
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g" SortSam \
    -R $genomes_dir/Pf3D7_human.fa \
    -I ${clean_bam} \
    -O ${pair_id}.sorted.bam \
    -SO coordinate \
    --CREATE_INDEX true \
    --TMP_DIR TMP

    # rm ${clean_bam}
    """
}

// mark duplicates
process sam_duplicates {
    
    tag "sam mark duplicates ${pair_id}"
    scratch true

    publishDir "${params.outputdir}/final_bams", mode:'copy'

    input:
    tuple val(pair_id), path(sorted_bam)
    path genomes_dir

    output:
    tuple val(pair_id), path("${pair_id}.sorted.dup.bam"),
    path("${pair_id}_dup_metrics.txt")

    afterScript "rm -rf TMP"

    script:
    """
    mkdir -p TMP

    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g" MarkDuplicates \
    -R $genomes_dir/Pf3D7_human.fa \
    -I ${sorted_bam} \
    -O ${pair_id}.sorted.dup.bam \
    -M ${pair_id}_dup_metrics.txt \
    -ASO coordinate \
    --TMP_DIR TMP

    # rm ${sorted_bam}
    """
}
// subsampling Plasmodium falciparum specific bams
process target_pf {
    
    tag "target Pf ${pair_id}"

    publishDir "${params.outputdir}/final_bams", mode:'copy'

    input:
    tuple val(pair_id), path(sorted_dup_bam), path(dup_metrics_txt)
    path genomes_dir

    output:
    tuple val(pair_id), path("${pair_id}.sorted.dup.pf.bam"), path("${pair_id}.sorted.dup.pf.bam.csi")

    script:
    """
    samtools view -b -h ${sorted_dup_bam} \
    -T $genomes_dir/Pf3D7.fasta \
    -L $genomes_dir/Pf3D7_core.bed > ${pair_id}.sorted.dup.pf.bam

    samtools index -bc ${pair_id}.sorted.dup.pf.bam
    """    
}

// index pf bam (create .bai file which is needed for pathweaver)
process index_pf_bam {
    
    tag "index Pf bam ${pair_id}"

    publishDir "${params.outputdir}/final_bams", mode:'copy'

    input:
    tuple val(pair_id), path(pf_bam), path(pf_bam_index)

    output:
    path("${pair_id}.sorted.dup.pf.bam.bai")

    script:
    """
    samtools index ${pf_bam}
    """    
}

// subsampling human specific bams
process target_human {
    
    tag "target Hs ${pair_id}"

    publishDir "${params.outputdir}/intermediate_files/align"

    input:
    tuple val(pair_id), path(sorted_dup_bam), path(dup_metrics_txt)
    path genomes_dir

    output:
    tuple val(pair_id), path("${pair_id}.sorted.dup.hs.bam")

    script:
    """
    samtools view -b -h ${sorted_dup_bam} \
    -T $genomes_dir/genome.fa \
    -L $genomes_dir/human.bed > ${pair_id}.sorted.dup.hs.bam
    """    
}

// insert size calculation
process insert_sizes {
    
    tag "insert sizes ${pair_id}"

    publishDir "${params.outputdir}", mode: 'copy',
        saveAs: {filename -> filename.indexOf("_histo.pdf") > 0 ? "histograms/$filename" : "insert_size_metrics/$filename"}
        
    input:
    tuple val(pair_id), path(pf_bam), path(pf_bam_index)

    output:
    tuple val(pair_id), path("${pair_id}.insert.txt"), path("${pair_id}.insert2.txt"), path("${pair_id}_histo.pdf")
    
    script:
    """
    gatk CollectInsertSizeMetrics \
    -I ${pf_bam} \
    -O ${pair_id}.insert.txt \
    -H ${pair_id}_histo.pdf \
    -M 0.05

    # Note: MODE_INSERT_SIZE (2) and PAIR_ORIENTATION (9) are excluded 
    awk 'BEGIN{OFS="\t"} FNR>=8 && FNR<=8 {print \$1,\$3,\$4,\$5,\$6,\$7,\$8,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$NF="${pair_id}"}' ${pair_id}.insert.txt > ${pair_id}.insert2.txt
    """
}

// insert summary
process insert_summary {

    tag "insert summary"

    publishDir "${params.outputdir}/final_qc_reports", mode:'copy'

    input:
    path(insert2_files)

    output:
    path('InsertSize_Final.tsv')

    script:
    """
    cat $insert2_files | sed '1iMEDIAN_INSERT_SIZE	MEDIAN_ABSOLUTE_DEVIATION	MIN_INSERT_SIZE	MAX_INSERT_SIZE	MEAN_INSERT_SIZE	STANDARD_DEVIATION	READ_PAIRS	WIDTH_OF_10_PERCENT	WIDTH_OF_20_PERCENT	WIDTH_OF_30_PERCENT	WIDTH_OF_40_PERCENT	WIDTH_OF_50_PERCENT	WIDTH_OF_60_PERCENT	WIDTH_OF_70_PERCENT	WIDTH_OF_80_PERCENT	WIDTH_OF_90_PERCENT	WIDTH_OF_95_PERCENT	WIDTH_OF_99_PERCENT	SAMPLE_ID' > InsertSize_Final.tsv
    """
}

// Total bam statistics by sample
process total_bam_stat_per_sample {
    
    tag "total bam stat ${pair_id}"

    publishDir "${params.outputdir}/intermediate_files/stat"
    
    input: 
    tuple val(pair_id), path(sorted_dup_bam), path(dup_metrics_txt)

    output:
    tuple val(pair_id), path("${pair_id}_bamstat_total_final.tsv"), path("${pair_id}_bamstat_total.tsv")

    conda 'bioconda::samtools bioconda::datamash'

    script:
    """
    samtools stats ${sorted_dup_bam} | grep ^SN | cut -f 2- | awk -F"\t" '{print \$2}' > ${pair_id}_bamstat_total.tsv

    datamash transpose < ${pair_id}_bamstat_total.tsv | awk -F"\t" -v OFS="\t" '{ \$(NF+1) = "${pair_id}"; print }' > ${pair_id}_bamstat_total_final.tsv
    """
}

// Total bam statistic summary
process total_stat_summary {
    
    tag "total stat summary"

    publishDir "${params.outputdir}/final_qc_reports", mode:'copy'

    input:
    path(bamstat_total) 

    output:
    path('Bam_stats_total_Final.tsv') 

    script:
    """
    cat $bamstat_total  | sed '1irow_total_reads_total	filtered_reads_total	sequences_total	is_sorted_total	1st_fragments_total	last_fragments_total	reads_mapped_total	reads_mapped_and_paired_total	reads_unmapped_total	reads_properly_paired_total	reads_paired_total	reads_duplicated_total	reads_MQ0_total	reads_QC_failed_total	non_primary_alignments_total	supplementary_alignments_total	total_length_total	total_first_fragment_length_total	total_last_fragment_length_total	bases_mapped_total	bases_mapped_(cigar)_total	bases_trimmed_total	bases_duplicated_total	mismatches_total	error_rate_total	average_length_total	average_first_fragment_length_total	average_last_fragment_length_total	maximum_length_total	maximum_first_fragment_length_total	maximum_last_fragment_length_total	average_quality_total	insert_size_average_total	insert_size_standard_deviation_total	inward_oriented_pairs_total	outward_oriented_pairs_total	pairs_with_other_orientation_total	pairs_on_different_chromosomes_total	percentage_of_properly_paired_reads_(%)_total	sample_id' > Bam_stats_total_Final.tsv
    """
}

// Pf bam statistics by sample
process pf_bam_stat_per_sample {
    
    tag "Pf bam stat ${pair_id}"

    publishDir "${params.outputdir}/intermediate_files/stat"
    
    input: 
    tuple val(pair_id), path(pf_bam), path(pf_bam_index)

    output:
    tuple val(pair_id), path("${pair_id}_bamstat_pf_final.tsv"), path("${pair_id}_bamstat_pf.tsv")

    conda 'bioconda::samtools bioconda::datamash'

    script:
    """
    samtools stats ${pf_bam} | grep ^SN | cut -f 2- | awk -F"\t" '{print \$2}' > ${pair_id}_bamstat_pf.tsv

    datamash transpose < ${pair_id}_bamstat_pf.tsv | awk -F"\t" -v OFS="\t" '{ \$(NF+1) = "${pair_id}"; print }' > ${pair_id}_bamstat_pf_final.tsv
    """
}

// Pf bam statistic summary
process pf_stat_summary {
    
    tag "Pf stat summary"

    publishDir "${params.outputdir}/final_qc_reports", mode:'copy'

    input:
    path(bamstat_pf) 

    output:
    path('Bam_stats_pf_Final.tsv') 

    script:
    """
    cat $bamstat_pf  | sed '1irow_total_reads_pf	filtered_reads_pf	sequences_pf	is_sorted_pf	1st_fragments_pf	last_fragments_pf	reads_mapped_pf	reads_mapped_and_paired_pf	reads_unmapped_pf	reads_properly_paired_pf	reads_paired_pf	reads_duplicated_pf	reads_MQ0_pf	reads_QC_failed_pf	non_primary_alignments_pf	supplementary_alignments_pf	total_length_pf	total_first_fragment_length_pf	total_last_fragment_length_pf	bases_mapped_pf	bases_mapped_(cigar)_pf	bases_trimmed_pf	bases_duplicated_pf	mismatches_pf	error_rate_pf	average_length_pf	average_first_fragment_length_pf	average_last_fragment_length_pf	maximum_length_pf	maximum_first_fragment_length_pf	maximum_last_fragment_length_pf	average_quality_pf	insert_size_average_pf	insert_size_standard_deviation_pf	inward_oriented_pairs_pf	outward_oriented_pairs_pf	pairs_with_other_orientation_pf	pairs_on_different_chromosomes_pf	percentage_of_properly_paired_reads_(%)_pf	sample_id' > Bam_stats_pf_Final.tsv
    """
}

// Hs bam statistics by sample
process hs_bam_stat_per_sample {

    tag "Hs bam stat ${pair_id}"

    publishDir "${params.outputdir}/intermediate_files/stat"
    
    input:
    tuple val(pair_id), path(hs_bam)

    output:
    tuple val(pair_id), path("${pair_id}_bamstat_hs_final.tsv"), path("${pair_id}_bamstat_hs.tsv")

    conda 'bioconda::samtools bioconda::datamash'

    script:
    """
    samtools stats ${hs_bam} | grep ^SN | cut -f 2- | awk -F"\t" '{print \$2}' > ${pair_id}_bamstat_hs.tsv

    datamash transpose < ${pair_id}_bamstat_hs.tsv | awk -F"\t" -v OFS="\t" '{ \$(NF+1) = "${pair_id}"; print }' > ${pair_id}_bamstat_hs_final.tsv
    """
}

// Hs bam statistic summary
process hs_stat_summary {
    tag "Hs stat summary"

    publishDir "${params.outputdir}/final_qc_reports", mode:'copy'

    input:
    path(bamstat_hs)

    output:
    path('Bam_stats_hs_Final.tsv')

    script:
    """
    cat $bamstat_hs | sed '1irow_total_reads_hs	filtered_reads_hs	sequences_hs	is_sorted_hs	1st_fragments_hs	last_fragments_hs	reads_mapped_hs	reads_mapped_and_paired_hs	reads_unmapped_hs	reads_properly_paired_hs	reads_paired_hs	reads_duplicated_hs	reads_MQ0_hs	reads_QC_failed_hs	non_primary_alignments_hs	supplementary_alignments_hs	total_length_hs	total_first_fragment_length_hs	total_last_fragment_length_hs	bases_mapped_hs	bases_mapped_(cigar)_hs	bases_trimmed_hs	bases_duplicated_hs	mismatches_hs	error_rate_hs	average_length_hs	average_first_fragment_length_hs	average_last_fragment_length_hs	maximum_length_hs	maximum_first_fragment_length_hs	maximum_last_fragment_length_hs	average_quality_hs	insert_size_average_hs	insert_size_standard_deviation_hs	inward_oriented_pairs_hs	outward_oriented_pairs_hs	pairs_with_other_orientation_hs	pairs_on_different_chromosomes_hs	percentage_of_properly_paired_reads_(%)_hs	sample_id' > Bam_stats_hs_Final.tsv
    """
}

// distribution of Pf read depth by chromosome
process pf_read_depth {
    
    tag "read depth Pf chroms ${pair_id}"
    scratch true

    publishDir "${params.outputdir}/intermediate_files/stat"

    input:
    tuple val(pair_id), path(pf_bam), path(pf_bam_index)
    path genomes_dir

    output:
    path("${pair_id}_read_coverage.tsv")

    afterScript "rm -rf TMP"

    script: 
    """
    mkdir -p TMP

    for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
        do
            gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g" DepthOfCoverage \
            -R "$genomes_dir/Pf3D7.fasta" \
            -O chr"\$i" \
            --output-format TABLE \
            -L Pf3D7_"\$i"_v3 \
            --omit-locus-table true \
            -I ${pf_bam} --tmp-dir TMP
            sed '\${/^Total/d;}' chr"\$i".sample_summary > tmp.chr"\$i".sample_summary
            awk -F"\t" -v OFS="\t" '{ print \$0, \$(NF) = "chr'\$i'" }' tmp.chr"\$i".sample_summary > chr"\$i".sample2_summary
        done

    cat *.sample2_summary | awk '!/sample_id/ {print \$0}' | sed '1isample_id	total	mean	third_quartile	median	first_quartile	bases_perc_above_15	chromosome' > ${pair_id}_read_coverage.tsv
    """
}

// concatenate Pf read depth by chromosome summaries
process pf_read_depth_summary {
    
    tag "read depth Pf chroms summary"

    publishDir "${params.outputdir}/final_qc_reports", mode:'copy'

    input:
    path(read_coverage)

    output:
    path("ReadCoverage_final.tsv")

    script: 
    """
    cat $read_coverage | awk '!/sample_id/ {print \$0}' | sed '1isample_id	total	mean	third_quartile	median	first_quartile	bases_perc_above_15	chromosome' > ReadCoverage_final.tsv
    """
}

// Rmd run quality report generation
process run_report_and_calculate_ratio {
 
    tag "run quality report and calculate Pf:Hs ratio"

    publishDir "${params.outputdir}/final_qc_reports", mode:'copy'

    input:
    path(pf_summary)
    path(hs_summary)
    path(coverage_summary)
    path(insert1_files)
    path rscript

    output:
    tuple path('Ratios_hs_pf_reads.tsv'), path('run_quality_report.html')

    afterScript "rm -rf TMP"

    script:
    """    
    mkdir -p TMP/INSERT_FILES
    cp $pf_summary TMP
    cp $hs_summary TMP
    cp $coverage_summary TMP
    cp $insert1_files TMP/INSERT_FILES/
    Rscript -e 'rmarkdown::render(input = "$rscript", output_dir = getwd(), params = list(directory = "TMP"))'
    """
}

workflow.onComplete { 
    println ( workflow.success ? "\nQC run complete!": "Oops .. something went wrong" )
}

workflow QC {
    main: 
        read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )
        // trim reads
        trimmed_reads_ch = trim_reads(read_pairs_ch, params.trim_adapter)

        // fastqc report 
        fastqc_ch = fastqc(trimmed_reads_ch)
        // multiqc report --
        multiqc(fastqc_ch.collect()) 

        // bwa alignment
        sam_ch = bwa_align(trimmed_reads_ch, params.genomes_dir)

        // sam format converter, clean, sort, and mark duplicates
        bam_ch = sam_convert(sam_ch, params.genomes_dir)
        bam_clean_ch = sam_clean(bam_ch, params.genomes_dir)
        bam_sort_ch = sam_sort(bam_clean_ch, params.genomes_dir)
        bam_dup_ch = sam_duplicates(bam_sort_ch, params.genomes_dir)

        // samtools sorting Pf and human reads
        pf_bam_ch = target_pf(bam_dup_ch, params.genomes_dir)
        hs_bam_ch = target_human(bam_dup_ch, params.genomes_dir)

        // index pf bam
        index_pf_bam(pf_bam_ch) 

        // distribution of Pf read depth by chromosome -- 
        pf_read_depth_ch = pf_read_depth(pf_bam_ch, params.genomes_dir) 
        coverage_summary_ch = pf_read_depth_summary(pf_read_depth_ch.collect())

        // insert size calculation
        inserts_ch = insert_sizes(pf_bam_ch) 
        insert1_ch = inserts_ch.map{T->[T[1]]} // select *.insert.txt
        insert2_ch = inserts_ch.map{T->[T[2]]} // select *.insert2.txt
        // insert summary -- 
        insert_summary(insert2_ch.collect()) 

        // Total bam statistics by sample
        total_bamstat_ch = total_bam_stat_per_sample(bam_dup_ch)
        total_final_bamstat_ch = total_bamstat_ch.map{T->[T[1]]} // select *_bamstat_total_final.tsv
        // Total bam statistic summary
        total_summary_ch = total_stat_summary(total_final_bamstat_ch.collect()) 

        // Pf bam statistics by sample
        pf_bamstat_ch = pf_bam_stat_per_sample(pf_bam_ch)
        pf_final_bamstat_ch = pf_bamstat_ch.map{T->[T[1]]} // select *_bamstat_pf_final.tsv
        // Pf bam statistic summary
        pf_summary_ch = pf_stat_summary(pf_final_bamstat_ch.collect()) 

        // Hs bam statistics by sample
        hs_bamstat_ch = hs_bam_stat_per_sample(hs_bam_ch)
        hs_final_bamstat_ch = hs_bamstat_ch.map{T->[T[1]]} // select *_bamstat_hs_final.tsv
        // Hs bam statistic summary
        hs_summary_ch = hs_stat_summary(hs_final_bamstat_ch.collect())

        // Rmd run quality report generation and Calculate Pf:Hs read ratio -- 
        run_report_and_calculate_ratio(pf_summary_ch, hs_summary_ch, coverage_summary_ch, insert1_ch.collect(), params.rscript)
        
    emit: pf_bam_ch

}