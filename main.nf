/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

params.rscript = "$projectDir/refs/run_quality_report.Rmd" 
params.reads = "${params.inputdir}/*_R{1,2}*.fastq.gz" // if start from QC 
params.bams = "${params.inputdir}/*.sorted.dup.pf.{bam,bam.csi}" // if start from GVCF

log.info """\
W G S - P I P E L I N E!
================================
inputdir        : $params.inputdir
outputdir       : $params.outputdir
qc_only         : $params.qc_only
gvcf_only       : $params.gvcf_only
filter_only     : $params.filter_only
trim_adapter    : $params.trim_adapter
genomes_dir     : $params.genomes_dir
output_zarr     : $params.output_zarr
sif_path        : $params.sif_path
ploidy          : $params.ploidy
"""

// workflows 
include { QC } from './workflows/qc.nf'
include { GVCF } from './workflows/gvcf.nf'
include { FILTER } from './workflows/filter.nf'

workflow {
    if(params.qc_only && params.gvcf_only){
        // check parameters
        error "Error: only one of (qc_only, gvcf_only, filter_only) can be enabled."
    } else if (params.qc_only){
        // qc only
        QC()
    } else if (params.gvcf_only) {
        // gvcf only
        pf_bam_ch = Channel.fromFilePairs(params.bams, checkIfExists: true).map{index, bam_index -> [index, *bam_index.flatten()]}
        GVCF(pf_bam_ch)
    } else if (params.filter_only) {
        // filter only - expects gVCF files
        // Input format: sample_id, chr, gvcf, gvcf_idx
        // Using {,.idx} glob pattern to pair each gVCF with its index
        gvcf_pattern = "${params.inputdir}/*.chr*.g.vcf{,.idx}"
        gvcf_ch = Channel.fromFilePairs(gvcf_pattern, size: 2, flat: true, checkIfExists: true)
            .map { key, gvcf, idx ->
                def matcher = key =~ /(.+)\.chr(\d+)/
                if (matcher.matches()) {
                    def sample_id = matcher[0][1]
                    def chr = matcher[0][2]
                    tuple(sample_id, chr, gvcf, idx)
                }
            }
            .filter { it != null }
        
        FILTER(gvcf_ch)
    }
    else {
        // full pipeline: qc -> gvcf -> filter
        QC()
        gvcf_out = GVCF(QC.out)
        
        // Transform gVCF output for FILTER workflow
        // gvcf_out emits: tuple(gvcf, gvcf_idx, log)
        // We need to parse the filenames to extract sample_id and chromosome
        gvcf_for_filter = gvcf_out
            .flatMap { gvcf, gvcf_idx, log ->
                def filename = gvcf.name
                def matcher = filename =~ /(.+)\.chr(\d+)\.g\.vcf/
                if (matcher.matches()) {
                    def sample_id = matcher[0][1]
                    def chr = matcher[0][2]
                    [[sample_id, chr, gvcf, gvcf_idx]]
                } else {
                    []
                }
            }
        FILTER(gvcf_for_filter)
    }    
}