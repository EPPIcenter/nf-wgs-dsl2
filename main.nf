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
vqsr_only       : $params.vqsr_only
trim_adapter    : $params.trim_adapter
genomes_dir     : $params.genomes_dir
cross_vcfs_dir  : $params.cross_vcfs_dir
sif_path        : $params.sif_path
"""

// workflows 
include { QC } from './workflows/qc.nf'
include { GVCF } from './workflows/gvcf.nf'
include { VQSR } from './workflows/vqsr.nf'

workflow {
    if(params.qc_only && params.gvcf_only){
        // check parameters
        error "Error: only one of (qc_only, gvcf_only, vqsr_only) can be enabled."
    } else if (params.qc_only){
        // qc only
        QC()
    } else if (params.gvcf_only) {
        // gvcf only
        pf_bam_ch = Channel.fromFilePairs(params.bams, checkIfExists: true).map{index, bam_index -> [index, *bam_index.flatten()]}
        GVCF(pf_bam_ch)
    } else if (params.vqsr_only) {
        // vqsr only - expects gVCF files
        // Input format: sample_id, chr, gvcf, gvcf_idx
        gvcf_pattern = "${params.inputdir}/*.chr*.g.vcf"
        gvcf_ch = Channel.fromPath(gvcf_pattern)
            .map { file ->
                def filename = file.name
                def matcher = filename =~ /(.+)\.chr(\d+)\.g\.vcf/
                if (matcher.matches()) {
                    def sample_id = matcher[0][1]
                    def chr = matcher[0][2]
                    def idx_file = file.parent.resolve("${filename}.idx")
                    tuple(sample_id, chr, file, idx_file)
                }
            }
        VQSR(gvcf_ch)
    }
    else {
        // full pipeline: qc -> gvcf -> vqsr
        QC()
        gvcf_out = GVCF(QC.out)
        
        // Transform gVCF output for VQSR workflow
        // gvcf_out emits: tuple(gvcf, gvcf_idx, log)
        // We need to parse the filenames to extract sample_id and chromosome
        gvcf_for_vqsr = gvcf_out
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
        VQSR(gvcf_for_vqsr)
    }    
}