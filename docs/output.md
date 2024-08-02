---
layout: default
title: Pipeline Output
nav_order: 6
has_children: false
---

# QC Workflow Output 

The QC workflow `--inputdir` is a directory of fastq files. 

After running the Pipeline the `outputdir` will have the following output: 

- `execution_report.html`: nextflow execution report 
- `execution_timeline.html`: nextflow timeline 
- `execution_trace.txt`: nextflow trace 
- **final_bams**: for each sample we will have the following files
  - `sample_dup_metrics.txt`: metrics file
  - `sample.sorted.dup.bam`: the mapped cleaned and sorted bam file 
  - `sample.sorted.dup.pf.bam`: the final pf bam file 
  - `sample.sorted.dup.pf.bam.bai`: index for pf bam file (this file is need if later running PathWeaver)
  - `sample.sorted.dup.pf.bam.csi`: index for pf bam file 
- **final_qc_reports**: contains summary report files for all samples
  - `Bam_stats_hs_Final.tsv`: stats for human bams
  - `Bam_stats_pf_Final.tsv`: stats for pf bams
  - `Bam_stats_total_Final.tsv`: stats for both pf and human bams
  - `InsertSize_Final.tsv`: insert size information for all samples
  - `multiqc_report.html`: multiqc/fastqc reports for all samples
  - `Ratios_hs_pf_reads.tsv`: table of ratios of human to pf reads
  - `ReadCoverage_final.tsv`: table of read coverage for each sample and chromosome
  - `run_quality_report.html`: a .Rmd script output summarizing run quality 
- **histograms**: for each sample we will have the following files
  - `sample_histo.pdf`: histogram intermediate file output
- **insert_size_metrics**: for each sample we will have the following files
  - `S3_WGS-batchC_4051126659-ST136_S5.insert2.txt`: insert size output
  - `S3_WGS-batchC_4051126659-ST136_S5.insert.txt`: insert size raw output


# GVCF Workflow Output 

The GVCF workflow `--inputdir` is the directory `final_bams` which is an output of the QC workflow. Most important are the `.dup.pf.bam` and `.dup.pf.bam.csi` files. 

After running the Pipeline the `outputdir` will have the following structure: 

- **chr1**: for each sample we will have the following files*
  - `sample.chr1.g.vcf`
  - `sample.chr1.g.vcf.idx`
  - `sample.chr1_log.txt`
- **chr11**: *
- **chr12**: *
- **chr13**: *
- **chr14**: *
- **chr2**: *
- **chr3**: *
- **chr4**: *
- **chr5**: *
- **chr6**: *
- **chr7**: *
- **chr8**: *
- **chr9**: *
- **chr14**: *
- `execution_report.html`: nextflow execution report 
- `execution_timeline.html`: nextflow timeline 
- `execution_trace.txt`: nextflow trace 





