---
layout: default
title: Pipeline Parameters
nav_order: 4
has_children: false
---

# Modifying Parameters

Pipeline specific parameters can be found in `nextflow.config` 

| Retuired Parameters | Description |
|---|---|
| inputdir | The folder that contains the input files (default 'data') |
| outputdir | The folder where you want the resulting data to be saved (default 'results/results') |

| Optional Parameters | Description |
|---|---|
| qc_only | If enabled, only QC workflow is run (default 'false') |
| gvcf_only | If enabled, only gVCF workflow is run (default 'false') |
| trim_adapter | The adapter used for initial trimming of reads (default 'NexteraPE-custom.fa') |


Some common nextflow parameters 

| Nextflow Parameters | Description |
|---|---|
| -profile | allows you to select a specific configuration profile defined in your Nextflow configuration files |
| -work-dir | specifies the directory where Nextflow will store intermediate files and results |
| -config | allows you to specify an additional configuration file |


Notice how there some parameters have two hypens in front of their name (ie. `--inputdir`) and some only have one (ie. `-profile`). The difference is because some parameters are *pipeline-specific* while others are defined by nextflow. Two hyphens ('--') indcate that the parameter is defined within the pipeline, and a single hypen (`-`) indicates that the parameter is nextflow defined. 


## Example

Below is an example of how you may use the above parameters on the command line:

Running QC Only with different trim adapters (default is a custom Nextera adapter) by adding --qc_only and --trim_adapter parameters as well as nextflow -workdir parameter
```
nextflow run main.nf \
-work-dir path/work \
-profile sge,apptainer \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory \
--trim_adapter path/adapters/TruSeq2-PE.fa 
--qc_only true
```






