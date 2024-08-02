---
layout: default
title: Getting Started
nav_order: 2
has_children: false
---

# Getting Started

There are three ways to download or access the wgs nextflow pipeline:


### 1. Clone the Repository from GitHub:
You can clone the repository to your local machine using Git. Use the following command:

```
git clone https://github.com/EPPIcenter/nf-wgs-dsl2.git
```

### 2. Copy from Wynton:
If you have access to the Wynton server, you can copy the pipeline files directly. Use the following command to copy the files to your desired location:

```
cp /wynton/home/eppicenter/shared/WGS_pipeline_nextflow /path/to/destination
```

Replace /path/to/destination with your desired local path.

### 3. Point to the Repository on Wynton:
If you are running the pipeline directly on Wynton, you can simply point to the repository's location without copying the files. This allows you to access and execute the pipeline directly from the server.


## Reference Genome Files

There are reference genome files which are required for running the wgs nextflow pipeline. 

If you clone the wgs nextflow pipeline repository from Github, you will need to separately download the reference genome files. 

You can copy these reference genome files from `/wynton/home/eppicenter/shared/WGS_pipeline_nextflow/refs/genomes`. 

Please store the `genomes` directory inside `refs` directory. 


## Pipeline Overview
- `main.nf`: WGS workflow 
- `nextflow.config`: config file
- `workflows` 
  - `qc.nf`: QC sub-workflow 
  - `gvcf.nf`: GVCF sub-workflow
- `config`
  - `Apptainer`: file used to build nf-wgs-dsl2.sif  
  - `Dockerfile`: file for building docker image 
  - `base.config`: base config file 
  - `envs`: conda envs (under construction  🚧)
- `refs`: reference files used by both `QC_workflow` and `gVCF_workflow`
  - `adapters`: folder containing trimmomatic adapter files
  - `genomes`: reference genome files and more
  - `run_quality_report.Rmd`: r script for quality report used in `QC_workflow`


## Software Dependencies 

The wgs nextflow pipeline uses [nextflow](https://www.nextflow.io/) and will need to be installed prior to using the pipeline. Information about how to [install](https://www.nextflow.io/) and use the [command line tool](https://www.nextflow.io/docs/latest/cli.html) can be found on their [website](https://www.nextflow.io/). The tool is also available from other package managers such as [conda](https://anaconda.org/bioconda/nextflow) if you would like an alternative installation pathway. 


{: .note }
Nextflow requires the Java 11 (or higher) Runtime Environment.
