---
layout: default
title: Getting Started
nav_order: 2
has_children: false
---

# Getting Started

There are three ways to download or access the wgs nextflow pipeline:


1. Clone the Repository from GitHub:
You can clone the repository to your local machine using Git. Use the following command:

```
git clone https://github.com/EPPIcenter/nf-wgs-dsl2.git
```

2. Copy from Wynton:
If you have access to the Wynton server, you can copy the pipeline files directly. Use the following command to copy the files to your desired location:

```
cp /wynton/home/eppicenter/shared/WGS_pipeline_nextflow /path/to/destination

```
Replace /path/to/destination with your desired local path.

3. Point to the Repository on Wynton:
If you are running the pipeline directly on Wynton, you can simply point to the repository's location without copying the files. This allows you to access and execute the pipeline directly from the server.


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
  - `envs`: conda envs (under construction :construction:)
- `refs`: reference files used by both `QC_workflow` and `gVCF_workflow`
  - `adapters`: folder containing trimmomatic adapter files
  - `genomes`: reference genome files and more
  - `run_quality_report.Rmd`: r script for quality report used in `QC_workflow`


## Reference Genome Files

There are reference genome files which are required for running the wgs nextflow pipeline. 

If you clone the wgs nextflow pipeline repository from Github, you will need to separately download the reference genome files. You can copy these reference genome files from `/wynton/home/eppicenter/shared/WGS_pipeline_nextflow/refs/genomes`. Please store the `genomes` directory inside `refs` directory. 


## Software Dependencies 

The wgs nextflow pipeline uses [nextflow](https://www.nextflow.io/) and will need to be installed prior to using the pipeline. Information about how to [install](https://www.nextflow.io/) and use the [command line tool](https://www.nextflow.io/docs/latest/cli.html) can be found on their [website](https://www.nextflow.io/). The tool is also available from other package managers such as [conda](https://anaconda.org/bioconda/nextflow) if you would like an alternative installation pathway. 


{: .note }
Nextflow requires the Java 11 (or higher) Runtime Environment.


## Runtime Profiles

Runtime profiles will provide all dependencies and setup needed for different computing environments. As an example, if you are using a cluster, grid or HPC environment, `apptainer` would be an appropriate profile as it supplies an image with all dependencies ready. If you are using a local computer, `docker` would be more appropriate. You can also choose to install the dependencies independently and run the pipeline that way if you choose to, but it is not recommended. 

Currently, [Sun of Grid Engine (SGE)](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html) is the only supported cluster environment.

The nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on.

### Apptainer

{: .note }
[Apptainer](https://github.com/apptainer/apptainer/releases) is a prerequisite.

Apptainer should be used if you are using a computing cluster or grid. You will first need to build the apptainer image before you can use the image (if it is not already built). 

To build the image, run the command below from the `conf` directory:

```bash
apptainer build nf-wgs-dsl2.sif Apptainer
```

And then include the `apptainer` profile on the command line. 

```bash
nextflow run main.nf \
-profile sge,apptainer \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory
```

{: .note }
You should also include the job scheduler you will be using. In this case, `sge` is the job scheduler that will be used. Contact your system administrator if you are unsure about this setting.

### Docker

{: .note }
[Docker](https://www.docker.com/) is a prerequisite.

The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

Follow the steps below to setup your docker image from the `conf` directory:

*Note: [docker](https://www.docker.com/) is a prerequisite.*

```bash
docker build -t eppicenter/nf-wgs-dsl2 .
```
And you're done! To run the pipeline, simply add `-profile docker` in your command. 

```bash
nextflow run main.nf \
-profile docker \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory
```

### Conda (less recommended, under construction :construction:)

To use conda, you must first install either [conda](https://docs.conda.io/en/latest/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html). Once installed, include the `conda` profile on the command line.

```bash
nextflow run main.nf \
-profile conda \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory
```