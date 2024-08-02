---
layout: default
title: Running the Pipeline
nav_order: 3
has_children: false
---

# Running the Pipeline

## A Basic Example

Below is an example of how to run the pipeline setting only the most necessary parameters on an HPC.

```bash
nextflow run main.nf \
-profile sge,apptainer \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory
```



## Runtime Profiles

How you run the Nextflow pipeline depends on your environment (locally or on a cluster). Due to the large size of whole genome sequencing data, it is recommended to run the pipeline on a cluster like Wynton.

Runtime profiles will provide all dependencies and setup needed for different computing environments. 

- If you are using a cluster, grid or HPC environment, `apptainer` would be an appropriate profile as it supplies an image with all dependencies ready. 
- If you are using a local computer, `docker` would be more appropriate. You can also choose to install the dependencies independently and run the pipeline that way if you choose to, but it is not recommended. 

Currently, [Sun of Grid Engine (SGE)](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html) is the only supported cluster environment.

The nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on.

### Apptainer (cluster)

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

### Docker (local)

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

### Conda (less recommended, under construction 🚧)

To use conda, you must first install either [conda](https://docs.conda.io/en/latest/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html). Once installed, include the `conda` profile on the command line.

```bash
nextflow run main.nf \
-profile conda \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory
```


