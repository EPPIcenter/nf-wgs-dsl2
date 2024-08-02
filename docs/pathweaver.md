---
layout: default
title: Running PathWeaver
nav_order: 11
has_children: false
---

# Running PathWeaver

After running the WGS nextflow pipeline or the QC workflow, you may want to run PathWeaver to interpret your results. 

## Set up 

To set up Pathweaver, copy the PathWeaver workflow from Arwen `/data/general/finterly/pathweaver_run` to the desired location. 

## PathWeaver Input

### 1. Run Name
The PathWeaver workflow will create a new subdirectory in the `pathweaver_run` directory with the given "Run Name" and output all results to that new subdirectory. 

### 2. Bam Files Directory
PathWeaver works with the QC outputs from `final_bams` directory: 
- `.dup.pf.bam` : the final pf bam file 
- `.dup.pf.bam.bai` : index for pf bam file 

Copy these files to a single directory 

### 3. Genes List
PathWeaver will run on a given genes list, simple text file list, example below: 

`gene_list.txt`
```
PF3D7_0304600
PF3D7_0315200
PF3D7_1335900
```

## Running PathWeaver

The following command will run PathWeaver

```
bash final_run_pathweaver_2024_07_30.sh <run_name> <bam_dir> <gene_list> 
```

For example: 
```
bash final_run_pathweaver_2024_07_30.sh \
test_run \
path/final_bam_dir \ 
gene_list.txt  
```

Recommend running PathWeaver on Arwen using tmux. A log will print out as PathWeaver runs.


## Running PathWeaver Analysis / Plots 

To generate plots for your pathweaver output 

1. open the `generate_plots_2024_07_30.Rmd` in RStudio and edit the `run_dir` params to point to your run output directory: 

```
---
title: "R Notebook"
output: 
params:
  run_dir: "/home/me/pathweaver_run/test_run"
---
```

2. simply run the `generate_plots_2024_07_30.Rmd` file. Some adjustments may be needed for different analysis!