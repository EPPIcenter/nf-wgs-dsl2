---
layout: default
title: Pipeline Resources
nav_order: 5
has_children: false
---

# Pipeline resources

Each process should be assigned a specific amount of RAM, processors and execution time. There is no one size fits all per se, but a good strategy to begin is requesting more resources than you think you will need, and then reduce as needed.

There are multiple resource tiers defined in `conf/base.config` and each module is assigned one of these tiers. You may update these tiers if you wish to modify resources broadly. The resources are described as follows:

|Resource|Description|
|:--:|:--:|
|time|The amount of time the process should use. After surpassing this time, the job will timeout.|
|cpus|The number of cores to use. A core is a processing unit that executes commands in code. Cores are crucial for parallel processing.|
|penv|The parallel environment to use. You will generally never need to change this value.|
|memory|The amount of memory to use for the process|

If you look you will notice that each process is assigned different resources: 

For example bwa_align is one of the more memory intensive processes and uses 8 cpus, 16 GB of memory, and up to 24 hours. If a task attempt fails, the resources are increased. 
```
 withName:bwa_align {
        cpus = { check_max( 8 * task.attempt, 'cpus') }
        memory = { check_max( 16.GB * task.attempt, 'memory' ) }
        time = { check_max( 24.h * task.attempt, 'time') }    
    }
```