# commands for building the docker image locally and converting it to singularity for use on the Crick HPC

# locally - ensure Docker Desktop is running and you have space! 
Docker build --platform linux/amd64 -t nf-wgs-dsl2:local -f conf/Dockerfile .
docker save nf-wgs-dsl2:local | gzip > ~/Desktop/nf-wgs-dsl2.tar.gz
scp ~/Desktop/nf-wgs-dsl2.tar.gz whittog@babs002.nemo.thecrick.org:/nemo/stp/babs/working/whittog/pipelines/nf-wgs-pf-dsl2/conf/

# ssh into babs002 and convert to singularity
ssh whittog@babs002.nemo.thecrick.org

cd /nemo/stp/babs/working/whittog/pipelines/nf-wgs-pf-dsl2/conf/
gunzip nf-wgs-dsl2.tar.gz

module load Singularity/3.11.3
export NXF_SINGULARITY_CACHEDIR="$PWD/singularity_cache"

singularity build nf-wgs-dsl2.sif docker-archive://nf-wgs-dsl2.tar
rm nf-wgs-dsl2.tar 