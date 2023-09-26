################################################################################
#######################           RiboDoc           ############################
################################################################################

################################################################################
### Tool for Ribosome Profiling Analysis
### Genome, Structure and Translation (GST) team - I2BC
################################################################################

if [ $# -eq 0 ]|| [ $# -eq 1 ]; then
    printf "PLEASE SPECIFY AVAILABLE RESOURCES BY ADDING CPUs AND MEMORY AMOUNT (in GB)\n"
    printf "Example (with 10 CPU and 30GB of memory) :\n\tsingularity run --bind /path/to/workdir/ ribodoc_v0.9.0.19.sif 10 30\n"
else
    if [ "$1" -eq "1" ]; then
        cpu_use=$(($1/2))
    else
        cpu_use=$1
    fi
    available_memory_mb=$(($2*1000))
    printf "Maximum RAM used for the analysis : %sMb\n" "${available_memory_mb}"
    printf "Maximum CPU used for the analysis : %s\n" $((${cpu_use}*2))

    used_memory="--resources mem_mb=${available_memory_mb}"

    # To be sure the hidden '.snakemake/' folder is in the working directory, snakemake must be invoqued in '/data/'
    cd /data/
    snakemake -s /RiboDoc/RiboDoc/Snakefile --rerun-incomplete -j --dag -np --forceall | dot -Tsvg > /data/dag_all.svg
    snakemake -s /RiboDoc/RiboDoc/Snakefile -k --rerun-incomplete -j ${cpu_use}
fi
