################################################################################
#######################           RiboDoc           ############################
################################################################################

################################################################################
### Tool for Ribosome Profiling Analysis
### Genome, Structure and Translation (GST) team - I2BC
################################################################################

if [ $# -eq 0 ]|| [ $# -eq 1 ]; then
    echo "PLEASE SPECIFY AVAILABLE RESOURCES BY ADDING CPUs AND MEMORY AMOUNT (in GB)"
else
    if [ "$1" -eq "1" ]; then
        cpu_use=$(($1/2))
    else
        cpu_use=$1
    fi
    available_memory_mb=$(($2*1000))
    echo "Maximum RAM used for the analysis : ${available_memory_mb}Mb"
    echo "Maximum CPU used for the analysis : "$((${cpu_use}*2))

    used_memory="--resources mem_mb=${available_memory_mb}"

    # conda list;
    snakemake -s /RiboDoc/RiboDoc/Snakefile --rerun-incomplete -j --dag -np --forceall | dot -Tsvg > /data/dag_all.svg
    snakemake -s /RiboDoc/RiboDoc/Snakefile -k --rerun-incomplete -j ${cpu_use}
fi
