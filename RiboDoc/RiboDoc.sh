################################################################################
#######################           RiboDoc           ############################
################################################################################

################################################################################
### Tool for Ribosome Profiling
### GST team
################################################################################

# If the user does not provide a minimum number of cpu available :
# a quarter of computer resources is assigned for the analysis (doubled on specific steps)
used_memory=""
if [ $# -eq 0 ];
then
    cpu=`nproc --all`;
    echo "Total CPU available : "${cpu};
    cpu_use=$((${cpu}/4));
    available_memory_mb=""
else
    cpu_use=$(($1/2));
    available_memory_mb=$(($2*1000))
    echo "Maximum RAM used for the analysis : ${available_memory_mb}Mb";
fi;
echo "Maximum CPU used for the analysis : "$((${cpu_use}*2));

used_memory="--resources mem_mb=${available_memory_mb}"

# conda list;
snakemake -s /work/RiboDoc/Snakefile -j --dag -np | dot -Tsvg > /data/dag_last-run.svg;
snakemake -s /work/RiboDoc/Snakefile -j --dag -np --forceall | dot -Tsvg > /data/dag_all.svg;
snakemake -s /work/RiboDoc/Snakefile -j ${cpu_use} -k ${used_memory};
# snakemake -s /ORFribo_RiboDoc/RiboDoc_BAM2Reads/Snakefile -j ${cpu_use} -k --resources mem_mb=${used_memory};
