################################################################################
#######################           RiboDoc           ############################
################################################################################

################################################################################
### Tool for Ribosome Profiling
### GST team
################################################################################

printf "RiboDoc version : 0.6.2\n"

# If the user does not provide a minimum number of cpu available :
# a quarter of computer resources is assigned for the analysis (doubled on specific steps)
if [ $# -eq 0 ];
then
    cpu=`nproc --all`;
    echo "Total CPU available : "${cpu};
    cpu_use=$((${cpu}/4));
else
    cpu_use=$(($1/2));
fi;
echo "Maximum CPU used for the analysis : "$((${cpu_use}*2));

available_memory_kb=`grep MemTotal /proc/meminfo | sed -E "s/^[^0-9]+([0-9]+)[^0-9]+$/\1/"`;
available_memory_mb=$((${available_memory_kb}/1000))
used_memory=$((${available_memory_mb}/4*3))

# conda list;
snakemake -s /work/RiboDoc/Snakefile -j --dag -np | dot -Tsvg > /data/dag_last-run.svg;
snakemake -s /work/RiboDoc/Snakefile -j --dag -np --forceall | dot -Tsvg > /data/dag_all.svg;
snakemake -s /work/RiboDoc/Snakefile -j ${cpu_use} -k --resources mem_mb=${used_memory};
