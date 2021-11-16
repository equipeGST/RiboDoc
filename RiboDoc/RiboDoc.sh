################################################################################
#######################           RiboDoc           ############################
################################################################################

################################################################################
### Tool for Ribosome Profiling
### GST team
################################################################################

cpu=`nproc --all`;
echo "Total CPU available : "${cpu};
cpu_use=$((${cpu}/4));
echo "Minimum CPU used for the analysis : "${cpu_use};

# conda list;
snakemake -s /RiboDoc/RiboDoc/Snakefile -j --dag -np | dot -Tsvg > /data/dag_last-run.svg;
snakemake -s /RiboDoc/RiboDoc/Snakefile -j --dag -np --forceall | dot -Tsvg > /data/dag_all.svg;
snakemake -s /RiboDoc/RiboDoc/Snakefile -j ${cpu_use} -k;
