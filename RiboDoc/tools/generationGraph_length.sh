#! /bin/bash

usage() { echo "Usage: $0 -p <Path to working directory> -P <Path to scripts> -N <Sample name>" 1>&2 ; echo -e "\n -p\tPath to working directory\n -P\tPath to scripts -N\tsample name\n" ; exit 1; }

while getopts ":p:P:N:" option; do
    case "${option}" in
        p)
            p=${OPTARG}
            ;;
        P)
            P=${OPTARG}
            ;;
        N)
            N=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [ -z "${p}" ] || [ -z "${P}" ] || [ -z "${N}" ]; then
    usage
fi;

########### Reads Length Repartition ###########
Rscript "${P}generationGraph_length.R" --work_dir "${p}" --sample_name "${N}";
