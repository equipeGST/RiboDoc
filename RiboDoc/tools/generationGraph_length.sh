#! /bin/bash

usage() { echo "Usage: $0 -p <Path to working directory> -N <Sample name>" 1>&2 ; echo -e "\n -p\tPath to working directory\n -N\tsample name\n" ; exit 1; }

while getopts ":p:N:" option; do
    case "${option}" in
        p)
            p=${OPTARG}
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
if [ -z "${p}" ] || [ -z "${N}" ]; then
    usage
fi;

path="${p}RESULTS/qualitativeAnalysis/";

########### Reads Length Repartition ###########
echo "#! bin/R" > "${path}graphes/readsLengthRepartition/${N}.tempoR.R";
echo "readsLength<-read.table(file = '"${path}"readsLengthRepartition/"${N}".readsLengthRepartition.txt')" >> "${path}graphes/readsLengthRepartition/${N}.tempoR.R";
echo "jpeg(filename = '"${path}"graphes/readsLengthRepartition/"${N}".readsLengthRepartition.jpeg')" >> "${path}graphes/readsLengthRepartition/${N}.tempoR.R";
echo "barplot(readsLength\$V2,names.arg = readsLength\$V1, xlab = 'Read lengths', ylab = 'Number of reads')" >> "${path}graphes/readsLengthRepartition/${N}.tempoR.R";
echo "dev.off()" >> "${path}graphes/readsLengthRepartition/${N}.tempoR.R";
R CMD BATCH "${path}graphes/readsLengthRepartition/${N}.tempoR.R";
rm -f "${path}graphes/readsLengthRepartition/${N}.tempoR.R";
