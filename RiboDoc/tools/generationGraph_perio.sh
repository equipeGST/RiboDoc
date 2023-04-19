#! /bin/bash

usage() { echo "Usage: $0 -p <Path to working directory> -P <Path to scripts> -N <Sample name> -l <Read length> -m <before start position> -M <after start position> -t <feature type in GFF3 annotation>" 1>&2 ; echo -e "\n -p\tPath to working directory\n -P\tPath to scripts\n -N\tsample name\n -l\tRead length\n -m <int>\tnumber of base before start codon\n -M <int>\tnumber of base after start codon\n -t <string>\tSpecify feature type in GFF3 annotation" ; exit 1; }

while getopts ":p:P:N:l:m:M:t:" option; do
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
        l)
            l=${OPTARG}
            ;;
        m)
            m=${OPTARG}
            ;;
        M)
			M=${OPTARG}
			;;
		t)
			t=${OPTARG}
			;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [ -z "${p}" ] || [ -z "${P}" ] || [ -z "${N}" ] || [ -z "${l}" ] || [ -z "${m}" ] || [ -z "${M}" ] || [ -z "${t}" ]; then
    usage
fi;

########### Periodicity ###########
for pos in start stop
do

    if [ ${pos} = "start" ]; then
        file="${N}.${l}.periodicity.start.${t}.-${m}+${M}";
        Rscript "${P}generationGraph_perio.R" --work_dir "${p}" --name "${file}" --feature_type "${t}" --start_or_stop "${pos}" --length_in_utr "${m}";
    else
        file="${N}.${l}.periodicity.stop.${t}.-${M}+${m}";
        Rscript "${P}generationGraph_perio.R" --work_dir "${p}" --name "${file}" --feature_type "${t}" --start_or_stop "${pos}" --length_in_utr "${M}";
    fi;

done;
