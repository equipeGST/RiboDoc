#!/bin/bash

usage() { echo "Usage: $0 -l <Read length> -D <path to BAM by length> -N <Sample name> -F <Fasta file> -O <output dir>" 1>&2 ; echo "\n -N\tsample name\n -D\tpath to directory with uniq and sort BAM files\n -l\tRead length\n -F\tGenome sequence in fasta format\n -O\tOutput directory" ; exit 1; }

while getopts ":D:N:l:F:O:" option; do
    case "${option}" in
        D)
            D=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
        N)
            N=${OPTARG}
            ;;
        F)
            F=${OPTARG}
            ;;
        O)
            O=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [ -z "${D}" ] || [ -z "${l}" ] || [ -z "${N}" ] || [ -z "${F}" ]; then
    usage
fi
if [ -z "${O}" ]; then
    O="."
fi

#Read length repartition

echo "Analyse read length " ${l};
bamToBed -i ${D}${N}.${l}.uniq.sort.bam > ${O}readsLengthRepartition/${N}.${l}.bed ;
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' ${O}readsLengthRepartition/${N}.${l}.bed | uniq -c | tr -s ' ' | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$2":"$3"-"$4"("$5")\t"$5}' | sort -k 5 > ${O}bedCount/${N}.${l}.count.bed;
bedtools getfasta -fi ${F} -bed ${O}bedCount/${N}.${l}.count.bed -s -tab -fullHeader -fo ${O}bedCount/${N}.${l}.fa;
uniq ${O}bedCount/${N}.${l}.fa > ${O}bedCount/${N}.${l}.uniq.fa;
join -1 5 -2 1 -o 1.1,1.2,1.3,1.4,2.2,1.6 ${O}bedCount/${N}.${l}.count.bed ${O}bedCount/${N}.${l}.uniq.fa | tr " " "\t" > ${O}sequenceBedCount/${N}.${l}.count.sequence.bed;
rm -f ${O}bedCount/${N}.${l}.fa ${O}bedCount/${N}.${l}.uniq.fa;
