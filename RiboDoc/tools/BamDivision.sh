#! /bin/bash

###
###This script allowed the division of SAM file according to read length and turns it into BAM
###

usage() { echo "Usage: $0 -l <Read length> -B <BAM file> -T <Threads> -N <sample name> -o <output dir>" 1>&2 ; echo -e "\n -B\tBAM files (with path) without multimapping and unalign reads\n -l\tRead length\n -T\tThreads number allows (default value : 4)\n -N\tSample name\n -o\tOutput directory" ; exit 1; }

#init variables
while getopts ":B:T:N:l:o:" option; do
    case "${option}" in
        B)
            B=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
		T)
			T=${OPTARG}
			;;
		N)
			N=${OPTARG}
			;;
		o)
			o=${OPTARG}
			;;
        *)
            usage
            ;;
    esac
done

#testing if arguments are non-empty
shift $((OPTIND-1))
if [ -z "${B}" ] || [ -z "${l}" ] || [ -z "${N}" ] ; then
    usage
fi
#if threads empty, default value = 4
if [ -z "${T}" ]; then
	T=4
fi
#if output directory is empty, write in working directory
if [ -z "${o}" ]; then
	o="."
fi

# Output directory
mkdir -p ${o}"bamDivision"

#Bam Division
# Division of reads depending of read lengths
gawk -v l=${l} '{if($1 !~ /^ *@/){if(length($10)==l){print $0}}else{print $0}}' $(samtools view -@ ${T} -h ${B}) > ${o}bamDivision/${N}.${l}.sam;

#Conversion as bam
samtools view -@ ${T} -b ${o}bamDivision/${N}.${l}.sam > ${o}bamDivision/${N}.${l}.uniq.bam ;
rm -rf ${o}bamDivision/${N}.${l}.sam
# Sort bam file
samtools sort -@ ${T} ${o}bamDivision/${N}.${l}.uniq.bam -o ${o}bamDivision/${N}.${l}.uniq.sort.bam ;
rm -rf ${o}bamDivision/${N}.${l}.uniq.bam
# Index bam file
samtools index ${o}bamDivision/${N}.${l}.uniq.sort.bam;
