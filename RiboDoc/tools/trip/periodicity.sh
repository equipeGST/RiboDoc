#!/bin/bash

usage() { echo "Usage: $0 -N <sample name> -l <read length> -m <before start position> -M <after start position> -G <GFF3 file> -D <path to bed count> -p <position> -t <feature type in GFF3 annotation> -r <by row> -O <output dir>" 1>&2 ; echo "\n -G <string>\tgff3 file with all annotations\n -D <string>\tbed count directory path\n -N <string>\tSample name\n -l <string>\tRead length of interest\n -p <string>\tcodon of interest (start or stop)\n -m <int>\tnumber of base before start codon\n -M <int>\tnumber of base after start codon\n -t <string>\tSpecify feature type in GFF3 annotation. ('CDS','five_prime_UTR','three_prime_UTR')\n -r <string>\tTwo option : 'metagene' (default value) to have all GFF3 genes in the same table ; 'byrow' to have a table by GFF3 gene\n -O\tOutput directory" ; exit 1; }

while getopts ":G:D:m:M:p:t:N:l:r:O:" option; do
    case "${option}" in
        G)
            G=${OPTARG}
            ;;
        D)
            D=${OPTARG}
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
		p)
			p=${OPTARG}
			;;
		t)
			t=${OPTARG}
			;;
		r)
			r=${OPTARG}
			;;
		O)
			O=${OPTARG}
			;;
        *)
            usage
            ;;
    esac
done

#testing if arguments are non-empty
shift $((OPTIND-1))
if [ -z "${G}" ] || [ -z "${D}" ] || [ -z "${m}" ] || [ -z "${M}" ] || [ -z "${p}" ] || [ -z "${t}" ] || [ -z "${N}" ] || [ -z "${l}" ]; then
    usage
fi
if [ -z "${r}" ]; then
	r="metagene"
fi
if [ -z "${O}" ]; then
	O="."
fi

#verification GFF only
testGFF=`grep "=" ${G}`
if [ -z "${testGFF}" ] ; then
	echo "Error, GFF format needed, not GTF or other";
	exit
fi


##Soft GFF3 - tempo file
#depend of option p
if [ "${p}" = "start" ]; then
	grep -v "#" ${G} | gawk '{if($3=="'${t}'" && $4>='"${m}"'){if($7=="+"){print $1"\t"$4-'"${m}"'-1"\t"$4+'"${M}"'"\t"$7}else if($7=="-"){print $1"\t"$5+'"${m}"'"\t"$5-'"${M}"'-1"\t"$7}}}' > "${O}${N}.${l}.databaseTempo.periodicity.txt";
elif [ "${p}" = "stop" ]; then
    grep -v "#" ${G} | gawk '{if($3=="'${t}'"){if($7=="+"){print $1"\t"$5-'"${m}"'-1"\t"$5+'"${M}"'"\t"$7}else if($7=="-"){print $1"\t"$4+'"${m}"'"\t"$4-'"${M}"'-1"\t"$7}}}' > "${O}${N}.${l}.databaseTempo.periodicity.txt";
fi

if [ "${r}" = "metagene" ]; then
    #Periodicity calcul
	for samp in `ls ${D} | grep "${N}.${l}" | grep "bed"`;
	do
		name="${samp%.*.*}"
		sample_temp=$(mktemp /tmp/sampleTempo.periodicity.XXX);
		echo "Sample "${samp}
		gawk '{if($6=="+"){print $4"\t"$1":"$2":+"}else if($6=="-"){print $4"\t"$1":"$3":-"}}' ${D}${samp} | sort -k2 > "${sample_temp}" ;
		total=$((${m}+${M})) ;
		for pos in `seq 0 ${total}`;
		do
			data_temp_plus=$(mktemp /tmp/dataplusTempo.periodicity.XXX);
			data_temp_minus=$(mktemp /tmp/dataminusTempo.periodicity.XXX);
			gawk -v num=${pos} '{if($4=="+"){print $1":"$2+num":"$4}}' "${O}${N}.${l}.databaseTempo.periodicity.txt" | sort > ${data_temp_plus} ;
			join -1 2 -2 1 -o 1.1 ${sample_temp} ${data_temp_plus} | gawk -v num=${pos} 'BEGIN{SUM=0}{SUM+=$1}END{print num"\t"SUM}' >> "${O}periodicity/${name}.periodicity.${p}.${t}.-${m}+${M}.+.txt" ;
			gawk -v num=${pos} '{if($4=="-"){print $1":"$2-num":"$4}}' "${O}${N}.${l}.databaseTempo.periodicity.txt" | sort > ${data_temp_minus} ;
			join -1 2 -2 1 -o 1.1 ${sample_temp} ${data_temp_minus} | gawk -v num=${pos} 'BEGIN{SUM=0}{SUM+=$1}END{print num"\t"SUM}' >> "${O}periodicity/${name}.periodicity.${p}.${t}.-${m}+${M}.-.txt" ;
			rm -f ${data_temp_plus} ${data_temp_minus}
		done;
		join -1 1 -2 1 -o 1.1,1.2,2.2 "${O}periodicity/${name}.periodicity.${p}.${t}.-${m}+${M}.+.txt" "${O}periodicity/${name}.periodicity.${p}.${t}.-${m}+${M}.-.txt" | awk 'BEGIN{FS==" "}{print $1"\t"$2+$3}' > "${O}periodicity/${name}.periodicity.${p}.${t}.-${m}+${M}.txt";
		rm -f "${O}periodicity/${name}.periodicity.${p}.${t}.-${m}+${M}.+.txt" "${O}periodicity/${name}.periodicity.${p}.${t}.-${m}+${M}.-.txt";
	done;
fi;

rm -f "${O}${N}.${l}.databaseTempo.periodicity.txt" "${O}${N}.${l}.ligneDatabaseTempo.periodicity.txt" "${sample_temp}";
