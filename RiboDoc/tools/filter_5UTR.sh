#!/bin/bash

# This script only saves transcripts containing 5-prime UTRs in GFF for metaprofiles, as they give a better view of the periodicity in complexe genomes
# Command example : > bash filter_5UTR.sh -g RESULTS/annex_database/NamedCDS_Homo_sapiens.GRCh38.104.gff3 -o Homo_sapiens.GRCh38.104.5UTR_genes_only.gff -r mRNA -u five_prime_UTR -t 0.25 -p RESULTS/annex_database/

usage() { echo "Usage: $0 -g <Input GFF> -o <Output GFF> -p <Path to features output folder> -r <transcript name in GFF> -u <5-prime-UTR name in GFF> -t <5-prime-UTR threshold>" 1>&2 ; echo -e "\n -g\tInput GFF\n -o\tPath and name for the output GFF\n -p\tPath to output folder for counts of each gff features\n -r\tFeature name for a transcript in GFF file\n -u\tFeature name for 5-prime-UTR name in GFF\n -t\tMinimum proportion of 5-prime-UTR lines compared to transcript lines\n" ; exit 1; }

#init variables
while getopts ":g:o:p:r:u:t:" option; do
    case "${option}" in
        g)
            g="${OPTARG}"
            ;;
        o)
            o="${OPTARG}"
            ;;
        p)
			p="${OPTARG}"
			;;
		r)
			r="${OPTARG}"
			;;
		u)
			u="${OPTARG}"
			;;
		t)
			t="${OPTARG}"
			;;
        *)
            usage
            ;;
    esac
done

#testing if arguments are non-empty
shift $((OPTIND-1))
if [ -z "${g}" ] || [ -z "${o}" ] || [ -z "${p}" ]; then
    usage
fi;

if [ -z "${r}" ]; then
    rna="mRNA"
fi;
if [ -z "${u}" ]; then
    utr="five_prime_UTR"
fi;
if [ -z "${t}" ]; then
    t="0.25";
fi;

# Look for the proportion of each feature in the GFF
awk '{printf("%s\n",$3)}' $g | sort | uniq -c > "${p}gff_features_counts.txt";

# Find the number of 5UTR and transcripts in the GFF
transcript_nbr=$(grep "$r" "${p}gff_features_counts.txt" | awk '{printf("%s\n",$1)}');
utr_nbr=$(grep "$u" "${p}gff_features_counts.txt" | awk '{printf("%s\n",$1)}');

# Calculate the minimum number of 5UTR lines compared to transcript lines to select only transcript with 5UTR
min_utr=$(echo "${transcript_nbr}*${t}" | bc);

# If there is a good proportion of 5-prime UTR features, filter transcript for metaprofiles
if [ ${utr_nbr} -ge $(printf "%.0f\n" ${min_utr/./,}) ]; then
    awk '
        BEGIN {utr="false"; i=0;}
        /^#/ {print}
        $3 ~ /^mRNA$|^gene$/ {
            if(utr=="true") {
                for(var in lines) {
                    print(lines[var]);
                }
                utr="false";
            }
            for(var in lines) {
                delete lines[var]
            }
            i=0;
        }
        {lines[i]=$0;i++;}
        $3=="five_prime_UTR" {utr="true";}
        END {
            if(utr=="true") {
                for(var in lines) {print(lines[var]);}
            }
        }
    ' "$g" > "$o";
else
    # If there are not enough lines defining a 5UTR, every transcript in the analysis is used
    cat "$g" > "$o";
fi;
