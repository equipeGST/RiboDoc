#!/bin/bash

# This script only saves transcripts containing 5-prime UTRs in GFF for metaprofiles, as they give a better view of the periodicity in complexe genomes
# Command example : > bash filter_5UTR.sh -g RESULTS/annex_database/NamedCDS_Homo_sapiens.GRCh38.104.gff3 -o Homo_sapiens.GRCh38.104.5UTR_genes_only.gff -r mRNA -u five_prime_UTR -t 0.25 -p RESULTS/annex_database/

Help()
{
    echo "This script filters transcripts containing 5-prime UTRs in GFF for metaprofiles, as they give a better view of the periodicity in complexe genomes"
    echo
    echo "report_table [-g|h|p|r|o|t|u]"
    echo "options:"
    echo "-g GFF            --gff           Input GFF File"
    echo "-h HELP           --help          Print this Help."
    echo "-p PATH           --path          Path to output folder for counts of each gff feature"
    echo "-r TRANSCRIPT     --transcript    Feature name for a transcript in GFF file. Default : 'mRNA'."
    echo "-o OUTPUT         --output        Name of output file"
    echo "-t THRESHOLD      --threshold     Minimum proportion of 5-prime-UTR lines compared to transcript lines. Default : 0.25"
    echo "-u UTR            --utr           Feature name for 5-prime-UTR name in GFF. Default : 'five_prime_UTR'"
    echo

    exit 1
}

while getopts ":g:o:p:r:u:t:" option; do
    case "${option}" in
        g) # Path to the GFF file
            g=${OPTARG}
            ;;
        o) # Path to the output file
            o=${OPTARG}
            ;;
        p) # Path to output folder
			p=${OPTARG}
			;;
		r) # Feature name for a transcript in GFF file. 
			r=${OPTARG}
			;;
		u) # 5-prime UTR threshold 
			u=${OPTARG}
			;;
		t) # 5-prime UTR name in GFF 
			t=${OPTARG}
			;;
        h) # display Help
            Help
            ;;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit
            ;;
    esac
done

if [ -z "${g}" ] || [ -z "${o}" ] || [ -z "${p}" ]; then
    Help
fi;

if [ -z "${r}" ]; then
    rna="mRNA"
fi;

if [ -z "${u}" ]; then
    utr="five_prime_UTR"
fi;

if [ -z "${t}" ]; then
    t="0.25"
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
        # Save header
        /^#/ {print}
        $3 ~ /^mRNA$|^gene$/ {
            # If a 5UTR has been found for previous gene/transcript
            if(utr=="true") {
                # Print all corresponding lines
                for(var in lines) {
                    print(lines[var]);
                }
                utr="false";
            }
            # Delete every line corresponding to the previous gene/transcript
            for(var in lines) {
                delete lines[var]
            }
            i=0;
        }
        # Save current line
        {lines[i]=$0;i++;}
        # If a 5UTR line is present, change flag to "true"
        $3=="five_prime_UTR" {utr="true";}
        END {
            if(utr=="true") {
                for(var in lines) {print(lines[var]);}
            }
        }
    ' "$g" > "$o";
else
    # If there are not enough lines defining a 5UTR (usually in simple genomes), every transcript in the analysis is used
    cat "$g" > "$o";
fi;
