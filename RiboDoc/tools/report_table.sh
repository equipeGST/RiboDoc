#!/bin/bash

####################################################################################################
#                                                                                                  #
#    This script is used to generate a csv file with all important stats from RNAseq analysis.     #
#                                                                                                  #
####################################################################################################

# Usage :

    # bash report_to_array.sh -d DATE -l LIBRARY -p PROTOCOLE -s SEQUENCER -r INPUT FILES -o OUTPUT

# Example :

    # bash report_to_array.sh -d 2021-09-05 -p RiboSeq -s NextSeq -D /home/enora.corler/Enora/FTSJ1/ -r /home/enora.corler/Enora/FTSJ1/stats/ -o FTSJ1_analysis-summary.csv

Help()
{
    echo "This script is used to generate a csv file with all important stats from RNAseq analysis"
    echo
    echo "report_table [-d|D|h|l|p|r|s|o]"
    echo "options:"
    echo "-d DATE           --date          Date of Sequencing"
    echo "-D DIRECTORY      --directory     Path to working directoy"
    echo "-h HELP           --help          Print this Help."
    echo "-l LIBRARY        --library       Library Name"
    echo "-p PROTOCOLE      --protocole     Protocole (RNAseq, RiboSeq, ...)"
    echo "-r PATH           --input         Path to stats directory"
    echo "-o FILE           --output        Name of output file"
    echo "-s SEQUENCER      --sequencer     Instrument Name"
    echo

    exit 1
}

while getopts ":r:d:D:p:s:l:o:" option; do
    case "${option}" in
        r) # Path to input files
            stats=${OPTARG}
            ;;
        o) # Name of output file
            output=${OPTARG}
            ;;
        d) # Date
            date=${OPTARG}
            ;;
        D) # Working directory
            dir=${OPTARG}
            ;;
        p) # Protocole
            protocole=${OPTARG}
            ;;
        s) # Sequencer used
            seq=${OPTARG}
            ;;
        l) # Library used for sequencing
            library=${OPTARG}
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

if [ -z "${stats}" ] || [ -z "${output}" ]; then # If no input or output is given, display help
    Help
fi

if [ -z "${date}" ]; then # If no date is given, a blank is written
    date=" "
fi

if [ -z "${dir}" ]; then # If no working directory is given, the current directory is used
    dir="./"
fi

if [ -z "${library}" ]; then # If no library is given, a blank is written
    library=" "
fi

if [ -z "${seq}" ]; then # If no sequencer name is given, a blank is written
    seq=" "
fi

if [ -z "${protocole}" ]; then # If no protocole is given, RiboSeq is written
    protocole="RiboSeq"
fi

# Create output file with header
echo -e "Protocole\tDate\tLibrary\tSequencer\tSample\tTotal Read\tReads W adapters\t% of Total Reads\tToo short\t% of reads W adapters\tToo long\t% of reads W adapters\tPass filters\t% of Total Reads\tReads mapped on rRNA\tContamination (%)\tUniquely Mapped\t% of raw data\t% of reads W adapters\tMapped more than once\tUnmapped" > "${output}"

dirlist=$(ls ${dir}/fastq)

for file in ${dirlist}; do
    file_name="$(basename ${file})"
    sample="${file_name%.*.*}"

    adapt_trimming_log="${stats}/${sample}_adapt_trimming.log" # Path to adapt trimming log file
    outRNA_log="${stats}/${sample}_bowtie2_run_outRNA.log" # Path to outRNA log file
    mapping_bowtie2_log="${stats}/${sample}_run_mapping_bowtie2.log" # Path to mapping bowtie2 log file
    mapping_hisat2_log="${stats}/${sample}_run_mapping_hisat2.log" # Path to mapping hisat2 log file

    Total_reads=$(grep "Total reads" ${adapt_trimming_log} | awk -F ':' '{{print $2}}' | tr -d '[:space:],')
    Reads_with_adapters=$(grep "Reads with adapters" ${adapt_trimming_log} | grep -oP '\d+,\d+,\d+' | tr -d ',')
    Percent=$(grep "Reads with adapters" ${adapt_trimming_log} | grep -oP '\(\K[^\)]+')
    Too_short=$(grep "too short" ${adapt_trimming_log} | grep -oP '\d+,\d+,\d+' | tr -d ',')
    Percent_too_short=$(grep "too short" ${adapt_trimming_log} | grep -oP '\(\K[^\)]+')
    Too_long=$(grep "too long" ${adapt_trimming_log} | grep -oP '\d+,\d+,\d+' | tr -d ',')
    Percent_too_long=$(grep "too long" ${adapt_trimming_log} | grep -oP '\(\K[^\)]+')
    Passing_filters=$(grep "Reads written (passing filters)" ${adapt_trimming_log} | grep -oP '\d+,\d+,\d+' | tr -d ',')
    Percent_passing_filters=$(bc <<< "scale=1; $Passing_filters * 100 / $Total_reads")

    uniq_reads_rRNA=$(grep "aligned exactly 1 time" ${outRNA_log} | awk '{{print $1}}')
    mult_reads_rRNA=$(grep "aligned >1 times" ${outRNA_log} | awk '{{print $1}}')
    rRNA=$(( uniq_reads_rRNA + mult_reads_rRNA ))
    Contamination=$(bc <<< "scale=1; $rRNA * 100 / $Passing_filters")

    Unique_reads_hisat=$(grep "aligned exactly 1 time" $mapping_hisat2_log | awk '{{print $1}}')
    Unique_reads_bowtie=$(grep "aligned exactly 1 time" $mapping_bowtie2_log | awk '{{print $1}}')
    Uniquely_mapped_genome=$(( Unique_reads_hisat + Unique_reads_bowtie ))
    Multi_reads_hisat=$(grep "aligned >1 times" $mapping_hisat2_log | awk '{{print $1}}')
    Multi_reads_bowtie=$(grep "aligned >1 times" $mapping_bowtie2_log | awk '{{print $1}}')
    Multi_mapped=$(( Multi_reads_hisat + Multi_reads_bowtie ))
    Unmapped_reads_bowtie=$(grep "aligned 0 time" $mapping_bowtie2_log | awk '{{print $1}}')
    Unmapped_reads_hisat=$(grep "aligned 0 time" $mapping_hisat2_log | awk '{{print $1}}')
    Unmapped_genome=$(( Unmapped_reads_bowtie + Unmapped_reads_hisat ))
    Raw_data=$(bc <<< "scale=1; $Uniquely_mapped_genome * 100 / $Total_reads")
    Reads_ok=$(bc <<< "scale=1; $Uniquely_mapped_genome * 100 / $Passing_filters")

    echo -e "$protocole\t$date\t$library\t$seq\t$sample\t$Total_reads\t$Reads_with_adapters\t$Percent\t$Too_short\t$Percent_too_short\t$Too_long\t$Percent_too_long\t$Passing_filters\t$Percent_passing_filters%\t$rRNA\t$Contamination%\t$Uniquely_mapped_genome\t$Raw_data%\t$Reads_ok%\t$Multi_mapped\t$Unmapped_genome" >> "${output}"
done