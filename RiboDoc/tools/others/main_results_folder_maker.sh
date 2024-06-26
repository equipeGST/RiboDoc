#!/bin/bash

# Create a folder containing the main interesting files for biologists from RiboDoc outputs

usage() {
    echo "Usage: $0 -p <Path to working directory> -q <Qualitative analysis> -s <Stats reports> -d <DESeq2 folder> -g <Graphs_folder> -r <riboWaltz folder> -t <report table>" 1>&2
    echo -e "\n -p\tPath to working directory\n -q\tQualitative analysis chosen in configuration file ('trip' or 'ribowaltz')\n -s\tText report containing filtering and alignment statistics\n -d\tDESeq2 folder path and name\n -g\tFolder containing graphs of periodicity\n -r\tFolder containing riboWaltz outputs\n -t\tReport table"
    exit 1
}

while getopts ":p:P:q:d:s:g:r:t:" option; do
    case "${option}" in
        p)
            p=${OPTARG}
            ;;
        P)
            P=${OPTARG}
            ;;
        q)
            q=${OPTARG}
            ;;
        d)
            d=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        g)
            g=${OPTARG}
            ;;
        r)
            r=${OPTARG}
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
if [ -z "${p}" ] || [ -z "${q}" ] || [ -z "${d}" ] || [ -z "${s}" ] || [ -z "${g}" ] || [ -z "${r}" ] || [ -z "${t}" ]; then
    usage
fi;

# Create the folder's structure
mkdir -p "${p}MAIN_RESULTS/Quantitative_analysis/"

# Copy the alignments and filtering stats
cp "${s}" "${t}" "${p}MAIN_RESULTS/"

# Copy the differential analysis results
cp -r "${d}" "${p}MAIN_RESULTS/Quantitative_analysis/"


if [ "${q}" == "ribowaltz" ]; then
    # Save the phasing figures from riboWaltz
    cp "${r}frame_psite.tiff" "${r}frame_psite_length.tiff" "${r}region_psite.tiff" "${p}MAIN_RESULTS/Qualitative_analysis_riboWaltz/"
    # Save periodicity graphs
    cp -r "${g}" "${p}MAIN_RESULTS/Qualitative_analysis_riboWaltz/"
else
    # Save periodicity graphs
    cp -r "${g}" "${p}MAIN_RESULTS/Qualitative_analysis_TRiP/"
fi
