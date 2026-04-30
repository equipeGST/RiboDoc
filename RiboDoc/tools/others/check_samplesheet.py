#!/usr/bin/env python

import os
import sys
import errno
import argparse
from pathlib import Path

def parse_args(args=None) -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Args:
        args (list) : List of arguments to parse.
        
    Returns:
        argparse.Namespace : Object containing parsed arguments.
    """
    Description = "Reformat input samplesheet file and check its contents."
    Epilog      = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "FILE_IN", 
        help = "Input samplesheet file.")
    parser.add_argument(
        "FILE_OUT",
        help = "Output file.")
    return parser.parse_args(args)

def make_dir(path):
    """
    Make a directory if it doesn't exist.
    """
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception

def read_head(handle, num_lines:int=10) -> str:
    """
    Read the specified number of lines from the current position in the file.
    
    Args:
        handle (str)   : File handle.
        num_lines (int) : Number of lines to read.
        
    Returns:
        str : Lines read from the file.
    """
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)

def print_error(error, context:str="Line", context_str:str=""):
    """
    Print an error based on context and end script
    
    Args:
        error (str)       : Error message.
        context (str)     : Context of error.
    """
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)

def check_samplesheet(file_in:str, file_out:str) -> None:
    """
    This function checks that the samplesheet follows the following structure:
    
    condition	replicate	type	fastq1	fastq2	strandedness	barcode1	barcode2	adapter1	adapter2	fasta	gff
    Mutant	1	riboseq	Mutant.1.fastq.gz		forward		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	genome.fa	genome.gff
    Args:
        file_in (str)  :        Input samplesheet file.
        file_out (str) :        Output samplesheet file.
    Returns:
        None
    """
    input_extensions = []
    sample_info_dict = {}
    with open(file_in, "r") as f:
        ## Check header
        min_cols = 3
        header = ["condition", "replicate", "type", "fastq1", "fastq2", "strandedness", "barcode1", "barcode2", "adapter1", "adapter2",
                  "fasta", "gff"]
        line = f.readline().strip().split("\t")
        if line[: len(header)] != header:
            print("ERROR: Please check samplesheet header -> {} != {}".format(
                ",".join(line), ",".join(header)))
            sys.exit(1)

        ## Check sample entries
        for line in f:
            lspl = [x.strip() for x in line.strip().split("\t")]

            ## Check valid number of columns per row
            if len(lspl) < len(header):
                print_error("Invalid number of columns (minimum = {})!".format(len(header)), "Line", line)

            num_cols = len([x for x in lspl if x])
            if num_cols < min_cols:
                print_error("Invalid number of populated columns (minimum = {})!".format(min_cols), "Line", line)

            ## Check group name entries
            group, replicate, seq_type, fastq1, fastq2, strandedness, barcode1, barcode2, adapter1, adapter2, fasta, gtf = lspl[: len(header)]
            if group:
                if group.find(" ") != -1:
                    print_error("Group entry contains spaces!", "Line", line)
            else:
                print_error("Group entry has not been specified!", "Line", line)

            ## Check replicate entry is integer
            if replicate:
                if not replicate.isdigit():
                    print_error("Replicate id not an integer!", "Line", line)
            else:
                print_error("Replicate id not specified!", "Line", line)
            replicate = int(replicate)

            ## Check type name entries
            if seq_type:
                if seq_type.find(" ") != -1:
                    print_error("Type entry contains spaces!", "Line", line)
                if seq_type not in ["riboseq", "rnaseq"]:
                    print_error("Type entry must be 'riboseq' or 'rnaseq'!", "Line",
                                line)
            else:
                print_error("Type entry not specified!", "Line", line)

            # Check input file(s)
            if fastq1:
                if fastq2:
                    input_file = [fastq1, fastq2]
                else:
                    input_file = [fastq1]
                
                for fastq in input_file:
                    file = Path(fastq)
                    if not file.is_file():
                        
                        print_error("One of your input file does not exist!", "Line", line)
                    if fastq.find(" ") != -1:
                        print_error("One of your input FASTQ file contains spaces!", "Line", line)
                    if fastq.endswith(".fastq.gz"):
                        input_extensions.append("*.fastq.gz")
                    elif fastq.endswith(".fq.gz"):
                        input_extensions.append("*.fq.gz")
                    else:
                        print_error("Path does not end with '.fastq.gz' or '.fq.gz' extension!", "Line", line)
            else:
                print_error("Fastq1 entry not specified!", "Line", line)

            ## Check strandedness entry
            if strandedness:
                if strandedness.find(" ") != -1:
                    print_error("Strandedness entry contains spaces!", "Line", line)
                if strandedness not in ["forward", "reverse", "auto"]:
                    print_error("Strandedness entry should be either forward, reverse or auto!", "Line", line)
            else:
                print_error("Strandedness entry not specified!", "Line", line)

            ## Check barcode entry
            if barcode1:
                # if not barcode1.isdigit():
                #     print_error("Barcode entry is not an integer!", "Line", line)
                # else:
                barcode1 = "barcode%s" % (barcode1.zfill(2))

            ## Check adapter entry
            if adapter1:
                if adapter1.find(" ") != -1:
                    print_error("Adapter entry contains spaces!", "Line", line)
                if not all(x in ["A", "T", "C", "G"] for x in adapter1):
                    print_error("Adapter entry contains others characters than A, T, C or G ! ", "Line", line)
            # else:
            #     print_error("adapter entry not specified!", "Line", line)

            ## Check genome entries
            if fasta:
                file = Path(fasta)
                if not file.is_file():
                    print_error("Fasta file does not exist!", "Line", line)
                if fasta.find(" ") != -1:
                    print_error("Genome entry contains spaces!", "Line", line)
                if len(fasta.split(".")) > 1:
                    if (
                        fasta[-6:] != ".fasta"
                        and fasta[-3:] != ".fa"
                        and fasta[-4:] != ".fna"
                        and fasta[-7:] != ".fna.gz"
                        and fasta[-9:] != ".fasta.gz"
                        and fasta[-6:] != ".fa.gz"
                    ):
                        print_error(
                            "Genome entry does not have extension '.fasta', '.fa', '.fna', 'fna.gz', '.fasta.gz' or '.fa.gz'!",
                            "Line",
                            line,
                        )
                
            ## Check transcriptome entries
            is_transcripts = "0"
            if gtf:
                file = Path(gtf)
                if not file.is_file():
                    print_error("GFF file does not exist!", "Line", line)
                if gtf.find(" ") != -1:
                    print_error("Genes annotation entry contains spaces!", "Line", line)
                if (
                    gtf[-4:] != ".gtf"
                    and gtf[-7:] != ".gtf.gz"
                    and gtf[-5:] != ".gff3"
                    and gtf[-8:] != ".gff3.gz"
                    and gtf[-4:] != ".gff"
                    and gtf[-7:] != ".gff.gz"
                ):
                    print_error("Genes annotation entry does not have extension '.gtf', 'gff', 'gtf.gz' or '.gff.gz'!", "Line", line)

            ## Create sample mapping dictionary = {group: {replicate : [ barcode, input_file, genome, gtf, is_transcripts, nanopolish_fast5 ]}}
            sample_info = [fastq1, fastq2, barcode1, barcode2, adapter1, adapter2, fasta, gtf, is_transcripts]
            
            if f"{group}_{seq_type}" not in sample_info_dict:
                sample_info_dict[f"{group}_{seq_type}"] = {}
            if replicate not in sample_info_dict[f"{group}_{seq_type}"]:
                sample_info_dict[f"{group}_{seq_type}"][replicate] = sample_info
            else:
                print_error("Same replicate id provided multiple times!", "Line", line)

    ## Check all input files have the same extension
    if len(set(input_extensions)) > 1:
        print_error(
            "All input files must have the same extension!",
            "Multiple extensions found",
            ", ".join(set(input_extensions)),
        )

    ## Write validated samplesheet with appropriate columns
    if len(sample_info_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                "\t".join(["sample", "fastq1", "fastq2", "barcode1", "barcode2", "adapter1", "adapter2", "fasta", "gtf", "is_transcripts"])
                + "\n"
            )
            for sample in sorted(sample_info_dict.keys()):
                ## Check that replicate ids are in format 1..<NUM_REPS>
                uniq_rep_ids = set(sample_info_dict[sample].keys())
                if len(uniq_rep_ids) != max(uniq_rep_ids):
                    print_error("Replicate ids must start with 1..<num_replicates>!", "Group", sample)

                ### Write to file
                for replicate in sorted(sample_info_dict[sample].keys()):
                    sample_id = "{}.{}".format(sample, replicate)
                    fout.write("\t".join([sample_id] + sample_info_dict[sample][replicate]) + "\n")

def main(args=None):
    """
    Run main function
    
    Args:
        args (list) : List of arguments to parse.
    """
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)

if __name__ == "__main__":
    sys.exit(main())
