#!/usr/bin/env python

import os
import sys
import errno
import argparse

def parse_args(args=None):
    Description = "Reformat input samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("UPDATED_PATH", help="Input path for test_nobc_nodx_rnamod")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)

def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception

def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)

def print_error(error, context:str="Line", context_str:str=""):
    """Print an error based on context and end script"""
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)

def check_samplesheet(file_in:str, updated_path:SyntaxError, file_out:str):
    """
    This function checks that the samplesheet follows the following structure:
    sample	replicate	fastq1	fastq2	strandedness	barcode	adapter	fasta	gff
    Mutant	1	Mutant.1.fastq.gz		forward		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	genome.fa	genome.gff
    Mutant	2	Mutant.1.fastq.gz		forward		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	genome.fa	genome.gff
    WT	1	WT.1.fastq.gz		forward		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	genome.fa	genome.gff
    WT	2	WT.2.fastq.gz		forward		AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	genome.fa	genome.gff
    """
    input_extensions = []
    sample_info_dict = {}
    with open(file_in, "r") as f:
        ## Check header
        min_cols = 3
        header = ["sample", "replicate", "fastq1", "fastq2", "strandedness", "barcode", "adapter",
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
            group, replicate, barcode, input_file, fasta, gtf = lspl[: len(header)]
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

            ## Check barcode entry
            if barcode:
                if not barcode.isdigit():
                    print_error("Barcode entry is not an integer!", "Line", line)
                else:
                    barcode = "barcode%s" % (barcode.zfill(2))

            ## Check input file extension
            nanopolish_fast5 = ""
            if input_file:
                if input_file.find(" ") != -1:
                    print_error("Input file contains spaces!", "Line", line)
                if input_file.endswith(".fastq.gz"):
                    input_extensions.append("*.fastq.gz")
                elif input_file.endswith(".fq.gz"):
                    input_extensions.append("*.fq.gz")
                elif input_file.endswith(".bam"):
                    input_extensions.append("*.bam")
                else:
                    if updated_path != "not_changed":
                        input_file = "/".join([updated_path, input_file.split("/")[-1]])
                    list_dir = os.listdir(input_file)
                    nanopolish_fast5 = input_file
                    if not all(fname.endswith(".fast5") for fname in list_dir):
                        if "fast5" in list_dir and "fastq" in list_dir:
                            nanopolish_fast5 = input_file + "/fast5"
                            ## CHECK FAST5 DIRECTORY
                            if not all(fname.endswith(".fast5") for fname in os.listdir(nanopolish_fast5)):
                                print_error("fast5 directory contains non-fast5 files.")
                            ## CHECK PROVIDED BASECALLED FASTQ
                            fastq_path = input_file + "/fastq"
                            basecalled_fastq = [
                                fn for fn in os.listdir(fastq_path) if fn.endswith(".fastq.gz") or fn.endswith(".fq.gz")
                            ]
                            if len(basecalled_fastq) != 1:
                                print_error("Please input one basecalled fastq per sample.")
                            else:
                                input_file = fastq_path + "/" + basecalled_fastq[0]
                                if not basecalled_fastq[0].endswith(".fastq.gz"):
                                    if not basecalled_fastq[0].endswith(".fq.gz"):
                                        print_error('basecalled fastq input does not end with ".fastq.gz" or ".fq.gz"')
                        else:
                            print_error(
                                'path does not end with ".fastq.gz", ".fq.gz", or ".bam" and is not an existing directory with correct fast5 and/or fastq inputs.'
                            )

            ## Check genome entries
            if fasta:
                if fasta.find(" ") != -1:
                    print_error("Genome entry contains spaces!", "Line", line)
                if len(fasta.split(".")) > 1:
                    if (
                        fasta[-6:] != ".fasta"
                        and fasta[-6:] != ".fasta"
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
            # gtf = ''
            is_transcripts = "0"
            if gtf:
                if gtf.find(" ") != -1:
                    print_error("Transcriptome entry contains spaces!", "Line", line)
                print(gtf[-4:])
                if gtf[-4:] != ".gtf" and gtf[-7:] != ".gtf.gz":
                    print_error("Transcriptome entry does not have extension '.gtf' or '.gtf.gz'!", "Line", line)
                # if transcriptome[-6:] != '.fasta' and transcriptome[-3:] != '.fa' and transcriptome[-9:] != '.fasta.gz' and transcriptome[-6:] != '.fa.gz' and transcriptome[-4:] != '.gtf' and transcriptome[-7:] != '.gtf.gz':
                #    print_error("Transcriptome entry does not have extension '.fasta', '.fa', '.fasta.gz', '.fa.gz', '.gtf' or '.gtf.gz'!", 'Line', line)
                # if transcriptome[-4:] == '.gtf' or transcriptome[-7:] == '.gtf.gz':
                #    gtf = transcriptome
                #    if not genome:
                #        print_error("If genome isn't provided, transcriptome must be in fasta format for mapping!", 'Line', line)
                # else:
                #    is_transcripts = '1'
                #    genome = transcriptome

            ## Create sample mapping dictionary = {group: {replicate : [ barcode, input_file, genome, gtf, is_transcripts, nanopolish_fast5 ]}}
            sample_info = [barcode, input_file, fasta, gtf, is_transcripts, nanopolish_fast5]
            if group not in sample_info_dict:
                sample_info_dict[group] = {}
            if replicate not in sample_info_dict[group]:
                sample_info_dict[group][replicate] = sample_info
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
                ",".join(["sample", "barcode", "input_file", "fasta", "gtf", "is_transcripts", "nanopolish_fast5"])
                + "\n"
            )
            for sample in sorted(sample_info_dict.keys()):
                ## Check that replicate ids are in format 1..<NUM_REPS>
                uniq_rep_ids = set(sample_info_dict[sample].keys())
                if len(uniq_rep_ids) != max(uniq_rep_ids):
                    print_error("Replicate ids must start with 1..<num_replicates>!", "Group", sample)

                ### Write to file
                for replicate in sorted(sample_info_dict[sample].keys()):
                    sample_id = "{}_R{}".format(sample, replicate)
                    fout.write(",".join([sample_id] + sample_info_dict[sample][replicate]) + "\n")

def main(args=None):
    """Run main function"""
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.UPDATED_PATH, args.FILE_OUT)

if __name__ == "__main__":
    sys.exit(main())
