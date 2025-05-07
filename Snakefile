import pandas as pd
from optparse import OptionParser
import gffutils
import re
import os

configfile: 
    'config/config.yaml'

# Load the configuration file into the config variable
config = config

ribodoc_version = "0.9.2"


# Sets the number of threads for multi-threading steps
multi_threads_nbr = 3
mem_mb_resources = 10000
# mem_mb_resources = (workflow.cores/3)*3000
utr_threshold = "0.25"


# Sets paths for inside the container
local_path      = "/data/Enora/RiboDoc/"
ribodoc_tools   = "/data/Enora/RiboDoc/RiboDoc/tools"

stats_path = local_path + "stats/"
logs_path = local_path + "logs/"
snakemake_log_path = local_path + ".snakemake/log/"

REFERENCE_COND  = str(config['reference_condition'])
# OTHERS_COND     = list(set([SAMPLE.split("_")[0] for SAMPLE in SAMPLES if SAMPLE.split("_")[0] != REFERENCE_COND]))
LENGTHS = list(map(str,range(int(config['readsLength_min']),int(config['readsLength_max'])+1)))

# Strings with minimum and maximum read lengths to be used in file names
frag_length_S = "." + LENGTHS[0]
frag_length_L = "." + LENGTHS[0] + "-" + LENGTHS[len(LENGTHS)-1]

paired_nums = [1, 2] if config['paired_end'] == True else [1]
folder_sortmerna = ["aligned", "non_aligned"] if config['save_rRNA'] == True else ["non_aligned"]


def extract_column(file, column_index):
    with open(file, "r") as f:
        lines = f.read().strip().split("\n")[1:]  # Skip the header
        return [line.split("\t")[column_index] for line in lines]

def get_samples(wildcards):
    checkpoint_output = checkpoints.check_samplesheet.get(**wildcards).output[0]
    return extract_column(checkpoint_output, 0)

def find_fastq_path(wildcards):
    for line in open(local_path + "RESULTS/samplesheet.tsv").read().strip().split("\n")[1:]:
        parts = line.strip().split("\t")
        if parts[0] == wildcards.sample:
            return [p for p in parts[1:3] if p != ""]
    return []

def get_adapter(wildcards):
    for line in open(local_path + "RESULTS/samplesheet.tsv").read().strip().split("\n")[1:]:
        parts = line.strip().split("\t")
        if parts[0] == wildcards.sample:
            if len(parts[4]) > 0:
                return parts[4:6]
            else:
                raise ValueError(f"Adapter sequence not found for sample {wildcards.sample}.")
    raise ValueError(f"Sample {wildcards.sample} not found in samplesheet.")

def get_fasta(wildcards):
    for line in open(local_path + "RESULTS/samplesheet.tsv").read().strip().split("\n")[1:]:
        parts = line.strip().split("\t")
        if parts[0] == wildcards.sample:
            if len(parts[6]) > 0:
                return parts[6]
            else:
                raise ValueError(f"Annotation file (FASTA format) not found for sample {wildcards.sample}.")
    raise ValueError(f"Sample {wildcards.sample} not found in samplesheet.")

def fastqc_before_trimming_outputs(wildcards):
    return expand(
        local_path + "RESULTS/fastqc/before_trimming/{sample}_R{num}_fastqc.html",
        sample=get_samples(wildcards),
        num=paired_nums
    )

def trimgalore_outputs(wildcards):
    return expand(
        local_path + "RESULTS/trimgalore/{sample}_R{num}_trimgalore.fastq.gz",
        sample=get_samples(wildcards),
        num=paired_nums
    )

def fastqc_after_trimming_outputs(wildcards):
    return expand(
        local_path + "RESULTS/fastqc/after_trimming/{sample}_R{num}_fastqc.html",
        sample=get_samples(wildcards),
        num=paired_nums
    )

def sortmerna_outputs(wildcards):
    return expand(
        local_path + "RESULTS/sortMeRNA/{folder}/{sample}_R{num}.fastq.gz",
        sample=get_samples(wildcards),
        folder=folder_sortmerna,
        num=paired_nums)

def star_index_outputs(wildcards):
    return expand(
        local_path + "RESULTS/annex_database/index/{fasta_basename}",
        fasta=get_fasta(wildcards),
        fasta_basename=os.path.splitext(os.path.basename(get_fasta(wildcards)))[0])

rule all:
    input:
        local_path + "RESULTS/samplesheet.tsv",
        fastqc_before_trimming_outputs,
        trimgalore_outputs,
        fastqc_after_trimming_outputs,
        sortmerna_outputs,
        star_index_outputs
        # local_path + "RESULTS/annex_database/index/" + "{fasta_basename}", 

checkpoint check_samplesheet:
    input:
        local_path + "config/samplesheet.tsv"
    output:
        local_path + "RESULTS/samplesheet.tsv"
    shell:
        "python3 {ribodoc_tools}/others/check_samplesheet.py {input} {output}"

rule fastqc_before_trimming:
    input:
        fastq_files = lambda wildcards: find_fastq_path(wildcards)
    output:
        r1_zip  = local_path + "RESULTS/fastqc/before_trimming/{sample}_R1_fastqc.zip",
        r1_html = local_path + "RESULTS/fastqc/before_trimming/{sample}_R1_fastqc.html",
        r2_zip  = local_path + "RESULTS/fastqc/before_trimming/{sample}_R2_fastqc.zip" if config['paired_end'] == True else [], 
        r2_html = local_path + "RESULTS/fastqc/before_trimming/{sample}_R2_fastqc.html" if config['paired_end'] == True else []
    log:
        logs_path + "fastqc_before_trimming/{sample}.log"
    benchmark:
        local_path + "benchmarks/fastqc_before_trimming/{sample}.benchmark.txt"
    params:
        outdir = local_path + "RESULTS/fastqc/before_trimming/tmp_fastqc_{sample}"
    shell:
        """
        mkdir -p {params.outdir}

        files=({input.fastq_files})

        r1=${{files[0]}}
        r2=${{files[1]:-""}}

        fastqc $r1 --outdir {params.outdir} >> {log} 2>&1
        mv {params.outdir}/$(basename $r1 .fastq.gz)_fastqc.zip {output.r1_zip}
        mv {params.outdir}/$(basename $r1 .fastq.gz)_fastqc.html {output.r1_html}

        if [ -n "$r2" ]; then
            fastqc $r2 --outdir {params.outdir} >> {log} 2>&1
            mv {params.outdir}/$(basename $r2 .fastq.gz)_fastqc.zip {output.r2_zip}
            mv {params.outdir}/$(basename $r2 .fastq.gz)_fastqc.html {output.r2_html}
        fi

        rm -r {params.outdir}
        """

# Removes/cuts potential adapters on the reads
rule trimgalore:
    input:
        fastq_files = lambda wildcards: find_fastq_path(wildcards),
    output:
        r1 = local_path + "RESULTS/trimgalore/{sample}_R1_trimgalore.fastq.gz",
        r2 = local_path + "RESULTS/trimgalore/{sample}_R2_trimgalore.fastq.gz" if config['paired_end'] == True else []
    log:
        trimming = logs_path + "trimming/{sample}_trimgalore.log",
    benchmark:
        local_path + "benchmarks/trimming/{sample}.benchmark.txt"
    resources:
        mem_mb = mem_mb_resources
    threads:
        multi_threads_nbr
    params:
        min_len             = config['readsLength_min'],
        max_len             = config['readsLength_max'],
        skip_trimming       = config['skip_trimming'],
        outdir              = local_path + "RESULTS/trimgalore/",
        adapter_sequence    = lambda wildcards: get_adapter(wildcards),

    shell:  
        """
        adapter_sequence=({params.adapter_sequence});
        a1="${{adapter_sequence[0]}}";
        a2="${{adapter_sequence[1]:-""}}";

        min_len={params.min_len};
        max_len={params.max_len};
        
        files=({input.fastq_files});
        r1=${{files[0]}};
        r2=${{files[1]:-""}};
        
        if [[ "{wildcards.sample}" == *"riboseq"* ]]; then
            len_args="--length ${{min_len}} --max_length ${{max_len}}";
        else
            len_args="";
        fi 
        
        if [ {params.skip_trimming} = "True" ]; then
            trim="";
        else
            if [ -n "$r2" ]; then
                trim="-a ${{a1}} -a2 ${{a2}}";
            else
                trim="-a ${{a1}}";
            fi
        fi 
        
        if [ -n "$r2" ]; then
            mode="--paired ${{r1}} ${{r2}}";
        else
            mode="${{r1}}";
        fi;
        
        trim_galore ${{trim}} -e 0.125 -j {threads} --max_n 1 ${{len_args}} -o {params.outdir} ${{mode}} --no_report_file &> {log.trimming};
        if [ -n "$r2" ]; then
            mv {params.outdir}$(basename $r1 .fastq.gz)_val_1.fq.gz {output.r1};
            mv {params.outdir}$(basename $r2 .fastq.gz)_val_2.fq.gz {output.r2};
        else
            mv {params.outdir}$(basename $r1 .fastq.gz)_trimmed.fq.gz {output.r1};
        fi
        """

rule fastqc_after_trimming:
    input:
        r1 = local_path + "RESULTS/trimgalore/{sample}_R1_trimgalore.fastq.gz",
        r2 = local_path + "RESULTS/trimgalore/{sample}_R2_trimgalore.fastq.gz" if config['paired_end'] == True else []
        # fastq_files = lambda wildcards: [local_path + f"RESULTS/{config['trimmer']}/{wildcards.sample}_trimmed.fastq.gz"]
    output:
        r1_zip  = local_path + "RESULTS/fastqc/after_trimming/{sample}_R1_fastqc.zip",
        r1_html = local_path + "RESULTS/fastqc/after_trimming/{sample}_R1_fastqc.html",
        r2_zip  = local_path + "RESULTS/fastqc/after_trimming/{sample}_R2_fastqc.zip" if config['paired_end'] == True else [], 
        r2_html = local_path + "RESULTS/fastqc/after_trimming/{sample}_R2_fastqc.html" if config['paired_end'] == True else []
    log:
        logs_path + "fastqc_after_trimming/{sample}.log"
    benchmark:
        local_path + "benchmarks/fastqc_after_trimming/{sample}.benchmark.txt"
    params:
        skip_depletion  = config['skip_depletion'],
        outdir          = local_path + "RESULTS/fastqc/after_trimming/tmp_fastqc_{sample}"
    shell:
        """
        mkdir -p {params.outdir}

        r1={input.r1}
        r2={input.r2}

        fastqc $r1 --outdir {params.outdir} >> {log} 2>&1
        mv {params.outdir}/$(basename $r1 .fastq.gz)_fastqc.zip {output.r1_zip}
        mv {params.outdir}/$(basename $r1 .fastq.gz)_fastqc.html {output.r1_html}

        if [ -n "$r2" ]; then
            fastqc $r2 --outdir {params.outdir} >> {log} 2>&1
            mv {params.outdir}/$(basename $r2 .fastq.gz)_fastqc.zip {output.r2_zip}
            mv {params.outdir}/$(basename $r2 .fastq.gz)_fastqc.html {output.r2_html}
        fi

        rm -r {params.outdir}
        """
        
rule rrna_depletion:
    input:
        r1 = local_path + "RESULTS/trimgalore/{sample}_R1_trimgalore.fastq.gz",
        r2 = local_path + "RESULTS/trimgalore/{sample}_R2_trimgalore.fastq.gz" if config['paired_end'] == True else [],
        rrna = config['rRNA']
    output:
        # local_path + "RESULTS/annex_database/idx/",
        r1_aligned      = local_path + "RESULTS/sortMeRNA/aligned/{sample}_R1.fastq.gz" if config['save_rRNA'] == True else [],
        r2_aligned      = local_path + "RESULTS/sortMeRNA/aligned/{sample}_R2.fastq.gz" if config['save_rRNA'] == True and config['paired_end'] == True else [],
        r1_non_aligned  = local_path + "RESULTS/sortMeRNA/non_aligned/{sample}_R1.fastq.gz",
        r2_non_aligned  = local_path + "RESULTS/sortMeRNA/non_aligned/{sample}_R2.fastq.gz" if config['paired_end'] == True else []
    log:
        logs_path + "sortmerna/{sample}.log"
    benchmark:
        local_path + "benchmarks/sortmerna/{sample}.benchmark.txt"
    params:
        workdir = local_path + "RESULTS",
    shell:
        """
        if [ -n "{input.r2}" ]; then
            reads="--reads {input.r1} --reads {input.r2} --out2"
        else
            reads="--reads {input.r1}"
        fi

        mkdir "{params.workdir}/annex_database/index/sortMeRNA/"

        sample_name="{wildcards.sample}"
        sortmerna \\
            --ref {input.rrna} \\
            --threads {threads} \\
            $reads \\
            --fastx \\
            --aligned {params.workdir}/sortMeRNA/aligned/${{sample_name}} \\
            --other {params.workdir}/sortMeRNA/non_aligned/${{sample_name}} \\
            --workdir {params.workdir}/ \\
            --idx-dir {params.workdir}/annex_database/index/sortMeRNA/ \\
            --kvdb {params.workdir}/sortMeRNA/${{sample_name}}/kvdb/ \\
            --readb {params.workdir}/sortMeRNA/${{sample_name}}/readb/ \\
            > {log} 2>&1

        if [ ! -d {params.workdir}/sortMeRNA/aligned ]; then
            mv {params.workdir}/sortMeRNA/aligned/${{sample_name}}.fq.gz {output.r1_aligned}
        fi

        mv {params.workdir}/sortMeRNA/non_aligned/${{sample_name}}.fq.gz {output.r1_non_aligned}

        if [ -n "{input.r2}" ]; then
            mv {params.workdir}/sortMeRNA/non_aligned/${{sample_name}}.fq.gz {output.r2_non_aligned}
            if [ ! -d {params.workdir}/sortMeRNA/aligned ]; then
                mv {params.workdir}/sortMeRNA/aligned/${{sample_name}}.fq.gz {output.r2_aligned}
            fi
        fi

        rm -rf {params.workdir}/sortMeRNA/${{sample_name}}/
    """

rule STAR_index:
    input:
        fasta = lambda wildcards: get_fasta(wildcards)
    output:
        directory(local_path + "RESULTS/annex_database/index/" + "{fasta_basename}")
    log:
        logs_path + "STAR/{fasta_basename}_index.log"
    benchmark:
        local_path + "benchmarks/STAR/{fasta_basename}_index.benchmark.txt"
    params:
        star_dir = local_path + "RESULTS/annex_database/index/STAR",
        fasta_basename=lambda wildcards: os.path.splitext(os.path.basename(get_fasta(wildcards)))[0]
    threads:
        multi_threads_nbr
    shell:
        """
        if [[ $file == *.gff3 ]]
        then
            annotation="--sjdbGTFtagExonParentTranscript Parent"
        else
            annotation="--sjdbGTFtagExonParentTranscript transcript"
        fi

        filename=$(basename {input.fasta})
        index_folder="${filename%.*}"
        if [ ! -d {params.star_dir}/$index_folder ];
        then
            mkdir -p {params.star_dir}/$index_folder
            
            STAR \\
            --runMode genomeGenerate \\
            --genomeDir {params.star_dir}/ \\
            --genomeFastaFiles {input.fasta} \\
            # --sjdbGTFfile {input.gtf} \\
            ${{annotation}} \\
            --runThreadN {threads} \\
            --outFileNamePrefix {output}
        fi
        """

# rule STAR_ALIGN:
#     input:
#         local_path + "RESULTS/trimgalore/{sample}.trimmed" + frag_length_L + ".fastq.gz"
#     output:
#         ...
#     log:
#         logs_path + "STAR/{sample}.log"
#     benchmark:
#         local_path + "benchmarks/STAR/{sample}.benchmark.txt"
#     shell:
#         """
#         STAR \\
#             --genomeDir $index \\
#             --readFilesIn {input} \\
#             --runThreadN {threads} \\
#             --outFileNamePrefix {output} \\
#             --outSAMtype 
#         """