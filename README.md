
## Introduction

**RiboDoc** is a bioinformatics pipeline for Ribosome sequencing (Ribo-seq) data made to perform quality control, trimming, alignment and downstream qualitative and quantitative analysis.

It can be used with multiple operating systems, and it's goal is to standardize the general steps that must be performed systematically in Ribo-seq analysis, together with the statistical analysis and quality control of the sample. The data generated can then be exploited with more specific tools.

RiboDoc is a tool designed to standardize bioinformatics analyses in the field of translation, following the [FAIR](https://www.go-fair.org/fair-principles/) guidelines to make installation and analysis meet principles of findability, accessibility, interoperability, and reusability. Thus, this pipeline is built using [Snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system to create reproducible and scalable data analyses. Additionally, it is a Docker-based package, which means it can be used by anyone. [Docker](https://www.docker.com/) is a container which packages up code and all its dependencies so RiboDoc can run quickly and reliably from one computing environment to another. 

If you want to easily understand how to launch RiboDoc on your own computer, you can check our video tutorial just here :
[![RiboDoc_Video](https://github.com/equipeGST/RiboDoc/assets/75135539/51fdeaf2-8a7c-4cd8-a5b3-d58d3c3adf5e)](http://www.youtube.com/watch?v=e9_SFz_YEK0 "Tutorial - RiboSeq analysis with RiboDoc pipeline")

## Pipeline summary

RiboDoc is designed to perform all classical steps of **ribosome profiling** (RiboSeq) data analysis from the FastQ files to the differential expression analysis with necessary quality controls.

1. Quality Control of raw reads with [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Adapter and quality trimming, read length filtering with [`Cutadapt`](https://cutadapt.readthedocs.io/en/stable/)
3. Quality Control of trimmed reads with [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
4. Removal of contaminants RNA (rRNA, tRNA, viral RNA, ...) with [`Bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
5. Quality Control of depleted reads with [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
6. Genome and transcriptome alignment of reads conjointly with [`Hisat2`](http://daehwankimlab.github.io/hisat2/) and [`Bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
7. Sort and index alignments with [`samtools`](https://sourceforge.net/projects/samtools/files/samtools/)
8. Reads Count with [`htseq-count`](https://htseq.readthedocs.io/en/release_0.11.1/count.html#)
9.Analysis of differential gene expression with ['DESeq2'](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
10. Offset prediction and periodicity graph creation with [`ribowaltz`](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006169) or `TRiP`

### Files

- The `MAIN_RESULTS` folder contains a copy of the most interesting files present in the `RESULTS` folder.
- The `RESULTS` folder contains these subfolders:  
    - `adapter_lists`: contains text files with adapters list for each sample that were found in the `config.yaml` file or determined from data is the user did not provide any adapter sequence.  
    - `annex_database`: contains the indexes for the alignments, the re-formatted FASTA and GFF files for the analysis and the GTF file for the riboWaltz pipeline.
    - `BAM.25-35`: contains a BAM format alignment file (`*.bam`) for each sample.  
    - `BAM_transcriptome.25-35`: contains a BAM file for each sample generated from the transcriptome GTF and FASTA files generated for riboWaltz.  
    - `DESeq2_CDS.25-35`: contains the differential analysis html report (`Yeast_RiboSeq.Final_report.html`), the count matrices, the tables and the images related to the differential analysis, regrouped by gene or by transcript.  
    - `fastqc`: contains data quality controls FastQC reports of pre-processed and trimmed read.
    - `HTSeq-counts`: contains HTseq output for CDS counts (and UTR regions if selected in the `config.yaml` file).
    - `qualitativeAnalysis` or `riboWaltz.25-35` and `periodicity_-30+90`: they contain all files related to qualitative test like metaprofiles and reads lengths repartition.
It contains also three files:  
    - `Yeast_RiboSeq.Analysis_Report.25-35.txt` gathers standard output of each analysis main tool. It allows to know numbers of reads at each step of the analysis.  
    - `Yeast_RiboSeq.Analysis_Table_summary.25-35.csv` summarizes the same standard outputs as previous file in a table.  
    - `config.yaml` is a copy of `config.yaml` file at the root of your project for reproductability.  
- `logs` folder groups together all the standard errors messages from tools used in RiboDoc pipeline. Thus, in the event of an error, it allows you to identify the problematic step.  
- `stats` groups all main statistics which the `Analysis_Report` and `Analysis_Table_summary` files are made from.  
- `dag_all.svg` file is a graphical representation of all jobs where the edges represent dependencies.  

> In case a sample is too variable against other replicates or if new sequenced samples are to be added to your study, you can delete/move or add them in the `fastq` subfolder. RiboDoc will only process necessary steps based on the fastq files list.

> Containers state saves and cache can take more and more space if you do not use the `--rm` option for docker. To clean everything related to containers, just right in your terminal : `docker system prune -af`

## Citations

If you use RiboDoc for your analysis, please cite it using the following doi: [10.1016/j.csbj.2021.05.014](https://doi.org/10.1016/j.csbj.2021.05.014).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
