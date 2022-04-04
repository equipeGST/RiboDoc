# Step by step tutorial

Welcome to the RiboDoc tool tutorial !  

>RiboDoc is designed to perform all classical steps of **ribosome profiling** (RiboSeq) data analysis from the FastQ files to the differential expression analysis with necessary quality controls.


## 1) Install Docker or Singularity  
First of all, Docker or Singularity must be installed. Singularity might be prefered as it does not need super user rights, especially interesting for the use of RiboDoc on a cluster.  
Docker Engine is available on different OS like macOS and Windows 10 through Docker Desktop and as a static binary installation for a variety of Linux platforms. All are available here : https://docs.docker.com/engine/install/   

>Tips:  
>&emsp;&emsp;&emsp;For Windows, WSL2 and Ubuntu from Microsoft store applications are needed too.  

## 2) Directory preparation  
RiboDoc does not need installation but a precise architecture in your project folder is required.  
The first step is the project folder creation. It is named as your project and will be the volume linked to the container.  
Then, two sub-folders and a file have to be created and filled.   

> **Caution, those steps are majors for the good course of the analysis.**  
> **The subfolders names do not have uppercase letters.**    

Subfolders :
### a) *fastq*  
This subfolder, as its name suggests, should contain your FastQ files compressed in *.gz*.  
**Format of file names must be as following:**  
&emsp;&emsp;&emsp;***biological_condition_name.replicat_number.fastq.gz***     
For example, for the first replicate of the wild-type condition the sample can be named *Wild_Type.1.fastq.gz* and the name of the second replicate for the mutant samples would be *Mutant.2.fastq.gz*.   

>Caution, for **Windows**, extensions can be hidden.    

Folder architecture at this step:  
Project_name  
└── fastq   
&emsp;&emsp;&emsp;├── Wild_Type.1.fastq.gz   
&emsp;&emsp;&emsp;├── Wild_Type.2.fastq.gz   
&emsp;&emsp;&emsp;├── Mutant.1.fastq.gz   
&emsp;&emsp;&emsp;└── Mutant.2.fastq.gz   

If you want to try RiboDoc on an example dataset, you can find our data with on GEO : [GEO GSE173856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173856)

### b) *database*  
In this subfolder, you must put at least the following three files:  
- Your reference genome fasta file: Whether it's the genome or the transcriptome, it must be your reference fasta file to wich the reads will be aligned. We advise you to download it from the [Ensembl](https://www.ensembl.org/index.html) database, as the files are maintained up to date and following the standard GFF format.  
For example, for an entire human genome, you can look for [Human genome](https://www.ensembl.org/Homo_sapiens/Info/Index) and download the *Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz* file in the genome assembly list of fasta files, unzip it and place it in the *database* subfolder of your project directory.  

> Note:  
>&emsp;&emsp; If the transcriptome is used, the computer needs less RAM capacity than for a genome but specific files and formats can be needed.    


- A GFF format annotation file corresponding to the reference genome.  
For example, for the human genome annotations, you can look for [Human genome](https://www.ensembl.org/Homo_sapiens/Info/Index) and download the *Homo_sapiens.GRCh38.103.gff3.gz* file in the genome annotation list of gff files, unzip it and place it in the *database* subfolder of your project directory.  

- out-RNA fasta file: This file must gather together DNA sequences you want to remove from the analysis. As a rule, these are at least ribosomal sequences (rRNA). You can also add mitochondrial DNA, non-chromosomal DNA or any other fasta sequence of your choice which you want to be discarded in the analysis. You can find the fasta files associated to non-chromosomal DNA in the [Ensembl](https://www.ensembl.org/index.html) database in the genome assembly list of fasta files. If you want to remove some specific sequences, you just have to create a file with, for each sequence, one line starting with a ">" where you can add a name for your sequence followed by one line containing the sequence to remove. It gives you this format :  
>&emsp; > Sequence_X
>&emsp; GCTGACACGCTGTCCTCTGGCGACCTGTCGTCGGAGAGGTTGGGCCTCCGGATGCGCGCGGGGCTCTGGC
>&emsp; CTCACGGTGACCGGCTAGCCGGCCGCGCTCCTGCCTTGAGCCGCCTGCCGCGGCCCGCGGGCCTGCTGTT  
>&emsp; > Sequence_Y
>&emsp; CTCTCGCGCGTCCGAGCGTCCCGACTCCCGGTGCCGGCCCGGGTCCGGGTCTCTGACCCACCCGGGGGCG  

If you look for specific sequences, you can find them on the [NCBI website](https://www.ncbi.nlm.nih.gov/). For example, you can find the rRNA fasta file sequences [here](https://www.ncbi.nlm.nih.gov/nuccore/U13369.1?report=fasta).  
***You _need_ to have a file, even if it's empty***. So, if you want everything to be used in the analysis, just put an empty file.

If you need them, files containing each annotation length (see next paragraph) are also to be dropped into this folder.  

Folder architecture at this step:  
Project_name  
├── fastq   
│&emsp;&emsp;├── Wild_Type.1.fastq.gz   
│&emsp;&emsp;├── Wild_Type.2.fastq.gz   
│&emsp;&emsp;├── Mutant.1.fastq.gz   
│&emsp;&emsp;└── Mutant.2.fastq.gz  
└── database   
&emsp;&emsp;&emsp;├── reference_genome_sequences.fa  
&emsp;&emsp;&emsp;├── reference_genome_annotations.gff3  
&emsp;&emsp;&emsp;├── RNA_to_remove.fa  
&emsp;&emsp;&emsp;└── annotation_length.txt (if needed and provided)  

### c) [config.yaml](https://raw.githubusercontent.com/equipeGST/RiboDoc/main/config.yaml) file  
The *config.yaml* file allows you to define some parameters to tell RiboDoc which data you want to process and how.  
You must download it [here](https://raw.githubusercontent.com/equipeGST/RiboDoc/main/config.yaml) and open it with a text editor as a text file.    
It must be carefully completed and be present in the project directory everytime you want to run RiboDoc. A copy of this file will be made in the *RESULTS/* folder to keep a trace of the parameters you chose for a specific analysis.   

>Caution  
>&emsp;&emsp;&emsp;Spaces and quotation marks **must not be changed** ! Your information must be entered between quotes and should not have spaces     

####How to fill the configuration file :
##### Project name  
First and easy step, the project name ! You can use the same as your folder.  
*project_name*: "Project_title"  
##### Name of database files  
You must enter the full name **with extensions** without the path of files added in the database subfolder previously created.   
*fasta*: "reference_genome_fasta_file.fa"  
*gff*: "corresponding_GFF_annotation_file.gff3"  
*fasta_outRNA*: "unwanted_DNA_sequences_fasta_file.fa"  
##### Pipeline option selection  
During the RiboDoc process, data is trimmed and selected depending on their length.   
*already_trimmed*: If your data contains reads already trimmed of their adapter, you can set this option on “yes”. Else, set it on "no".   
*adapt_sequence*: If they are not trimmed, you should specify the sequence of the adapter in quotes on the line here like "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA". If you do not put anything between the quotes, RiboDoc will try to fond the adapter itself but this can sometimes lead to a wrong adapter sequence.

You also have to define the range for read length selection. Default values select reads from 25 to 35 nucleotides long.  
*readsLength_min*: minimum read length.   
*readsLength_max*: maximum read length.   

You might also need to specify features keywords in the GFF file to fit your GFF file format :   
*gff_element_cds*: feature corresponding to CDS in the annotation file. "CDS" is the default value (can sometimes be "ORF").     
*gff_attribut*: attribut to regroup reads during counting. "Parent" is the default value to regroup counts by transcripts (Parent of CDS features) for the differential analysis. You can also put "ID" if you want to count each CDS independently.
*gff_name_attribut*: Name of the genes features in the GFF. Default is "Name" but it can sometimes be "gene_name" or else.     
##### Statistical settings  
To be able to perform statistical analyzes, you must define a reference condition as well as your thresholds.   
*reference_condition*: it correspond to the reference biological_condition_name in your FastQ files name. Ex : "Wild_Type" (as in *Wild_Type.1.fastq.gz*).
*transcript_or_gene*: choose whether to perfom the differential analysis on features grouped by "transcripts" or by "gene"  
*p-val*: p-value threshold for the differential analysis. Default value is 0.01.  
*logFC*: logFC threshold for the differential analysis. Defaut value is 0 to keep all the genes without logFC filtering.  
##### Window for qualitative test  
During the quality analysis, the periodicity is observed on nucleotides around start and stop codons.     

>The periodicity must be calculated using a metagene profile (metaprofile). It provides the amount of footprints relative to all annotated start and stop codons.   

2 pipelines dedicated to quality controls are available in RiboDoc. The first one uses the [riboWaltz tool](https://github.com/LabTranslationalArchitectomics/riboWaltz) which can need high RAM resources depending on your data. The second pipeline is a series of scripts called TRiP which use a specific gff file format as annotation file. For more details on this format, please check the RiboDoc article.
*qualitative_analysis*: Choose between the 2 qualitative analysis pipeline. Default value is "ribowaltz" as it is more precise and does not need files with specific formats. "trip" is more complicated to use but asks for less resources.  
The window selected by default is -50/+100 nts and -100/+50 nts around start and stop codons respectively.   
*window_utr*: Define your window before start and after stop. Default is "50"  
*window_cds*: Define your window after start and before stop. Default is "100"  

################################################################
###### Optional and only for qualitative_analysis: "trip" ######
################################################################
##### UTR covering option  
*UTR*: This option has to be turned on if you want to compare UTRs coverage against CDS. However, to be realized, this option requires a file with the name of each gene and the length of the associated annotation (one file by annotation: CDS, 3’-UTR and 5’-UTR).  

The full name (**with extensions**) without the path of each file has to be report in the configfile.  
*CDSlength*: complete name of file with CDS length.  
*5primelength*: complete name of file with 5’-UTR length  
*3primelength*: complete name of file with 3’-UTR length  

Elements for UTRs have to be specified only if you put 'yes' for the 'UTR' option  
*gff_element_three_prime_utr*: feature corresponding to 3'UTR in the annotation file. Default value is "three_prime_UTR"  
*gff_element_five_prime_utr*: feature corresponding to 5'UTR in the annotation file. Default value is "five_prime_UTR"  
################################################################  


Folder architecture at this step:  
Project_name  
├── fastq   
│&emsp;&emsp;├── Wild_Type.1.fastq.gz   
│&emsp;&emsp;├── Wild_Type.2.fastq.gz   
│&emsp;&emsp;├── Mutant.1.fastq.gz   
│&emsp;&emsp;└── Mutant.2.fastq.gz  
├── database   
│&emsp;&emsp;├── reference_genome_sequences.fa  
│&emsp;&emsp;├── reference_genome_annotations.gff3  
│&emsp;&emsp;├── RNA_to_remove.fa  
│&emsp;&emsp;└── annotation_length.txt (if needed and provided)  
└── config.yaml  


## 3) Pull RiboDoc  
To get RiboDoc, open a terminal. We know this can look scary for some but do not worry, there are only 2 easy lines to write.  
If you never used RiboDoc on your workstation, you must pull it from Docker hub.  
Copy and paste the following command line:    
&emsp;&emsp;&emsp;`docker pull equipegst/ribodoc`  
If you have any error, it might come from a rights problem so you should try to copy and paste this command:    
&emsp;&emsp;&emsp;`sudo docker pull equipegst/ribodoc`     
### Use of singularity    
If you want to use a singularity image instead of a docker image (for the use on a cluster for example), you must pull it from dockerhub. You can name the image as you want but the extension should be *.sif* :    
&emsp;&emsp;&emsp;`singularity pull /path/to/your/singularity/image.sif docker://equipegst/ribodoc`    


## 4) Run RiboDoc  
Now that your folder's architecture is ready and the image is on your computer, it's time to start ! If you have pulled RiboDoc, you can run it with the following command:  
### Use of docker  
&emsp;&emsp;&emsp;`docker run --rm -v /path/to/your/project/folder/:/data/ equipegst/ribodoc CPU_NUMBER MEMORY_AMOUNT`  
### Use of singularity  
&emsp;&emsp;&emsp;`singularity run --disable-cache --bind /path/to/your/project/folder/:/data/ equipegst/ribodoc CPU_NUMBER MEMORY_AMOUNT`  
*/path/to/your/project/folder/* corresponds to the **full path** of the directory with all the files you prepared for the analysis (usually starting with a "/"). To get this full path, you can usually drag and drop your folder to the terminal or find it in the properties.    
*:/data/* **must not be modified in any way** as it corresponds to the path to the project directory inside the container.   
*CPU_NUMBER* is an integer corresponding to the total number of threads you want to use for your analysis (Ex : 4).   
*MEMORY_AMOUNT* is an integer corresponding to the maximum amount of memory in Gigabytes (RAM) you want to use for your analysis (Ex : 30). Be careful to have enough memory available for alignment steps or the pipeline will be stopped without a specific error written in the "logs" folder.   


>Caution:   
>&emsp;&emsp;&emsp;The path to your project folder and the "/data/" path must start and finish with a slash "/" if you are using docker.   


>For Windows users  
>&emsp;&emsp;&emsp;The path to your project folder has to start at the local disks C or D: C:\your\path\   
>&emsp;&emsp;&emsp;This path to your project folder has to be composed and finished with backslashes "\" (instead of slashes "/")   
>&emsp;&emsp;&emsp;/data/ path does not change in any way !   

## 5) In case of any error   
Managed by snakemake, the pipeline will finish all jobs unrelated to the rule/job that failed before exiting. You can still force the container to stop with Docker Desktop or with the following command lines (might need the "sudo" keyword at the beginning) :   
>&emsp;> docker container ls   

Which provides you the container's ID (*e.g.* 9989909f047d), then :   
>&emsp;> docker stop ID  
Where "ID" is the id obtained with the previous command

If you have issues with the use of Docker, you must refer to their [website](https://docs.docker.com/).     
If the error happens during the use of RiboDoc, the rule (job) which failed is indicated in your terminal. You can then find the error output in the *logs* folder. Each rule have a precise name and a folder related to it with files corresponding to the different steps of this rule. If you can not solve the problem by exploring those files, you can contact us through the ["issues"](https://github.com/equipeGST/RiboDoc/issues) part of our github.
>If you want to send us a request because of an error, the easiest way for us to help you is if you send us an archive named after the rule which failed in your analysis with your config.yaml, logs folder and logsTmp folder to it so we can help you finding what happened and why.

## 6) Understand the results  
Here is the project_name folder architecture after RiboDoc run.  
Initial folders and files are still present and highligth in bold in the tree architecture below.  
**Project_name**  
├── **fastq/**   
│&emsp;&emsp;├── **Wild_Type.1.fastq.gz**   
│&emsp;&emsp;├── **Wild_Type.2.fastq.gz**   
│&emsp;&emsp;├── **Mutant.1.fastq.gz**   
│&emsp;&emsp;└── **Mutant.2.fastq.gz**  
├── database/   
│&emsp;&emsp;├── **reference_genome_sequences.fa**  
│&emsp;&emsp;├── **reference_genome_annotations.gff3**  
│&emsp;&emsp;├── **RNA_to_remove.fa**  
│&emsp;&emsp;└── **annotation_length.txt (if provided)**   
├── **config.yaml**   
├── logs/   
│&emsp;&emsp;└── *one_log_folder_by_job*  
├── logsTmp/   
│&emsp;&emsp;└── *one_file_by_steps_of_interest_for_alignment_stats*  
├── RESULTS/  
│&emsp;&emsp;├── config.yaml     
│&emsp;&emsp;├── BAM/  
│&emsp;&emsp;│&emsp;&emsp;├── *one_bam_by_sample.bam*  
│&emsp;&emsp;│&emsp;&emsp;└── *one_bai_by_bam.bai*  
│&emsp;&emsp;├── DESeq2_by_gene or DESeq2_by_transcript/  
│&emsp;&emsp;│&emsp;&emsp;├── count_matrix.txt  
│&emsp;&emsp;│&emsp;&emsp;├── count_matrix_by_gene.txt  
│&emsp;&emsp;│&emsp;&emsp;├── complete.txt  
│&emsp;&emsp;│&emsp;&emsp;├── up.txt  
│&emsp;&emsp;│&emsp;&emsp;├── down.txt  
│&emsp;&emsp;│&emsp;&emsp;├── project_name.Final_report.html  
│&emsp;&emsp;│&emsp;&emsp;├── Images/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *all_quantitative_analysis_graphs*  
│&emsp;&emsp;├── fastqc/  
│&emsp;&emsp;│&emsp;&emsp;├── *one_html_by_sample.html*  
│&emsp;&emsp;│&emsp;&emsp;└── *one_zip_by_sample.zip*  
│&emsp;&emsp;├── htseqcount_CDS/  
│&emsp;&emsp;│&emsp;&emsp;└── *one_file_by_sample.txt*  
│&emsp;&emsp;├── riboWaltz/  # If you chose *qualitative_analysis: "ribowaltz"* in the config.yaml file  
│&emsp;&emsp;│&emsp;&emsp;├── *riboWaltz's qualitative analysis results*  
│&emsp;&emsp;├── qualitativeAnalysis/  # If you chose *qualitative_analysis: "trip"* in the config.yaml file  
│&emsp;&emsp;│&emsp;&emsp;├── *TRiP's qualitative analysis results*  
├── dag_all.svg    
└── dag_last_run.svg    
- The *logs/* folder groups together all the error output messages from tools used in RiboDoc analysis pipeline. Thus, in the event of an error, it allows you to identify the problematic step to give us feedback.   
- *RESULTS/* folder contains 5 subfolders:   
&emsp;&emsp;i) *BAM/*: it contains a BAM file for each sample (allows visualization on tools such as IGV).  
&emsp;&emsp;ii) *DESeq2/*: it contains the differential analysis html report (*PROJECT_NAME.Final_report.txt*), the count_matrix, the tables and the images related to the DESeq2 analysis.   
&emsp;&emsp;iii) *fastqc/*: it contains raw data quality controls.   
&emsp;&emsp;iv) *htseqcount_CDS/*: it contains htseq output for CDS counts.    
&emsp;&emsp;vii) *qualitativeAnalysis/* or *riboWatlz/*: it contains all files related to qualitative test like periodicity and reads length repartition   
It contains also two files:  
&emsp;&emsp;i) *PROJECT_NAME.Analysis_report.html* gathers standard output of each analysis pipeline tool. It allows to know numbers of reads at each step a)raw reads b)reads after trimming and length selection c)after out RNA depletion d)after double alignment on the reference genome.  
&emsp;&emsp;iii) *config.yaml* to have a parameters backup.     

- The *dag files* which represents the analysis steps with your samples.  

>Last big tip:  
In case a sample is too variable against other replicats or if new sequenced samples are to be added to your study, you can delete/move or add them in the *fastq* subfolder. RiboDoc will only process necessary steps based on the fastq files list. If you want to keep the previous results of your differential analysis, delete/move/rename the subfolder *RESULTS/DESeq2*. Run again RiboDoc on the same *project_name* folder and it only creates missing files to complete the analysis.  

Thank you for using RiboDoc !   
We wish you the best results for your analysis !  
