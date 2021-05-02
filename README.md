# Step by step tutorial

Welcome to the RiboDoc tool tutorial !  

>RiboDoc is designed to perform all classical steps of **ribosome profiling** (RiboSeq) data >analysis from the FastQ files to the differential expression analysis.


## 1) Install Docker  
First of all, Docker must be present in version 19 or higher.  
If you don’t already have it, now is the time to fix it !  
Docker Engine is available on different OS like macOS and Windows 10 through Docker Desktop and as a static binary installation for a variety of Linux platforms. All are available here : https://docs.docker.com/engine/install/   

>Tips:  
>&emsp;&emsp;&emsp;For Windows, WSL2 and Ubuntu from Microsoft store applications are needed too.  

## 2) Directory preparation  
RiboDoc does not need installation (yipee) but a precise folder architecture is required (boo).  
The first step is the project folder creation. It is named as your project and will be the volume linked to Docker.  
Then, two sub-folders and a file have to be created and completed respectively.   

> **Caution, those steps are majors for the good course of the analysis**  
> **Subfolders don’t have uppercase**    

Folder architecture at this step:  
Project_name  

### a) *fastq* subfolder  
This subfolder, as its name suggests, should contain your FastQs. These must be compressed in .gz.  
**Format of file name must be as following:**  
&emsp;&emsp;&emsp;biological_condition_name.replicat.fastq.gz     
For example, for the first replicat of the wild-type condition, sample will be named *WT.1.fastq.gz*   

>Caution, for **Windows**, extensions can be hidden.    

Folder architecture at this step:  
Project_name  
└── fastq   
&emsp;&emsp;&emsp;├── condA.1.fastq.gz   
&emsp;&emsp;&emsp;├── condA.2.fastq.gz   
&emsp;&emsp;&emsp;├── condB.1.fastq.gz   
&emsp;&emsp;&emsp;└── condB.2.fastq.gz   


### b) *database* subfolder  
In this subfolder, you must put at least the following three files:  
- your genome fasta file: Whether it's the genome or the transcriptome, it must be your reference fasta file where reads will be aligned. It must, like the other files, be downloaded from the [Ensembl](https://www.ensembl.org/index.html) database.  

> Note:    
>&emsp;&emsp;For complexe genomes, transcriptome is preferable (less complexity and better analysis).   
>&emsp;&emsp;If the genome is used, the computer needs a larger RAM capacity than for transcriptome.    


- GFF3 file corresponding to the reference genome dropped.  
- out-RNA fasta: This fasta file must gather together RNA sequences you want to remove. As a rule, these are at least rRNA. You can also add mitochondrial RNA or any other RNA. You need to have a file, even if it's empty. So, if you want everything to be used in the analysis, just put an empty file.

If you have them, files containing each annotation length (see next paragraph) also be dropped into this folder.  

Folder architecture at this step:  
Project_name  
├── fastq   
│&emsp;&emsp;├── condA.1.fastq.gz   
│&emsp;&emsp;├── condA.2.fastq.gz   
│&emsp;&emsp;├── condB.1.fastq.gz   
│&emsp;&emsp;└── condB.2.fastq.gz  
└── database   
&emsp;&emsp;&emsp;├── reference_transcriptome.fa  
&emsp;&emsp;&emsp;├── reference_transcriptome.gff3  
&emsp;&emsp;&emsp;├── RNA_to_remove.fa  
&emsp;&emsp;&emsp;└── annotation_length.txt (if possible)  

### c) [config.yaml](https://raw.githubusercontent.com/equipeGST/RiboDoc/main/config.yaml) file  
Config.yaml file is used to define parameters to tell RiboDoc how to process your data.  
You must download it [here](https://raw.githubusercontent.com/equipeGST/RiboDoc/main/config.yaml) and open it with a text editor.    
It must be carefully completed and be present in the project directory everytime you want to run RiboDoc.  

>Caution  
>&emsp;&emsp;&emsp;Spaces and quotation marks **must not be changed**! Your information must be entered in quotes     

#### Project name  
First and easy step, the project name ! You could (and it is recommended) use the same that your folder.  
*project_name*: principal directory name  
#### Name of database subfolder file  
You must enter the full name (**with extensions**) without the path of files added in the database subfolder previously created.   
*fasta*: reference transcriptome (or genome) complete name.  
*gff*: corresponding GFF3 name.  
*fasta_outRNA*: not interesting sequence fasta name.  
#### UTR covering option  
*UTR*: This option has to be turned on if you want to compare UTRs coverage against CDS. However, to be realized, this option requires a file with the name of each gene and the length of the associated annotation (one file by annotation: CDS, 3’-UTR and 5’-UTR).  
The full name (**with extensions**) without the path of each file has to be report in the configfile.  
*CDSlength*: complete name of file with CDS length.  
*5primelength*: complete name of file with 5’-UTR length  
*3primelength*: complete name of file with 3’-UTR length  
#### Pipeline option selection  
During the RiboDoc process, data is trimmed and selected depending on their length.   
*already_trimmed*: In case your data was already trimmed, you can set this option on “yes”.   
*adapt_sequence*: Else, you must specify the sequence adapters in quotes on the line below.  
You also have to define the range for read length selection. Default values are already filled.  
*readsLength_min*: minimum read length.   
*readsLength_max*: maximum read length.   
You can also parameters alignments:   
*times_mapped*: you can put 0 or 1. 0 for uniquely mapped reads and 1 for multi-mapped reads.    
#### Attribut setting     
*gff_attribut*: attribut to regroup reads during counting. 'ID' is default value.     
#### Statistical settings  
To be able to perform statistical analyzes, you must define a reference condition as well as your thresholds.   
*reference_condition*: it correspond to the reference biological_condition_name
We have pre-define them:  
*p-val*: defined at 0.01 to keep only genes with a high confidence.  
*logFC*: defined at 0 to keep all the genes.  
#### Window for qualitative test  
During the quality analysis, the periodicity is observed on bases around start and stop.     

>The periodicity must be calculated using a metagene profile. It provides the amount of footprints relative to all annotated start and stop codons in a selected window.   

The window selected by default is -50/+100 nts and -100/+50 nts around start and stop codons respectively.   
*window_bf*: define your window before start and after stop  
*window_af*: define your window after star and before stop  
#### Number of threads  
Thanks to the use of Snakemake, RiboDoc can analyse many samples at the same time. We define that ¼ of available CPUs are necessarily requisitioned for this multiple tasks in parallel.  
As the majority of tools used in the analysis pipeline have a multithreaded option, you can choose the number of threads you allow additionally.  
*threads*: number of threads you allowed.    

>Caution:  
>&emsp;&emsp;&emsp;You could attribute maximum 4 threads.  
>&emsp;&emsp;&emsp;Indeed:  
>&emsp;&emsp;&emsp;Allow 1 thread = Remain on a ¼ of threads used  
>&emsp;&emsp;&emsp;Allow 2 threads = Each sample will have 2 threads at its disposal. 2/4 of the threads will therefore be used.   
>&emsp;&emsp;&emsp;Allow 3 threads = ¾ of CPUs are used. We advise not to go beyond so as not to saturate the computer and to be able to continue to use it.     

Folder architecture at this step:  
Project_name  
├── fastq   
│&emsp;&emsp;├── condA.1.fastq.gz   
│&emsp;&emsp;├── condA.2.fastq.gz   
│&emsp;&emsp;├── condB.1.fastq.gz   
│&emsp;&emsp;└── condB.2.fastq.gz  
├── database   
│&emsp;&emsp;├── reference_transcriptome.fa  
│&emsp;&emsp;├── reference_transcriptome.gff3  
│&emsp;&emsp;└── RNA_to_remove.fa  
└── config.yaml  

Don't forget to save files !        

## 3) Pull RiboDoc  
When the folder architecture is ready, it’s time to take a RiboDoc !  
First, open a terminal.  
If you never had use RiboDoc on your workstation, you must pull it from Docker hub.  
Copy and past the following command line:    
&emsp;&emsp;&emsp;`docker pull equipegst/ribodoc`  
If you have rights problem, copy and past this command:    
&emsp;&emsp;&emsp;`sudo docker pull equipegst/ribodoc`     

## 4) Run RiboDoc  
Then, or if you already have used RiboDoc, you can run it thanks to the following command:  
&emsp;&emsp;&emsp;`docker run --rm -v /path/to/working/directory/:/data/ equipegst/ribodoc`  
*/path/to/working/directory/* corresponds to the project_name directory full path.   
*:/data/* **must not be modified in any way**   

>Caution:   
>&emsp;&emsp;&emsp;All paths must start and finish with a slash “/”   


>For Windows users  
>&emsp;&emsp;&emsp;Path has to start at the local disks C or D: C:\your\path\   
>&emsp;&emsp;&emsp;Path has to be composed and finished with backslashes “\”   
>&emsp;&emsp;&emsp;/data/ path do not change !   

## 5) In case of any error   
If you have issues with the use of Docker, you must refer to their [website](https://docs.docker.com/).     
If the error happens during the use of RiboDoc, the rule (job) which failed in your terminal. You can then find the error output in the *logs* folder. Each rule have a precise name and a file related to it. If you can not solve the problem, contact us either through ["issues" part of our github](https://github.com/equipeGST/RiboDoc/issues).
>Please note:
>&emsp;&emsp;&emsp;If you want to send us a request, join an archive with your config.yaml, logs folder and logsTmp folder to it so we could help you finding what happened and why.

## 6) Understand results  
Here the project_name folder architecture after RiboDoc run.  
Initial folders and files are still present and highligth in bold in the tree architecture below.  
**Project_name**  
├── **fastq/**   
│&emsp;&emsp;├── **condA.1.fastq.gz**   
│&emsp;&emsp;├── **condA.2.fastq.gz**   
│&emsp;&emsp;├── **condB.1.fastq.gz**   
│&emsp;&emsp;└── **condB.2.fastq.gz**  
├── database/   
│&emsp;&emsp;├── **reference_transcriptome.fa**  
│&emsp;&emsp;├── reference_transcriptome.index.ht2  
│&emsp;&emsp;├── reference_transcriptome.index.bt2  
│&emsp;&emsp;├── **reference_transcriptome.gff3**  
│&emsp;&emsp;├── **RNA_to_remove.fa**  
│&emsp;&emsp;├── RNA_to_remove.index.bt2  
│&emsp;&emsp;└── **annotation_length.txt (if possible)**   
├── logs/   
│&emsp;&emsp;└── *one_log_by_step*  
├── RESULTS/  
│&emsp;&emsp;├── **config.yaml**     
│&emsp;&emsp;├── BAM/  
│&emsp;&emsp;│&emsp;&emsp;├── *one_bam_by_sample.bam*  
│&emsp;&emsp;│&emsp;&emsp;└── *one_bai_by_bam.bai*  
│&emsp;&emsp;├── DESeq2/  
│&emsp;&emsp;│&emsp;&emsp;├── count_matrix.txt  
│&emsp;&emsp;│&emsp;&emsp;├── complete.txt  
│&emsp;&emsp;│&emsp;&emsp;├── up.txt  
│&emsp;&emsp;│&emsp;&emsp;├── down.txt  
│&emsp;&emsp;│&emsp;&emsp;└── Images/  
│&emsp;&emsp;├── fastqc/  
│&emsp;&emsp;│&emsp;&emsp;├── *one_html_by_sample.html*  
│&emsp;&emsp;│&emsp;&emsp;└── *one_zip_by_sample.zip*  
│&emsp;&emsp;├── htseqcount_CDS/  
│&emsp;&emsp;│&emsp;&emsp;└── *one_file_by_sample.txt*  
│&emsp;&emsp;├── htseqcount_fiveprime/  
│&emsp;&emsp;│&emsp;&emsp;└── *one_file_by_sample.txt*  
│&emsp;&emsp;├── htseqcount_threeprime/  
│&emsp;&emsp;│&emsp;&emsp;└── *one_file_by_sample.txt*  
│&emsp;&emsp;├── qualitativeAnalysis/  
│&emsp;&emsp;│&emsp;&emsp;├── bamDivision/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *one_bam_by_readsLength_by_sample.bam*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *one_bai_by_bam.bai*  
│&emsp;&emsp;│&emsp;&emsp;├── graphes/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── readsLengthRepartition/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *one_jpeg_by_sample.jpeg*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── periodicity/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;&emsp;&emsp;&emsp;├── *one_jpeg_by_readsLength_by_sample.start.jpeg*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;&emsp;&emsp;&emsp;└── *one_jpeg_by_readsLength_by_sample.stop.jpeg*  
│&emsp;&emsp;│&emsp;&emsp;├── readsLengthRepartition/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *one_file_by_sample.txt*  
│&emsp;&emsp;│&emsp;&emsp;└── periodicity/  
│&emsp;&emsp;│&emsp;&emsp;&emsp;&emsp;&emsp;├── *one_file_by_readsLength_by_sample.start.txt*  
│&emsp;&emsp;│&emsp;&emsp;&emsp;&emsp;&emsp;└── *one_file_by_readsLength_by_sample.stop.txt*  
│&emsp;&emsp;├── PROJECT_NAME.Analysis_report.txt  
│&emsp;&emsp;└── PROJECT_NAME.Final_report.html  
└── dag_all.svg    
- *logs* folder groups together all the output messages from tools used in RiboDoc analysis pipeline. Thus, in the event of an error, it allows you to identify the problematic step to give us feedback.   
- *RESULTS* folder contains 7 subfolders:   
&emsp;&emsp;i) *BAM*: it contains a BAM file for each sample (allows visualization on tools such as IGV).  
&emsp;&emsp;ii) *DESeq2*: it contains the count matrix, differential analysis tables and images related to the DESeq2 analysis.   
&emsp;&emsp;iii) *fastqc*: it contains raw data quality controls.   
&emsp;&emsp;iv) *htseqcount_CDS*: it contains htseq output for CDS counts.    
&emsp;&emsp;v) *htseqcount_fiveprime*: if the UTR covering option is on, this folder contains htseq output for five prime UTR counts.   
&emsp;&emsp;vi) *htseqcount_threeprime*: if the UTR covering option is on, this folder contains htseq output for three prime UTR counts.    
&emsp;&emsp;vii) *qualitativeAnalysis*: it contains all files related to qualitative test like periodicity and reads length repartition   
It contains also three files:  
&emsp;&emsp;i) *PROJECT_NAME.Analysis_report.html* gathers standard output of each analysis pipeline tool. It allows to know numbers of reads at each step a)raw reads b)reads after trimming and length selection c)after out RNA depletion d)after double alignment on the reference genome.  
&emsp;&emsp;ii) *PROJECT_NAME.Final_report.txt* presents all figures and explanation link with the differential analysis.  
&emsp;&emsp;iii) *config.yaml* to have a parameters backup.     

- *a dag file* which represents all the analysis pipeline steps with your samples.  

>Last big tip:  
In case that a sample is too variable against other replicats or if new samples sequencing are added to your study, you can delete/add them in the *fastq* subfolder, delete the subfolder *RESULTS/DESeq2*. Run again RiboDoc on the same *project_name* folder and it only (re)create missing files (complete analysis for added samples, new differential analysis with all samples available in *fastq* subfolder).  
