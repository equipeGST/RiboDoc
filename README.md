# Step by step tutorial

Welcome to the RiboDoc tool tutorial !  

>RiboDoc is designed to perform all classical steps of **ribosome profiling** (RiboSeq) data analysis from the FastQ files to the differential expression analysis with necessary quality controls.

If you want to easily understand how to launch RiboDoc on your own computer, you can check our video tutorial just here :
[![RiboDoc_Video](https://github.com/equipeGST/RiboDoc/assets/75135539/51fdeaf2-8a7c-4cd8-a5b3-d58d3c3adf5e)](http://www.youtube.com/watch?v=e9_SFz_YEK0 "Tutorial - RiboSeq analysis with RiboDoc pipeline")

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
&emsp;&emsp;&emsp;***biological_condition_name.replicate_number.fastq.gz***     
For example, a replicate of the wild-type condition the sample could be named *Wild_Type.56.fastq.gz* and the name of a replicate for the mutant samples could be *Mutant.42.fastq.gz*.    
**Please avoid having dashes "-", parentheses "(" and/or ")", blank spaces " " and more generally any special characters except underscores "_" in file names, as it may break RiboDoc pipeline**    
For example, avoid having a name like *Mutant(actin).1.fastq.gz*, *fastq from 2020.1.fastq.gz*, etc... And prefer using underscores like *My_awesome_example.1.fastq.gz*    

>Caution, for **Windows**, extensions can be hidden.    

Example of folder architecture at this step:  
Project_name  
└── fastq   
&emsp;&emsp;&emsp;├── Wild_Type.1.fastq.gz   
&emsp;&emsp;&emsp;├── Wild_Type.2.fastq.gz   
&emsp;&emsp;&emsp;├── Mutant.1.fastq.gz   
&emsp;&emsp;&emsp;└── Mutant.2.fastq.gz   

If you want to try RiboDoc on an example dataset, you can find one on GEO : [GEO GSE173856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173856). If you do so, be aware that this dataset does not align enough material for a use of the riboWaltz pipeline and TRiP should be selected (see later for details).  

### b) *database*  
In this subfolder, you must put at least the following three files:  
- Your reference genome fasta file: Whether it's the genome or the transcriptome, it must be your reference fasta file to wich the reads will be aligned. We advise you to download it from the [Ensembl](https://www.ensembl.org/index.html) database, as the files are maintained up to date and following the standard GFF format.  
For example, for an entire yeast genome, you can look for [S. cerevisiae genome](https://http://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index), then click on "*Download DNA sequence (FASTA)*" and download the *Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz* file in the genome assembly list of fasta files, decompress it and place it in the *database* subfolder of your project directory.  

- A GFF format annotation file corresponding to the reference genome.  
For example, for the yeast genome annotations, you can look for [S. cerevisiae genome](https://http://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index), then click on "*Download GFF3*" and download the *Saccharomyces_cerevisiae.R64-1-1.110.gff3* file (the '110' part would change depending on the version of annotation) in the genome annotation list of gff files, unzip it and place it in the *database* subfolder of your project directory.  

- out-RNA fasta file: This file must gather together DNA sequences you want to remove from the analysis. As a rule, these are at least ribosomal RNA sequences (rRNA). You can also add mitochondrial DNA, non-chromosomal DNA or any other fasta sequence of your choice which you want to be discarded in the analysis. You can generate the fasta files associated to the orgnanism rRNAs from the [Silva](https://www.arb-silva.de) database for example. You could also find your desired sequences in the database from [NCBI website](https://www.ncbi.nlm.nih.gov/) or design them yourself. It is up to you. If you want to remove some specific sequences, you just have to create a file with, for each sequence, one line starting with a ">" where you can add a name for your sequence followed by one line containing the sequence to remove. It gives you this format :  
>&emsp; > Sequence_X
>&emsp; GCTGACACGCTGTCCTCTGGCGACCTGTCGTCGGAGAGGTTGGGCCTCCGGATGCGCGCGGGGCTCTGGCCTCACGGTGACCGGCTAGCCGGCCGCGCTCCTGCCTTGAGCCGCCTGCCGCGGCCCGCGGGCCTGCTGTT  
>&emsp; > Sequence_Y
>&emsp; CTCTCGCGCGTCCGAGCGTCCCGACTCCCGGTGCCGGCCCGGGTCCGGGTCTCTGACCCACCCGGGGGCG  

If you do not want to remove any sequence, just leave the outRNA parameter of the configuration file empty (fasta_outRNA: "") **but we strongly advise to remove rRNA sequences**.

Example of folder architecture at this step:  
Project_name  
├── fastq   
│&emsp;&emsp;├── Wild_Type.1.fastq.gz   
│&emsp;&emsp;├── Wild_Type.2.fastq.gz   
│&emsp;&emsp;├── Mutant.1.fastq.gz   
│&emsp;&emsp;└── Mutant.2.fastq.gz  
└── database   
&emsp;&emsp;&emsp;├── reference_genome_sequences.fa  
&emsp;&emsp;&emsp;├── reference_genome_annotations.gff3  
&emsp;&emsp;&emsp;└── RNA_to_remove.fa  

### c) [config.yaml](https://raw.githubusercontent.com/equipeGST/RiboDoc/main/config.yaml) file  
The *config.yaml* file allows you to define some parameters to tell RiboDoc which data you want to process and how.  
You must download it [here](https://raw.githubusercontent.com/equipeGST/RiboDoc/main/config.yaml) and open it with a text editor as a text file.    
It must be carefully completed and be present in the project directory everytime you want to run RiboDoc. A copy of this file will be made in the *RESULTS/* folder to keep a trace of the parameters you chose for your last analysis.   

>Caution  
>&emsp;&emsp;&emsp;Spaces and quotation marks **must not be changed** ! Your information must be entered between quotes and should NOT have spaces     


####How to fill the configuration file :
##### Project name  
First and easy step, the project name ! You can use the same as your folder.  
*project_name*: "Project_title"  

##### Name of database files  
You must enter the full name **with extensions** without the path of files added in the database subfolder previously created.   
*fasta*: "reference_genome_fasta_file.fa"  
*gff*: "corresponding_GFF_annotation_file.gff3"  
*fasta_outRNA*: "unwanted_DNA_sequences_fasta_file.fa" (can be empty if no specific RNA sequence is to be removed from the analysis)  

##### Trimming information  
During RiboDoc process, reads/RPF are trimmed and selected depending on their length.   
*already_trimmed*: If your data contains reads already trimmed of their adapter, you can set this option on “yes”. Else, set it on "no".   
*adapt_sequence*: If they are not trimmed, you should specify the sequence of the adapter in quotes on the line here like "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA". If you do not put anything between the quotes, RiboDoc will try to fond the adapter itself but this can sometimes lead to a wrong adapter sequence then an error later in the pipeline. When this parameter is left empty, be careful and check on what is written in the "RESULTS/adapter_lists/" files.

##### RPF lengths selection  
You also have to define the range for read length selection. Default values select reads from 25 to 35 nucleotides long as RPF are usually around 30 nucleotides long.  
*readsLength_min*: Minimum reads length. Default : "25"   
*readsLength_max*: Maximum reads length. Default : "35"   

##### Annotations vocabulary  
You might also need to specify features keywords in the GFF file to fit your GFF file format :   
*gff_cds_feature*: Feature corresponding to CDS/ORF in the annotation file (is usually CDS or ORF). Default : "CDS".     
*gff_mRNA_feature*: Feature corresponding to mRNA/transcript in the annotation file (is usually mRNA or transcript). Default : "mRNA".     
*gff_5UTR_feature*: Feature corresponding to 5'UTR in the annotation file (is usually five_prime_UTR or 5UTR). Default : "five_prime_UTR".     
*gff_parent_attribut*: Attribut to regroup reads during counting (is usually Parent (stands for Parent of the CDS features) to regroup counts by transcripts for the differential analysis). Default : "Parent".     
*gff_name_attribut*: Name of the genes features in the GFF (is usually Name or gene_name). Default : "Name".     

##### Differential analysis settings  
To be able to perform statistical analyses, you must define a reference condition as well as your thresholds.   
*reference_condition*: Corresponds to the reference biological_condition_name in your FastQ files name. Ex : "Wild_Type" (as in *Wild_Type.1.fastq.gz*).
*p-val*: p-value threshold for the differential analysis. Default : 0.01.  
*logFC*: logFC threshold for the differential analysis. Defaut : 0 (keep all the genes without logFC filtering).  

##### Window for qualitative test  
During the qualitative analysis, the periodicity is observed on nucleotides around start and stop codons.     

>The periodicity is displayed as a metagene profile (metaprofile). It provides the amount of footprints  on relative coordinates around annotated start and stop codons.   

##### Trimming information  
2 pipelines dedicated to quality controls are available in RiboDoc. The first one uses the [riboWaltz tool](https://github.com/LabTranslationalArchitectomics/riboWaltz) which can need high RAM resources depending on your input data volume and reference organism. The second pipeline is a series of scripts called TRiP. The main difference between the two is that riboWaltz calculates a P-site offset (length from the beginning of the read to the first base of its P-site) for each length of RPF. It allows a more precise determination of the ribosome decoding sites coordinates and more accurate metaprofiles whereas TRiP makes the metaprofiles with the coordinates of the first nucleotide of the reads (5' end), which do not really correspond to the codons decoded by the ribosome.
*qualitative_analysis*: Choose between the 2 qualitative analysis pipelines. Default : "ribowaltz"
The window selected by default is -30/+90 nucleotides and -90/+30 nucleotides around start and stop codons respectively.   
*window_utr*: Define your window before start and after stop. Default is "30"  
*window_cds*: Define your window after start and before stop. Default is "90"  

######################
###### Optional ######
######################
##### UTR covering option  
*UTR*: Put "yes" if you want to count reads aligned on UTR regions. (Default value is "no")
*gff_3UTR_feature*: Feature corresponding to 3'UTR in the annotation file. Useless if *UTR: "yes"*. "three_prime_UTR" is the default value (can sometimes be "3UTR").     
################################################################  


Example of folder architecture at this step:  
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
To get RiboDoc, open a terminal. We know this can look scary for some but do not worry, there are only 2 easy lines to write and then you are done !  
If you never used RiboDoc on your workstation, you must pull it from Docker Hub.  

### Use of singularity    
If you want to use a singularity image instead of a docker image (for the use on a cluster for example), you must pull it from Docker Hub. You can name the image as you want but the extension should be *.sif* :    
&emsp;&emsp;&emsp;`singularity pull /path/to/your/singularity/image.sif docker://equipegst/ribodoc`    
### Use of docker
Copy and paste the following command line:    
&emsp;&emsp;&emsp;`docker pull equipegst/ribodoc`  
If you have any error, it might come from a rights/permission problem as docker needs administrator privilegies to work. If you do have admin privilegies, you should try to copy and paste this command:    
&emsp;&emsp;&emsp;`sudo docker pull equipegst/ribodoc`     


## 4) Run RiboDoc  
Now that your folder's architecture is ready and the image is on your computer, it's time to start ! If you have pulled RiboDoc, you can run it with the following command:  
### Use of singularity  
&emsp;&emsp;&emsp;`singularity run --bind /path/to/your/project/folder/:/data/ /path/to/your/singularity/image.sif CPU_NUMBER MEMORY_AMOUNT`  
### Use of docker  
&emsp;&emsp;&emsp;`docker run --rm -v /path/to/your/project/folder/:/data/ equipegst/ribodoc CPU_NUMBER MEMORY_AMOUNT`  
*/path/to/your/project/folder/* corresponds to the **full path** of the directory with all the files you prepared for the analysis (usually starting with a "/"). To get this full path, you can usually drag and drop your folder to the terminal to display its path, or find it in the folder properties.    
*:/data/* **must not be modified in any way** as it corresponds to the path of the project directory inside the container.   
*CPU_NUMBER* is an integer/number corresponding to the total number of threads you want to use for your analysis (Ex : 4).   
*MEMORY_AMOUNT* is an integer/number corresponding to the maximum amount of memory in Gigabytes (RAM) you want to use for your analysis (Ex : 30). Be careful to have enough memory available for alignment steps or the pipeline will be stopped without a specific error written in the "logs" folder. If the pipeline crashes without easily findable error, try lowering this value to enable less jobs at the same time.    
>You can also add the "--disable-cache" option after "singularity run" in order to delete singularity image and cache after the run is done.  

>Caution:   
>&emsp;&emsp;&emsp;The path to your project folder and the "/data/" path must start and finish with a slash "/".   


>For Windows users  
>&emsp;&emsp;&emsp;The path to your project folder has to start at the local disks C or D: C:\your\path\   
>&emsp;&emsp;&emsp;This path to your project folder has to be composed and finished with backslashes "\" (instead of slashes "/")   
>&emsp;&emsp;&emsp;/data/ path does not change in any way !   

Once this is written, press 'Enter' and your job is done ! You can leave and chill while RiboDoc does the rest for you.  

## 5) In case of any error   
Managed by snakemake, the pipeline will finish all jobs unrelated to the rule/job that failed before exiting. You should let it finish to be sure to avoid any error in future runs but you can still force the container to stop with Docker Desktop or with the following command lines (might need the "sudo" keyword at the beginning for docker) :  
>&emsp;> docker container ls   

Which provides you the container's ID (For example : 9989909f047d), then :   
>&emsp;> docker stop ID  
Where "ID" is the id obtained with the previous command

With singularity you can just terminate the job by pressing Ctrl+C **once** and snakemake will delete all files wich might be corrupted before stopping. If you press it again, this might lead to an error in future runs.
If you have issues with the use of Docker, you must refer to their [website](https://docs.docker.com/).     
If the error happens during the use of RiboDoc, the rule (job) which failed is indicated in your terminal. You can then find the error outputs in the *logs* folder. Each rule have a precise name and a folder related to it with files corresponding to the different steps of this rule.  
In most cases, a problem occurs because the memory provided is insufficient and no precise error message is shown in the *logs* folder. This usually happens during the index builds (for the alignement tools), the alignments or the riboWaltz tool script, as they are the steps asking for the more resources in RiboDoc pipeline. If you cannot find why the pipeline stopped, try resuming it with more memory available if you are on a cluster and/or reduce the *MEMORY_AMOUNT* in your launch command to reduce potential multi-threading.   
If you can not solve the problem by yourself, you can contact us through the ["issues"](https://github.com/equipeGST/RiboDoc/issues) part of our github or by mail. We will be happy to provide all the help we can.  
>If you want to send us a request because of an error, the easiest way for us to help you is if you send us an archive containing your "config.yaml" file, "logs/" folder and "stats/" folder so we can help you finding what happened and why.

## 6) Understand the results  
Here is the example folder architecture after a RiboDoc run.  
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
│&emsp;&emsp;└── **RNA_to_remove.fa**  
├── **config.yaml**   
├── logs/   
│&emsp;&emsp;└── *one_log_folder_by_job*  
├── stats/   
│&emsp;&emsp;└── *one_file_by_steps_of_interest_for_alignment_stats*  
├── MAIN_RESULTS/    
│&emsp;&emsp;├── *project_name.Analysis_Report.readsLength_min-readsLength_max.txt*     
│&emsp;&emsp;├── *project_name.Analysis_Table_summary.readsLength_min-readsLength_max.csv*    
│&emsp;&emsp;├── Quantitative_analysis/
│&emsp;&emsp;│&emsp;&emsp;├── DESeq2_CDS.readsLength_min-readsLength_max/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *count_matrix_by_transcript_IDs.csv*   
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *names_correspondence_list.csv*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── DESeq2_by_gene/ and DESeq2_by_transcript/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *count_matrix.txt*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *complete.txt*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *up.txt*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *down.txt*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *project_name.Final_report.html*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── Images/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *all_quantitative_analysis_graphs.png*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *all_quantitative_analysis_graphs.tiff*  
│&emsp;&emsp;└── Qualitative_analysis/
│&emsp;&emsp;│&emsp;&emsp;├── *frame_psite.tiff*
│&emsp;&emsp;│&emsp;&emsp;├── *frame_psite_length.tiff*
│&emsp;&emsp;│&emsp;&emsp;├── *region_psite.tiff*
│&emsp;&emsp;│&emsp;&emsp;└── periodicity_-window_utr+window_cds/  # If you chose *qualitative_analysis: "ribowaltz"* in the "config.yaml" file  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── transcriptome_elongated.one_folder_by_sample/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *one_file_by_read_length_start.tiff*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *one_file_by_read_length_stop.tiff*  
├── RESULTS/  
│&emsp;&emsp;├── *config.yaml*     
│&emsp;&emsp;├── *project_name.Analysis_Report.readsLength_min-readsLength_max.txt*     
│&emsp;&emsp;├── *project_name.Analysis_Table_summary.readsLength_min-readsLength_max.csv*     
│&emsp;&emsp;├── adapter_lists/  
│&emsp;&emsp;│&emsp;&emsp;└── *one_file_by_sample.txt*     
│&emsp;&emsp;├── annex_database/  
│&emsp;&emsp;│&emsp;&emsp;├── *gff_features_counts.txt* (summary of features present in GFF file)     
│&emsp;&emsp;│&emsp;&emsp;├── *NamedCDS_genomic_gff_files.gff*     
│&emsp;&emsp;│&emsp;&emsp;├── *transcriptome_elongated.nfasta* (only CDS+elongation sequences)     
│&emsp;&emsp;│&emsp;&emsp;├── *transcriptome_elongated.gff* (only CDS+elongation annotations)     
│&emsp;&emsp;│&emsp;&emsp;├── *transcriptome_elongated.exons_gtf_file_for_riboWaltz.gtf*     
│&emsp;&emsp;│&emsp;&emsp;├── *transcriptome_elongated.exons_fasta_file_for_riboWaltz.fa*     
│&emsp;&emsp;│&emsp;&emsp;└── index_files/
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *index_files_for_bowtie2.bt2*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *index_files_for_aligners.ht2*  
│&emsp;&emsp;├── BAM.readsLength_min-readsLength_max/  
│&emsp;&emsp;│&emsp;&emsp;├── *one_bam_by_sample.bam*  
│&emsp;&emsp;│&emsp;&emsp;└── *one_bai_by_bam.bai*  
│&emsp;&emsp;├── BAM_transcriptome.readsLength_min-readsLength_max/ (for riboWaltz use)  
│&emsp;&emsp;│&emsp;&emsp;├── *one_bam_by_sample.bam*  
│&emsp;&emsp;│&emsp;&emsp;└── *one_bai_by_bam.bai*  
│&emsp;&emsp;├── DESeq2_CDS.readsLength_min-readsLength_max/  
│&emsp;&emsp;│&emsp;&emsp;├── *count_matrix_by_transcript_IDs.csv*   
│&emsp;&emsp;│&emsp;&emsp;├── *names_correspondence_list.csv*  
│&emsp;&emsp;│&emsp;&emsp;└── DESeq2_by_gene/ and DESeq2_by_transcript/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *count_matrix.txt*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *complete.txt*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *up.txt*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *down.txt*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *project_name.Final_report.html*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── Images/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *all_quantitative_analysis_graphs.png*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *all_quantitative_analysis_graphs.tiff*  
│&emsp;&emsp;├── fastqc/  
│&emsp;&emsp;│&emsp;&emsp;├── fastqc_before_trimming/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *one_html_by_sample.html*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *one_zip_by_sample.zip*  
│&emsp;&emsp;│&emsp;&emsp;└── fastqc_after_trimming/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *one_html_by_sample.html*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *one_zip_by_sample.zip*  
│&emsp;&emsp;├── HTSeq-counts/  
│&emsp;&emsp;│&emsp;&emsp;└── htseqcount_CDS/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *count_matrix_by_transcript.csv*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *one_file_by_sample.txt*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *one_zip_by_sample.zip*  
│&emsp;&emsp;├── periodicity_-window_utr+window_cds/  # If you chose *qualitative_analysis: "ribowaltz"* in the "config.yaml" file  
│&emsp;&emsp;│&emsp;&emsp;└── transcriptome_elongated.one_folder_by_sample/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *read_length_start.tiff*  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;└── *read_length_stop.tiff*  
│&emsp;&emsp;├── riboWaltz.readsLength_min-readsLength_max/  # If you chose *qualitative_analysis: "ribowaltz"* in the "config.yaml" file  
│&emsp;&emsp;│&emsp;&emsp;└── *riboWaltz's qualitative analysis results*  
│&emsp;&emsp;└── qualitativeAnalysis/  # If you chose *qualitative_analysis: "trip"* in the "config.yaml" file  
│&emsp;&emsp;│&emsp;&emsp;└── *TRiP's qualitative analysis results*  
├── *dag_all.svg*    
└── *dag_last_run.svg*    

- The *MAIN_RESULTS/* folder contains a copy of the most interesting files present in the *RESULTS/* folder.    
- The *RESULTS/* folder contains these subfolders:   
&emsp;&emsp;i) *adapter_lists/*: it contains a text file with the adapters list for each sample that were found in the "config.yaml" file or determined from data is the user did not provide any adapter sequence in the "config.yaml" file.  
&emsp;&emsp;ii) *annex_database/*: it contains the indexes for the alignments, the re-formatted fasta and gff for the analysis and the gtf for the riboWaltz pipeline.  
&emsp;&emsp;iii) *BAM.readsLength_min-readsLength_max/*: it contains a BAM file for each sample (allows visualization on tools such as IGV).  
&emsp;&emsp;iiib) *BAM_transcriptome.readsLength_min-readsLength_max/*: it contains a BAM file for each sample generated from the transcriptome/exome gtf and fasta files generated for riboWaltz.  
&emsp;&emsp;iv) *DESeq2_CDS.readsLength_min-readsLength_max/*: it contains the differential analysis html report (*project_name.Final_report.html*), the count_matrices, the tables and the images related to the differential analysis, regrouped by gene or by transcript.   
&emsp;&emsp;v) *fastqc/*: it contains data quality controls.   
&emsp;&emsp;vi) *HTSeq-counts/*: it contains HTseq output for CDS counts (and UTR regions if selected in the "config.yaml" file).    
&emsp;&emsp;vii) *qualitativeAnalysis/* or *riboWatlz.readsLength_min-readsLength_max/* and *periodicity_-window_utr+window_cds*: they contain all files related to qualitative test like metaprofiles and reads lengths repartition   
It contains also three files:  
&emsp;&emsp;i) *project_name.Analysis_Report.readsLength_min-readsLength_max.txt* gathers standard output of each analysis main tool. It allows to know numbers of reads at each step a)raw reads b)reads after trimming and length selection c)after out RNA depletion d)after double alignment on the reference genome.  
&emsp;&emsp;ii) *project_name.Analysis_Table_summary.readsLength_min-readsLength_max.csv* summarizes the same standard outputs as previous file in a table.  
&emsp;&emsp;iii) *config.yaml* to have a parameters backup.     
- The *logs/* folder groups together all the error output messages from tools used in RiboDoc analysis pipeline. Thus, in the event of an error, it allows you to identify the problematic step to give us feedback.   
- The *stats/* folder groups the main tools statistics which the *Analysis_Report* and *Analysis_Table_summary* files are made from.   
- The *dag file* which represents the analysis steps with your samples.  

>Last big tips:
In case a sample is too variable against other replicates or if new sequenced samples are to be added to your study, you can delete/move or add them in the *fastq* subfolder. RiboDoc will only process necessary steps based on the fastq files list.  
Also, if the pipeline crashes but you cannot fin any obvious reason why, it usually means that it is due to a lack of memory. Re-run the pipeline specifying less MEMORY_AMOUNT is the command line before destroying your computer, or try using it with more resources if you are on a cluster ;-)  
Containers state saves and cache can take more and more space if you do not use the "--rm" option for docker. To clean everything related to containers, just right in your terminal :
>&emsp;> docker system prune -af

Thank you for using RiboDoc !   
We wish you the best results for your analyses !  
