# TERRA-QUANT
[![DOI](https://zenodo.org/badge/730564607.svg)](https://doi.org/10.5281/zenodo.15598540)

This repo is a RNA-seq pipeline that can quantify TERRA expression level from different RNA-seq datasets downloaded from GEO database. 

## Description
This pipeline is a combination of various bioinfomatic tools for performing reads dump, reads QC, alignment, total gene count and TERRA count. This pipeline inputs Sequence Read Archive (SRA) files and generates total gene and TERRA count tables for later expression analysis (For example, plotting heatmaps and boxplots).

## TERRA region and CHM13 gene annotation
The TERRA regions were defined based on RNA-seq reads mapping on CHM13 genome. Detail information is written in: CHM13_v2.0_add_TERRA_ITS_AS_v5.gtf.gz


3. For total human gene counts, we used gene annotation file provided by [CHM13 genome assembly](https://github.com/marbl/CHM13) (Original cytoband file was downloaded from this page as well.)

## Workflow

![TERRA_QUANT_pipeline](https://raw.githubusercontent.com/ls807terra/TERRA-QUANT/main/TERRA_QUANT_pipeline.svg)

## System Requirement

### Hardware Requirement
The minimum requirements that can perfform STAR aligner should be enough to execute all tools in this pipeline.

The requirements for STAR aligner are according to [STAR Github page](https://github.com/alexdobin/STAR).

1. Operation System: CentOS 7 or other later versions of Linux distrubutions.
2. CPU: x86-64 compatible processors
3. RAM: At least 16GB
4. Hard drive storage space: At lease 1TB

### Software Requirement
1. SRAToolkit (v2.11.1)
2. TrimGalore (v0.6.3)
3. Cutadapt (v2.3)
4. fastqc (v0.12.1)
5. STAR (2.7.9a)
6. samtools (v1.13)
7. deeptools (v3.3.1)
8. htseq-count (2.0.3)
10. Python version >= 3.8
11. R version >= 4.2


# Installation Guide 

## 1. Download this pipeline
Directly download from this link: [TERRA_RNA-seq_pipeline.zip](https://codeload.github.com/ls807terra/TERRA_RNA-seq_pipeline/zip/refs/heads/main)

or

Download by git clone command: `git clone https://github.com/ls807terra/TERRA_RNA-seq_pipeline.git`

## 2. Install required softwares
**Please download and install softwares on your system following their instructions.**

| Tool | version | Guide link |
|-------|:-----:|------:|
| SRAToolkit   | v2.11.1 | [Github page](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) |
| TrimGalore   | v0.6.3  | [Github page](https://github.com/FelixKrueger/TrimGalore) |
| Cutadapt     | v2.3    | [Github page](https://github.com/marcelm/cutadapt) |
| fastqc       | v0.12.1 | [Download page](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)|
| STAR         | 2.7.9a  | [Github page](https://github.com/alexdobin/STAR) |
| samtools     | v1.13   | [Github page](https://github.com/samtools/samtools) |
| deeptools    | v3.3.1  | [Official Website](https://deeptools.readthedocs.io/en/develop/) |
| htseq-count  | 2.0.3   | [Official Website](https://htseq.readthedocs.io/en/master/install.html) |

## 3. Setup conda enviroment
We created different conda enviroments for separating conflict tools.

To create a conda enviroment as same as ours, use the yaml file provided in the `install_env/`

**Two separated enviroment should be created!**

   1. `RNAseq_quantTERRA_env` for general RNA-seq analysis.
 
**_Command line_:**

   `conda env create -n "RNAseq_quantTERRA_env" -f RNAseq_quantTERRA_env.yml`

### Install R packages for `RNAseq_quantTERRA_env`

All the package you need for the RNA-seq analysis (mainly DEseq2) can be install by executing the `install_R_packages_RNAseq_quantTERRA.R`

To install, first activate the conda enviroment:

   `conda activate RNAseq_quantTERRA_env`

Then execute the **"install_R_packages_RNAseq_quantTERRA.R"** by Rscript.

   `Rscript install_R_packages_RNAseq_quantTERRA.R`


# Run this pipeline

### 1. Download SRA files from NCBI database

   First, you need to download the SRA files by SRAToolKit manually. Do this by the `prefetch` command.
   
   ```prefetch SRR_ID -o SRR_ID.sra```
   
   _**※ This is because NCHC nodes are not allowed for connecting to the Internet.**_
    
   _**※ Otherwise, it can be downloaded and exreacted to fastq by simply doing the command:**_
   
   ``` fastq-dump --split-files SRR_ID -O /path/to/output/dir/ ```

### 2. Edit the pipeline configure file

   Next, configure the `c.0_RNAseq_QuantTERRA.cfg` file.

   Edit the software path and other settings in the `c.0_RNAseq_QuantTERRA.cfg` file.

   After editing settings, run `sh 0_Configure_Setting.sh c.0_RNAseq_QuantTERRA.cfg`.

   **All pipeline scripts, `c.0_RNAseq_QuantTERRA.cfg` and `0_Configure_Setting.sh` should be under the same directory.**

### 3. Execute pipeline (This will run on computation nodes)

   **_A series of jobs, including all steps in the pipeline, will be sent to computation nodes of NCHC Taiwania3 by Slurm queueing System._**
   
   Finally, sent a single SRA file to execute this pipeline by the following commands:**

   *Usage:* `sh 0.0_ps_Pipeline_v3.sh SRR_ID`

   *Example:* `sh 0.0_ps_Pipeline_v3.sh SRRSRR3304509`

# Demo and Data reproduction
Here is a demonstration of using this pipeline, you can start from downloading SRA files from NCBI database or from processing fastq raw files. 

In the **_Demo_** folder, ```demo_HGPS_sra_list.txt``` records a series of SRA accession number that were used in the TERRA quantification of Figure 6E in our paper. The corresponding annotation file, ```HGPS_annotation.xlsx``` , for these sample is also desposited in **_Demo_** folder. This demo describes steps to reproduce the result.

### Start from downloading SRA files
First, download these data by SRAToolkit ```prefetch``` command or execute the script provided in this repository 

  **_command:_** ```sh a.0_ps_SRA_download.sh Demo/SRA list.txt SRA```

This will download all SRA listed in the text file to **_SRA_** folder then downloaded SRA files can be fed into the pipeline. 

**Please named the folder as ```SRA``` and store SRA file in it to run the pipeline, current version of script does not recognize SRA folder in other name!**

Next, you can execute the pipeline script ```0.0_ps_Pipeline_v4.sh``` to process the SRA file one by one, for example:

  **_command:_** ```sh 0.0_ps_Pipeline_v4.sh SRR7689162```

or execute by a loop:

   **_command:_** 
   ```
   while read ln
    do
      sh 0.0_ps_Pipeline_v4.sh $ln
    done < SRA list.txt
   ```
## Expect outputs
The raw data will process through the steps show in the **_Workflow_** diagram.

Outputs from each step will be store at individule folder:

| folder | file | description |
|-------|:-----:|------:|
| fastq    | *.fastq.gz | fastq files extract from a SRA file |
| trimed_fq| *.fq.gz    | Adaptor and quailty trimmed fastq and QC files |

# Configure file format

**The format of the configure file is as follow:**

   ```
   ## Genome
   CHM13GENOME:/staging/biology/ls807terra/0_genomes/star_index/CHM13_human
   
   ## SLURM setting
   PARTITION:186
   MEMORY:186
   NCORE:28
   
   ## Pipeline workdir, Please specifiy this!
   WORKDIR:../
   
   # Human CHM13 GTF and TERRA regions
   QGTF:/staging/biology/ls807terra/0_genomes/genome_gtf/CHM13/CHM13_v2.0.gtf
   
   # CHM13 cytoband file with TERRA regions for telomerehunter.
   CYTOBAND:/staging/biology/ls807terra/0_bedfiles/hTERRA/chm13v2.0_cytobands_allchrs_Add_TERRA_v3.bed
   
   ## Software path
   
   # SRAToolKit fastq-dump
   FASTQDUMP:/opt/ohpc/Taiwania3/pkg/biology/SRAToolkit/sratoolkit_v2.11.1/bin/fastq-dump
   
   # SRAToolKit prefetch
   PREFETCH:/opt/ohpc/Taiwania3/pkg/biology/SRAToolkit/sratoolkit_v2.11.1/bin/prefetch
   
   # Trimgalore
   TRIMGALORE:/opt/ohpc/Taiwania3/pkg/biology/TrimGalore/TrimGalore_v0.6.3/trim_galore
   
   # Cutadapt
   CUTADAPT:/opt/ohpc/Taiwania3/pkg/biology/Cutadapt/Cutadapt_v2.3/bin/cutadapt
   
   # samtools
   SAMTOOLS:/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools
   
   # STAR alignment
   STARALIGN:/opt/ohpc/Taiwania3/pkg/biology/STAR/STAR_v2.7.9a/bin/Linux_x86_64/STAR
   
   # deeptools bamcoverage
   BAMCOV:/opt/ohpc/Taiwania3/pkg/biology/deepTools/deepTools_v3.3.1/bin/bamCoverage
   
   # UCSC tools bigwig to bedgraph
   BW2BGD:/staging/biology/ls807terra/0_Programs/UCSC_tools/bigWigToBedGraph
   
   # HTseq-count
   HTSEQ:/staging/biology/ls807terra/0_Programs/anaconda3/envs/RNAseq_quantTERRA/bin/htseq-count

   
   # BAM tools for BAM merge
   BAMTOOL:/staging/biology/ls807terra/0_Programs/bamtools/build/bin/bamtools
   ```
### Configure file's Rules
   
1. Edit to change the softwares' paths, genome or Slurm settings, following the format of TOOL_NAME:/path/to/Software , separated by ":"
  
2. To register a new tool in configure file, Add a new line follow the format.

3. To add a new tool and its path to pipeline:

   3.1 In the script of a new software (file name should follow this format: *_ps_*.sh), write:
    	   ''' software="TOOL_ID"
    	      # To use the software
    	      ${software} -options [files] ... '''

   3.2 In the configure file, write:
    	   ''' TOOL_ID:/path/to/this/tool

   3.3 Run the script: ```sh 0_Configure_Setting.sh this_configure_file.cfg```
   	   This will reconginze the "TOOL_ID" in a script and replace it by "/path/to/this/tool"

4. To comment or annotate, please add a "#" at the start of a line.
   
5. Please update the software path if they were changed.
   
6. Make sure that all pipeline script (file name: *_ps_*.sh) are in the same directory with 0_Configure_Setting.sh and this configure file.

**Notice!**

**Your SRA file need to be locate at ${workdir}/SRA!**

**Check the configure file `c.0_RNAseq_QuantTERRA.cfg`.**

**The workdir is usually at "../" .**


# Potential Problems during installation and run time

### **Notice! Very important!**

You will have some error when installing the R package `devtools` and `xml2`.

For solving  `xml2`, you simply need to install `libxml2` by conda, which is automaticlly installed by the yaml file.

To handle `devtools` problem, you need to fix the libpng in the conda enviroment `lib/`

### Problem1 libpng15 error:

   libpng15 not found, when installing devtools and some others.
   
   It possible that you have already install libpng, but R can not find it.
   
   Because it usually recongizes the libpng15.so, and you probably have another version. 

**Solution:**

   Before running install_R_packages_RNAseq_quantTERRA.R
   
   1. Make sure that there is libpng installed in your conda enviroment.
   
    conda install libpng 
   
   2. Find it under /your/path/to/anaconda3/env/your_env/lib/
   
    For example. I have these file under my env path 
    
    ls ~/anaconda3/env/TEST_env/lib/ | grep png
    
    libpng16.a
    libpng16.so
    libpng16.so.16
    libpng16.so.16.39.0
    libpng.a
    libpng.so
   
    The version of my libpng is 16 not 15, so R can not find it.
    
   3. Create a softlink that named "libpng15.so" but were direct to libpng16.so
   
    Build a softlink 
   
    ln -s libpng16.so.16.39.0 libpng15.so
    
   4. Now R should find libpng15 properly.

Make sure the libpng were fixed. Then you can execute `install_R_packages_RNAseq_quantTERRA.R` by Rscript:

   `Rscript install_R_packages_RNAseq_quantTERRA.R`

*This will run a verrrrry long time.*
