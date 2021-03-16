
#  Module #3: Reference-guided assembly and low frequency variant calling

***
Leads: Yunxi Liu, Nicolae Sapoval, and Todd Treangen 
Rice University
***

## Goals of this module
* This hands on tutorial will teach you how to FIXME

## Learning Objectives
* Objective 1
* Objective 2
* Objective 3

***

## Background On reference-guided assembly and low frequency variant calling

## Low frequency variant calling

### Workspace configuration and software installation

*Connect to your DNAnexus instance and open up a shell prompt.*

First let's make sure we have `conda-forge` and `bioconda` channels set up in our conda.
```
conda config --add channels conda-forge
conda config --add channels bioconda
```

Next, we will need to install three tools: BWA MEM (for read mapping), samtools (for sam/bam file processing), and LoFreq (for variant calling). To make the process easier we will create a new conda environment and install all three tools into it.
```
conda create --name workshop-env python=3.8
conda activate workshop-env
conda install samtools=0.1.18
conda install lofreq
conda install bwa
```

**Note:** the specified version of samtools is the one that currently works fine on MacOS Catalina. In case if you have a Linux instance (for example one from DNAnexus) you can ommit the `=0.1.18` portion of the `install samtools` command.

### Obtaining sequencing read files

Stub

### Read mapping

Stub

### SAM/BAM processing

Stub

### Variant calling

Stub

## Reference-guided assembly

The first step to starting your assembly --> context to DNAnexus

*Connect to your DNAnexus instance and open up a shell prompt.*

Install XYZ
```
conda install -c bioconda -y XYZ
```

Accesws the data: FIXME
```
cd awesome_data

```
This dataset contains paired end reads. 

Next: [module4!](module4.rst)
