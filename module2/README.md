
#  Module #2: Read mapping and variant calling

***
Lead: Fritz J Sedlazeck
TA: ‪Medhat Mahmoud + Sairam Behera
Baylor College of Medicine
***

## Goals of this module
The goal of this module is to get you familiarized with the mapping of Illumina reads to the SARS-CoV-2 genome and to identify variations. For the latter we will focus on SNV (point mutations) and Structural Variations (SV). In addition to identfying these two types of variations we will further show you one way to assess the quality and summary statistics across the mapping, SNV and SV calling. 

### The main steps in this Module are:
1. Align the reads from Module 1 (BWA-mem)
2. Obtain alignment statistics (samtools)
3. Identifying single-nucleotide variations (LoFreq)
5. Identify Structural Variations (Manta)


The tools and methods we introduce here are clear snapshots of best practices for this analysis and there are multiple other methods available for your own anlysis.


## Learning Objectives
In the end of the Module 2 you should have a clear understanding how to map Illumina reads and identify SNV and SV together. In addition, you shoudl be able to identify problems and obtain preliminary insights into your data set at hand. 

We will further discuss standard file formats for this application: BAM and VCF files.

***

## Alignment and quality control of aligning short read data

## Mapping reads 

As discussed in the lecture the purpose of a mapping / alignment is to identify the most likely region a given read (ie. sequenced segment) was orginating from a given genome. 

First we will align the reads to the reference genome. 

Navigate to the folder and create a folder for your mapping results

```
cd Module2
mkdir mapping
cd mapping
```

Next we want to use the reads to start the BWA mem alignments. BWA mem is currently one of the standart short read based mapper and that is why we are using it here. 


The first step to starting read mapping --> context to DNAnexus

*Connect to your DNAnexus instance and open up a shell prompt.*

# Prepare for running the analysis

download all data using  
`dx download -r *`  

What data/ directories you can see now?  
`ls -l`  

Move to Module2  
`cd Module2`  

Make sure you add other conda branches:
```
conda config --add channels conda-forge
conda config --add channels bioconda
```

Create environement:  
`conda create -n SVanalysis python=2.7 bwa=0.7.17 manta=1.6.0 lofreq=2.1.5`

Activate the environement:  
`conda activate SVanalysis`


Export Manta  
`MANTA=/opt/conda/envs/SVanalysis/share/manta-1.6.0-0/bin/"`

Try it now:  
`python $MANTA/configManta.py`


## Time to run analysis
# **First Align Reads**
```
bwa mem
```

<!-- Accesws the data: FIXME
```
cd awesome_data

```
This dataset contains paired end reads. -->

<!-- Next: [module3!](module3.rst) -->
