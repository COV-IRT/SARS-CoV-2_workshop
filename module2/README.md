
#  Module #2: Read mapping and variant calling

***
Lead: Fritz J Sedlazeck
Baylor College of Medicine
***

## Goals of this module
* This hands on tutorial will teach you how to:
1. Align the given reads using BWA-mem
2. samtools for aligning staistics
3. Identify Structural Variations using Manta
4. Detect single-nucleotide variations using LoFreq

## Learning Objectives
*   Identify Structural variations (SVs) (>=50bp)
*   Detect SNVs
*   Summary statistics for the alignment
*   How do SVs look like?
    *   Size distribution
    *   Types
*   SNVs
    *   How many we identified
    *   How many substitutions vs indels in SNVs?

***

## Background On Read mapping and variant calling

## Mapping reads and more

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
