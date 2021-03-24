
#  Module #2: Read mapping and variant calling

***
Lead: Fritz J Sedlazeck
Baylor College of Medicine
***

## Goals of this module
* This hands on tutorial will teach you how to:
1. Align the given reads using BWA-mem
2. Identify Structural Variations using Sniffles
3. Detect single-nucleotide variations using LoFreq

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

Install BWA aligner
```
conda install -c bioconda -y bwa
```

Accesws the data: FIXME
```
cd awesome_data

```
This dataset contains paired end reads.

Next: [module3!](module3.rst)
