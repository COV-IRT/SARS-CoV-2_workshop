
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

Next, we will need to install two tools: samtools (for sam/bam file processing), and LoFreq (for variant calling).
```
conda activate workshop-env
conda install samtools=0.1.18
conda install lofreq
```

**Note:** the specified version of samtools is the one that currently works fine on MacOS Catalina. In case if you have a Linux instance (for example one from DNAnexus) you can ommit the `=0.1.18` portion of the `install samtools` command.

### SAM/BAM processing

If after read mapping you only have a SAM file follow the instructions below. If you have a *sorted* BAM file skip to the *Variant calling* subsection, otherwise if you have an unsorted BAM file proceed to step 2.

1. Let's convert our mapped reads from SAM to BAM format using samtools.
```
samtools view -S -b my_reads_file.sam > my_reads_file.bam
```

2. Now, we need to sort our BAM file as follows.
```
samtools sort my_reads_file.bam -o my_reads_file.sorted.bam
```

Finally, we can proceed to the variant calling.

### Variant calling

If you want to call both indels and SNPs then you will first have to assign indel quality scores in your BAM. LoFreq makes this easy to do for Illumina data, for other use cases you might have to consider an external tool to assign the indel quality scores before calling indels.

1. Since our dataset consists of paired-end Illumina reads we can use the built-in implementation of Dindel ([PMID 20980555](https://pubmed.ncbi.nlm.nih.gov/20980555/)) inside the LoFreq to assign the indel qualities as follows.
```
lofreq indelqual --dindel -f SARS-CoV-2-reference.fasta -o your_bam_file.sorted.indelqual.bam your_bam_file.sorted.bam
```

2. Now, we can call the variants using the following LoFreq command.
```
lofreq call -f SARS-CoV-2-reference.fasta --call-indels -o output_name.lofreq.indel.vcf your_bam_file.sorted.indelqual.bam
```

**Note:** The first step is optional if you don't plan to call indels. To just call SNPs you would omit the `--call-indels` flag to `lofreq call`, and you can skip the first stpe in that case.

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
