
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

An example of BAM file with indel qualities for the SRR12447392 sample of SARS-CoV-2 visualized using IGV browser.
![](Figures/IGV-SRR12447392-indel-quals.png)

## Reference-guided assembly

The first step to starting your assembly --> context to DNAnexus

*Connect to your DNAnexus instance and open up a shell prompt.*

The tools we need to perform a reference guided assembly are bcftools and lofreq. First we still want to make sure that we have `conda-forge` and `bioconda` channels set up in our conda.
```
conda config --add channels conda-forge
conda config --add channels bioconda
```
Next, activate the working environment using 
```
conda activate workshop-env
```
Install tools using 
```
conda install bcftools
```
If you have already install lofreq in the variant calling step, you can skip the following command
```
conda install lofreq
```
Bcftools contains a collection of tools for variant calling and manipulating VCFs and BCFs, the specific command we are interested in is `bcftools consensus`, which will create consensus sequence by applying VCF variants to a reference fasta file. By default, the program will apply all alternative variants no matter what the allele frequencies are. To generate a valid reference-guided assembly, first we want to filter the variants by their allele frequencies (and some other characteristics if additional quality control is required). 

### Filtering Variants by the allele frequencies

In the previous step, we have generated the vcf file for the targeted alignment file (BAM file). Such vcf file will contain all variants with allele frequencies ranging from 0 to 1. By default lofreq would perform some default filtering to make sure the variant calls are accurate. For reference-guided assembly, we are going to mainly focus on variants with allele frequencies >= 0.5, which indicates that more than half of reads supports the allele base at certain location. We use the following command with option `-a` or `--af-min` to filter vcf files by the allele frequencies of the variants
```
lofreq filter -a 0.5 -i your_input.vcf -o your_output.filtered.af50.vcf
```

### Construction of consensus sequence
Now we can generate the consensus sequence using bcftools, currently bcftools only supports compressed gz file format, so we have to first compress the vcf file, and then index it.
```
bgzip your_output.filtered.af50.vcf
bcftools index your_output.filtered.af50.vcf.gz
cat SARS-CoV-2-reference.fasta | bcftools consensus your_output.filtered.af50.vcf.gz > your_output.consensus.fasta
```

Next: [module4!](module4.rst)
