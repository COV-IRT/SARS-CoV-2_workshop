
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

### Mapping reads 

As discussed in the lecture the purpose of a mapping / alignment is to identify the most likely region a given read (ie. sequenced segment) was orginating from a given genome. 

First we will align the reads to the reference genome. 

Navigate to the folder and create a folder for your mapping results

```
cd Module2
mkdir mapping
cd mapping
```

Next we want to use the reads to start the BWA mem alignments. BWA mem is currently one of the standart short read based mapper and that is why we are using it here.  To see the available options of bwa just execute its commad:

```
bwa
```

We first need to index the reference genome itself by using bwa index:

```
bwa index reference.fasta 
```
This should only take a few seconds since the SARS-Cov-2 genome is very small.


Next we are ready to start mapping the fastq reads to the genome itself. For this we want to use the bwa mem option that is best capable to handle the Illumina paired end reads:

```
 bwa mem -t 2 reference.fasta SRR12447392_1.fastq SRR12447392_2.fastq > our_mapped_reads.sam
```

This executes bwa mem with 2 threads (-t parameter) give our previously indexed reference.fasta and the two fastq files representing the Illumina paired end reads. 

After a few seconds the program ends and we have our first result as : our_mapped_reads.sam. This is a standard text file and we can take a look. As highlighted in the lecture we have a header in this file indicated with @ and then entries per read per line. 

For subsequent analysis we need to compress (SAM -> BAM) the file. For this we are using samtools with the option: view

```
samtools view -hb our_mapped_reads.sam > our_mapped_reads.bam
```

The options -h ensures that the header is kept for the output file and the option -b tells samtools that we want to obtain the compressed (BAM) version. 
Next we need to sort the file according to read mapping locations. For this we again are using samtools but this time the sort option. 

```
samtools sort our_mapped_reads.bam > our_mapped_reads.sort.bam
```
The output file our_mapped_reads.sort.bam is now a sorted and mapped read file that is necassary for subsequent analysis. Keep in mind that this file includes the same information as the previous files but sorted. 

You can see based on the file sizes that the compression and sorting significantly reduced the file size:
```
ll -h our_mapped_reads.*
```

This should show you that the sam file (our_mapped_reads.sam) is the largest file with a few hundred MB and the compressed and sorted bam file (our_mapped_reads.sort.bam) is actually the smalles file. 

Since these files contain all the same inforation we dont need to keep the larger files anymore. To remove them from your disk we run:

```
rm our_mapped_reads.bam 
rm our_mapped_reads.sam 
```

The last step that is necassary for a subsequent analysis is to index the sorted and compressed read file:
```
samtools index our_mapped_reads.sort.bam
```

Thus in the end you should have 2 files: our_mapped_reads.sort.bam and our_mapped_reads.sort.bam.bai . The latter is the index file. 

### Mapping QC: 

Before we move on to the variant analysis we want to inspect the mapped read file a little to see if all steps worked as expected. 

First we want to count the mapped and unmapped reads. This can be simply done using samtools view. For this we need to query /filter the reads based on their sam tag. You can here refresh what each tag stands for: https://broadinstitute.github.io/picard/explain-flags.html 

To compute the number of mapped reads we run:

```
 samtools view -c -F 4 our_mapped_reads.sort.bam
 ```
 The parameter -c tells samtools view to olny count and report that number to you. The parameter -F 4 tells it to only use reads that are in disagreement with the flag:4 . You can see based on the above URL that this flag represents unmapped reads. Thus we are querying not unmapped reads , which is the count of mapped reads. 
 
Often we want to restrict this given a certain mapping quality threshold. This can be done like this:

```
 samtools view -q 20 -c -F 4 our_mapped_reads.sort.bam
```
Here the addition of -q 20 restricts the reads to have mapping quality of 20 or more that is typically indicative of highly trustful mappings. 
Next we want to compute the number of unmapped reads:

```
 samtools view -c -f 4 our_mapped_reads.sort.bam
 ```
 Note we only need to change the -F to a -f to query for the Tag: 4. In our pre run this resulted in a mapping rate of 70.56%.
 
 Using samtools view we can also do other manipulations/ query. For example we can query the reads that are split reads, but only the supplement splits. These reads are goign to be useful later for the detection of SV. 
 
 ```
samtools view -c -f 2048  our_mapped_reads.sort.bam 
 ```
 Again feel free to check out what the tag 2048 means over at  https://broadinstitute.github.io/picard/explain-flags.html 
 
 
 Lastly we can also compute the reads that are mapped on the + vs. - strand. For some type of library preparation this is an important metric:
 ```
 samtools view -c -q 20 -f 16  our_mapped_reads.sort.bam
 samtools view -c -q 20 -f 0  our_mapped_reads.sort.bam
```
This time the -f 16 filters for reads on the - strand and the -f 0 for reads that mapped to the + strand. 


## Variant calling: 

Now that we have confidence in our mapped read file and we know its the right format and sorte we can continue with the variant calling. Firt we will call variants for SNV and subsequently for SV. 

### SNV calling: 


Now lets take some time to explore the mapped read file. 
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
