
#  Module #2: Read mapping and variant calling

***
Lead: Fritz J Sedlazeck

TA: ‪Medhat Mahmoud + Sairam Behera

Baylor College of Medicine
***


## Goals of this module
The goal of this module is to get you familiarized with the mapping of Illumina reads to the [SARS-CoV-2 genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta) and to identify variations. For the latter we will focus on SNV (point mutations) and Structural Variations (SVs). In addition to identifying these two types of variations we will further show you one way to assess the quality and summary statistics across the mapping, SNVs and SVs calling.

### The main steps in this Module are:
1. Align the reads from Module2/raw_reads ([BWA-mem](https://github.com/lh3/bwa))
2. Obtain alignment statistics ([samtools](https://github.com/samtools/samtools))
3. Identifying single-nucleotide variations (SNPs) ([LoFreq](https://github.com/andreas-wilm/lofreq3))
5. Identify Structural Variations (SVs) ([Manta](https://github.com/Illumina/manta))


*Note*: The tools and methods we introduce here are clear snapshots of best practices for this analysis and there are multiple other methods available for your own analysis.


## Learning Objectives
In the end of Module 2 you should have a clear understanding how to map Illumina reads and identify SNVs and SVs together. Furthermore, you should be able to identify problems and obtain preliminary insights into your data set at hand.

We will further discuss standard file formats for this application: [BAM](http://samtools.github.io/hts-specs/SAMv1.pdf) and [VCF](http://samtools.github.io/hts-specs/VCFv4.1.pdf) files.

***

## Alignment and quality control (QC) of aligning short-read data

### Mapping reads

As discussed in the lecture the purpose of a mapping / alignment is to identify the most likely region a given read (i.e. sequenced segment) was originating from a given genome location.

First we will align the reads to the reference genome.

Navigate to the folder and create a folder for your mapping results

```
cd Module2
mkdir -p mapping
cd mapping
```

Next we want to use the reads to start the [BWA](https://github.com/lh3/bwa) mem alignments. `BWA mem` is currently one of the standard short-read based mapper and that is why we are using it here. To see the available options of `bwa` just execute its command (or see [manual](http://bio-bwa.sourceforge.net/bwa.shtml)):

```
bwa
```

#### 1. Index reference genome:
We first need to index the reference genome itself by using `bwa index`.

```
bwa index reference.fasta
```
This should only take a few seconds since the SARS-Cov-2 genome is tiny.

#### 2. Map reads to genome:
Next we are ready to start mapping the fastq reads to the genome itself. For this we want to use the `bwa mem` option that is best capable to handle the Illumina paired end reads.

```
 bwa mem -t 2 reference.fasta SRR12447392_1.fastq SRR12447392_2.fastq > our_mapped_reads.sam
```

This executes bwa mem with `2` threads (`-t` parameter) give our previously indexed `reference.fasta` and the two fastq files representing the Illumina paired end reads.

After a few seconds the program ends and we have our first result as : `our_mapped_reads.sam`. This is a standard text file (see SAM file format specification [here](https://samtools.github.io/hts-specs/SAMv1.pdf)) and we can take a look. As highlighted in the lecture we have a header in this file indicated with `@` and then entries per read per line.

#### 3. Converting a SAM file to a BAM file

For subsequent analysis we need to compress (SAM -> BAM) the file. For this we are using [samtools](https://github.com/samtools/samtools) with the option: `view`

```
samtools view -hb our_mapped_reads.sam > our_mapped_reads.bam
```

The options `-h` ensures that the header is kept for the output file and the option `-b` tells `samtools` that we want to obtain the compressed (BAM) version.


#### 4. Sorting a BAM file

Next we need to sort the file according to read mapping locations. For this we again are using `samtools` but this time the `sort` option.

```
samtools sort our_mapped_reads.bam > our_mapped_reads.sort.bam
```
The output file `our_mapped_reads.sort.bam` is now a sorted and mapped read file that is necessary for subsequent analysis. Keep in mind that this file includes the same information as the previous files but sorted.

You can see based on the file sizes that the compression and sorting significantly reduced the file size:
```
ll -h our_mapped_reads.*
```

This should show you that the sam file (`our_mapped_reads.sam`) is the largest file with a few hundred MB and the compressed and sorted bam file (`our_mapped_reads.sort.bam`) is actually the smallest file.

Since these files contain all the same information we don't need to keep the larger files anymore. To remove them from your disk we run:

```
rm our_mapped_reads.bam
rm our_mapped_reads.sam
```

*bonus*    
You can get sorted bam direct output from `bwa mem` as follow:
```
bwa mem -t 2 reference.fasta SRR12447392_1.fastq SRR12447392_2.fastq | samtools sort -@2 -o our_mapped_reads.sort.bam -
```
Where  `-@2` means use two threads, and `-` means take input from stdout  


#### 5. Creating a BAM index file

The last step that is necessary for a subsequent analysis is to index the sorted and compressed read file:
```
samtools index our_mapped_reads.sort.bam
```

Thus in the end you should have 2 files: `our_mapped_reads.sort.bam` and `our_mapped_reads.sort.bam.bai`. The latter is the index file.

### Mapping QC

Before we move on to the variant analysis we want to inspect the mapped read file a little to see if all steps worked as expected.

First we want to count the mapped and unmapped reads. This can be simply done using `samtools view`. For this we need to query / filter the reads based on their sam tag. You can here refresh what each tag stands for: https://broadinstitute.github.io/picard/explain-flags.html

To compute the number of mapped reads we run:

```
samtools view -c -F 4 our_mapped_reads.sort.bam
```
 The parameter `-c` tells samtools view, to only count and report that number to you. The parameter `-F 4` tells it to only use reads that are in disagreement with the flag:4 . You can see based on the above [URL](https://broadinstitute.github.io/picard/explain-flags.html) that this flag represents unmapped reads. Thus we are querying not unmapped reads, which is the count of mapped reads.

Often we want to restrict this given a certain mapping quality threshold. This can be done like this:

```
 samtools view -q 20 -c -F 4 our_mapped_reads.sort.bam
```
Here the addition of `-q 20` restricts the reads to have mapping quality of 20 or more that is typically indicative of highly trustful mappings.
Next we want to compute the number of unmapped reads:

```
 samtools view -c -f 4 our_mapped_reads.sort.bam
 ```
 Note we only need to change the `-F` to a `-f` to query for the Tag: 4. In our pre run this resulted in a mapping rate of 70.56%.

 Using samtools view we can also do other manipulations/queries. For example we can query the reads that are split reads, but only the supplement splits. These reads are going to be useful later for the detection of SV.

 ```
samtools view -c -f 2048  our_mapped_reads.sort.bam
 ```
 Again feel free to check out what the tag 2048 means over at  https://broadinstitute.github.io/picard/explain-flags.html


 Lastly we can also compute the reads that are mapped on the `+` vs. `-` strand. For some type of library preparation this is an important metric:
 ```
 samtools view -c -q 20 -f 16  our_mapped_reads.sort.bam
 samtools view -c -q 20 -f 0  our_mapped_reads.sort.bam
```
This time the `-f 16` filters for reads on the `-` strand and the `-f 0` for reads that mapped to the `+` strand.

***

## Variant calling

Now that we have confidence in our mapped reads file and we know that it's the right format and reads are sorted, we can continue with the variant calling. First, we will call variants for SNV and subsequently for SV.  

To keep everything nice and tight we will change directory and create one called SNV:
```
cd ..
mkdir SNV
cd SNV
```

### SNV calling
For SNV calling we are going to use [LoFreq](https://github.com/andreas-wilm/lofreq3), which was first published in 2014: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526318/

Given our mapped read file and our reference fasta file we can execute Lofreq as follow:

```
lofreq call  -f reference.fasta -o our_snv.vcf --min-mq 10 our_mapped_reads.sort.bam
```
Overall this step will run for a couple of minutes so feel free to drink something:tea::coffee:  or stretch:walking::running:! :smile:

This will first start to index our `reference.fasta` and subsequently use our mapped reads to call SNVs. Note we have specified a mapping quality of minimum 10 (`--min-mq 10`).

In the end the program `lofreq` has produced a VCF file as its output: `our_snv.vcf`. We can open this file like:

```
less -S our_snv.vcf
```
Note to terminate this process press `q` to close less.
As we can see the VCF file follows a certain standard as it has first specified meta information as part of the header (`#`). This is then followed by each line showing a single variant.

Take your time to look into this file. Some of the important tags that are defined are `AF` (allele frequency within the sample), `DP4` list of supporting reads for reference and alternative split up over `+/-` strand. What is important to note is that each of these tags have to be defined in the header. Go and look up: `DP` and compare it to `DP4`.

To get a feeling about our file we want to query it a little to summarize our SNV calls.
First we want to count the total number of SNV in this file:
```
grep -vc '#' our_snv.vcf
```
This will count the number of lines that don't have an `#` in it. `-v` is inverting the match and `-c` is counting the number of these matches.

If we want to know if there is an imbalance in the nucleotides that has been changed we could use something simple like this:
```
grep -v '#' our_snv.vcf | cut -f 5 |sort | uniq -c
```
Again we are selecting against the header (`-v '#'`) then extracting column 5 (the alternative nucleotide) and sorting and counting the occurrence of each nucleotide (uniq ) with the `-c` option to count. You should see a clear preference for an `T` and `A` nucleotide that has been inserted.

We can also roughly and quickly see if there are hotspots for SNV along the genome:
```
grep -v '#' our_snv.vcf | cut -f 2 | awk '{print int($1/100)*100}'  | sort | uniq -c  | awk '$1 > 5 {print $0 }' | sort -n -k 2,2 | less
```
Here we extract similar to before the 2nd column (SNV position) and bin it by 100bp. Next we sort and count the occurrences and filter to have only regions that have more than 5 SNV within their 100bp. Lastly we make sure that the positions of the bins are sorted in the output.

A set of very useful methods are [bcftools](http://samtools.github.io/bcftools/) and [vcftools](https://vcftools.github.io/man_latest.html) to further filter and manipulate these files.

### SV calling
In the end we want to also identify Structural Variations (SVs). Here we are simply using [Manta](https://github.com/Illumina/manta), which was mainly designed to identify SV across a human genome.

Manta requires two steps:

#### 1. Initiate the run:
```
configManta.py --bam=our_mapped_reads.sort.bam --referenceFasta=reference.fasta --runDir=Out_Manta
```
This should just take seconds as it initiates the folder structure and specifies for the subsequent process to use our mapped reads and our reference file. In addition, we specify the output to be written in `Out_Manta`

#### 2. Run the analysis:
```
python Out_Manta/runWorkflow.py -j 2 -m local -g 10
```

This will launch the Manta pipeline that we previously configured. `-j` specifies the number of CPU threads (2 in our case), `-m` local indicates that it should not try to run things on different nodes or instances and `-g 30` specifies the available memory for the process in GB.

Manta now searches for abnormal paired-end reads and split-reads across our mapped reads. These will be analyzed together and clustered to identify SV in these samples. After ~2-3 minutes you should see that the program has finished.

Our SV calling results can be found here:
```
ls Out_Manta/results/variants
```
As you can see we have multiple VCF files. These represent the different stages of Manta and the confidence level for the SV calls. We typically use the `diploidSV.vcf.gz` file.




#PREVIOUSLY:


Now let's take some time to explore the mapped read file.
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

Create environment:  
`conda create -n SVanalysis python=2.7 bwa=0.7.17 manta=1.6.0 lofreq=2.1.5`

Activate the environment:  
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

<!-- Access the data: FIXME
```
cd awesome_data

```
This dataset contains paired end reads. -->

<!-- Next: [module3!](module3.rst) -->
