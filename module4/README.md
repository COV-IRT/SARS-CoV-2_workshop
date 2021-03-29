
#  Module #4: Phylogenetic analysis

***
Lead: Nídia S. Trovão

TA: James R. Otieno

Fogarty International Center, National Institutes of Health
***

## Goals of this module
* This hands on tutorial will teach you how to investigate the evolutionary history of SARS-CoV-2 genomes, through bioinformatics and phylogenetic analysis, and interpret the results to inform public health interventions.

## Learning Objectives
* Building a background dataset by retrieving sequences from public databases
  * Alignment viewer and editing software
  * Identify the study sequences’ lineage
  * Access sequences on public genetic databases
  * Assembly of the background dataset
* Bioinformatics towards phylogenetics
  * Selection of a reference genome
  * Alignment software
  * Inspection and manual curation of the alignment
  * Appropriate subsample of the background dataset for computational efficiency
* Phylogenetic analysis
  * Components of a phylogeny
  * Inference algorithms and software
  * Phylogenetic signal consideration for advances phylodynamic modeling
  * Understanding the Bootstrap
  * Visualizing Phylogenetic Trees in FigTree
  * Basic Interpretation of Phylogenetic Trees

***

## Background On phylogenetic analysis

Phylogenetics investigates the species evolutionary processes at different biological scales ranging from within-host evolution of persistently infecting viruses, such as HIV and HBV, to the global epidemic spread of SARS-CoV-2 for instance. To this purpose, we aim at integrating molecular biology and computational biology models to understand the underlying mechanisms of disease in populations. 

## Slide Presentations

[[Download]](https://github.com/COV-IRT/SARS-CoV-2_workshop/tree/main/module4/Slide%20Presentations)

## Supporting Materials

[[Download]](https://github.com/COV-IRT/SARS-CoV-2_workshop/tree/main/module4/Supporting%20Materials)

## Software required

[[Download]](https://drive.google.com/drive/folders/1imHg2WHql35rGzXpOyGZ6gKtUGYdcvNu?usp=sharing)

<!--
**Windows software (download and install all)**  
[[Download]](https://www.oracle.com/java/technologies/javase-jre8-downloads.html) Java JRE x64 v8 _pick the 1st "Windows x64" on the list_
[[Download]](https://ormbunkar.se/aliview/downloads/windows/windows-version-1.26/AliView-Setup.exe) AliView v1.26
[[Download]](https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree.v1.4.4.zip) FigTree v1.4.4
[[Download]](https://www.megasoftware.net/) MEGA v10 _pick Windows|Graphical(GUI)|MEGA X(64bit)_
[[Download]](https://github.com/beast-dev/beast-mcmc/releases/download/v1.5.3-tempest/TempEst.v1.5.3.zip) TempEst v1.5.3
[[Download]](/python.html) Python v3.8.6 and Biopython
[[Download]](https://github.com/notepad-plus-plus/notepad-plus-plus/releases/download/v7.9/npp.7.9.Installer.x64.exe) Notepad Plus Plus v7.9 (Text editor)
**macOS software (download and install all)**  
[[Download]](https://www.oracle.com/java/technologies/javase-jre8-downloads.html) Java JRE x64 v8 _pick "macOS x64 Installer"_
[[Download]](https://ormbunkar.se/aliview/downloads/mac/AliView-1.26-app.zip) AliView v1.26
[[Download]](https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree.v1.4.4.dmg) FigTree v1.4.4
[[Download]](https://www.megasoftware.net/) MEGA v10 _pick macOS|Graphical(GUI)|MEGA X(64bit)_
[[Download]](https://github.com/beast-dev/beast-mcmc/releases/download/v1.5.3-tempest/TempEst.v1.5.3.dmg) TempEst v1.5.3
[[Download]](/python.html) Python v3.8.6 and Biopython
[[Download]](https://s3.amazonaws.com/BBSW-download/BBEdit_13.5.dmg) BBEdit v13.5 (Text editor)
-->

**Web tools (installation not required)**  

[[Access]](https://pangolin.cog-uk.io/) PANGOLIN server web

[[Access]](https://www.epicov.org/epi3/start) GISAID website web

[[Access]](https://mafft.cbrc.jp/alignment/server/add_fragments.html) MAFFT server web

[[Access]](http://iqtree.cibiv.univie.ac.at/) IQTree server 1 web

[[Access]](https://www.hiv.lanl.gov/content/sequence/IQTREE/iqtree.html) IQTree server 2 web


<!--
The first step to starting phylogenetic analysis -> context to DNAnexus
-->
<!--
*Connect to your DNAnexus instance and open up a shell prompt.*
-->
<!--
Install XYZ
```
conda install -c bioconda -y XYZ
```
-->
<!--
Accesws the data: FIXME
```
cd awesome_data

```
This dataset contains paired end reads. 
-->
Next: [done!](../README.md)
