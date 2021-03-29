
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
  * Common Phylogenetic Tree Building Pitfalls


## Background On phylogenetic analysis

Phylogenetics investigates the species evolutionary processes at different biological scales ranging from within-host evolution of persistently infecting viruses, such as HIV and HBV, to the global epidemic spread of SARS-CoV-2 for instance. To this purpose, we aim at integrating molecular biology and computational biology models to understand the underlying mechanisms of disease in populations. 

***

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

***

## Commands 

* Run mafft to align genetic dataset using a reference sequence
```
user$ mafft --auto --thread 4 --addfragments Study_and_Background.fas Wuhan-Hu-1.fasta > Study_and_Background_al.fas
```

* Run fastatool.py to remove duplicate sequences
```
user$ python3 fastatool.py 
---Starting Fasta Modifcation Tool---

What operation (subsample,clean,tag,extension,review,genome): clean
*** Starting Fasta Clean ***
** Remove duplicate sequences, remove sparse columns, and remove sparse sequences **
Filepath for input (must be .fa (or .fas) or .csv format): Study_and_Background_aled.fas
INITIAL COUNT: 1848
Percent under max length for which columns will be cut (enter .5 for 50%): .5
REMOVED COLUMNS COUNT: 0
Percent under max length for which lines will be cut (enter .75 for 75%): .75
REMOVED LENGTH COUNT: 0
Example tag values:

0  >DC                             >hCoV-19/USA/TX-DSHS-2612/2021
1  MT646069                        EPI_ISL_983609                
2  B.1.1.33                        B.1.1.33                      
3  USA                             USA                           
4  2020-03-07                      2021-01-23                    
5  ----------------AGATCTGTTCT...  ----------------AGATCTGTTCT...

Enter the tag numbers for which duplicates will be compared by (starts at 0, comma-separated): 3,4,5
REMOVED DUP COUNT: 125

Filepath for output (will create a new file if none exist): Study_and_Background_aled_woDups.fas
```

* Run subsample_covid.py to subsample large genetic datasets based on metadata
```
user$ python3 subsample_covid.py 

Enter FASTA-format file for input: Study_and_Background_aled_noDups.fas
Enter name for FASTA-format output file: Study_and_Background_aled_noDups_sub.fas
Enter delimiter, or return for '|': 
Example tag values:

0  DC                              hCoV-19/Uruguay/UY-NYUMC857...
1  MT646069                        EPI_ISL_457953                
2  B.1.1.33                        B.1.1.33                      
3  USA                             Uruguay                       
4  2020-03-07                      2020-03-23                    

Note: only the year will be used if the date field is included

Enter one or more field numbers separated by commas or whitespace, or S to load spreadsheet: 3,4

35 categories
Mean category size: 49.23
Minimum: 1
Maximum: 1053

Save a spreadsheet showing category sizes? y(es)/N(o) (RETURN for No): y
Enter name for output spreadsheet (CSV) file: Study_and_Background_aled_noDups.csv

Wrote CSV file Study_and_Background_aled_noDups.csv

Launch CSV file Study_and_Background_aled_noDups.csv now? y(es)/N(o) (RETURN for No): y

Launched CSV file

Total samples for some choices of samples per category:

   1  35                    60  401
   2  59                    70  435
   3  75                    80  465
   4  91                    90  495
   5  104                  100  517
   6  117                  200  717
   7  129                  300  917
   8  141                  400  1070
   9  151                  500  1170
  10  161                  600  1270
  20  228                  700  1370
  30  278                  800  1470
  40  321                  900  1570
  50  361                 1000  1670

Number to sample from each category: 30

278 sequences will be sampled

Save a spreadsheet showing category sizes and number sampled for editing and reloading? y(es)/N(o) (RETURN for No): y
Enter name for output spreadsheet (CSV) file: Study_and_Background_aled_noDups_Sub.csv

Wrote CSV file Study_and_Background_aled_noDups_Sub.csv

Launch CSV file Study_and_Background_aled_noDups_Sub.csv now? y(es)/N(o) (RETURN for No): y

Launched CSV file

Specify the criteria for choosing within a category, 
in order of importance.  The first criterion takes priority, 
so the others only matter when candidates are "tied" for the
first.  Similarly, the second takes priority over all but 
the first, and so on.  Randomization applies in all cases.

Specify one or more letters (either case), separated by whitespace, commas, or nothing:

N  No preferences.  Purely random choice.  Must be only letter specified.
U  Seek uniformity of date distribution.
   If used, this must come first in priorities,
   and the date (in effect, year) must be part of
   the category definition.
L  Maximize sequence length (includes internal gaps, but not leading and trailing gaps)
D  Completeness of date.  YMD > YM > Y
M  Presence of month; no preference for YMD over YM.
   It can be meaningful and useful to provide both D and M,
   provided that M comes first and something comes between, e.g., MLD.
G  Minimize number of internal gaps with lengths not divisible by 3 and less than 5
A  Minimize number of ambiguity characters (N, etc.).

Specify preferences for isolate selection in order of priority: UDGAL

Chose 30 (or all) from each category using fields [3, 4] and preferences UDGAL
278 sequences written to Study_and_Background_aled_noDups_sub.fas
```

Next: [done!](../README.md)
