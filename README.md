# Ankara
Population and Adaptation Genomics Workshop, Turkey, 2016

## Material

Slides can be found [here](https://github.com/mfumagalli/Ankara/tree/master/Slides).

The data has been already downloaded and it is provided in `/truba/home/egitim/Data`.
These instructions, including all relevant files and scripts, can be found at `/truba/home/egitim/Ankara`.

## Data

As an illustration, we will use 80 BAM files of human samples (of African, European, East Asian, and Native American descent), a reference genome, and putative ancestral sequence.
*To make things more interesting, we have downsampled our data to an average mean depth of 2X!*

We will also use VCF files for 120 individuals from the same populations.
The human data represents a small genomic region (1MB on chromosome 2) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).

## Case study

*MOTIVATION*

Detecting signatures of natural selection in the genome has the twofold meaning of (i) understanding which adaptive processes shaped genetic variation and (ii) identifying putative functional variants.
In case of humans, biological pathways enriched with selection signatures include pigmentation, immune-system regulation and metabolic processes.
The latter may be related to human adaptation to different diet regimes, depending on local food availability (e.g. the case of lactase persistence in dairy-practicing populations).

The human Ectodysplasin A receptor gene, or EDAR, is part of the EDA signalling pathway which specifies prenatally the location, size and shape of ectodermal appendages (such as hair follicles, teeth and glands).
EDAR is a textbook example of positive selection in East Asians (Sabeti et al. 2007 Nature) with tested phenotypic effects (using transgenic mice).

Recently, a genome-wide association study found the same functional variant in EDAR associated to several human facial traits (ear shape, chin protusion, ...) in Native American populations (Adhikari et al. Nat Commun 2016).

*HYPOTHESIS*

- Is the functional allele in East Asian at high frequency in other human populations (e.g. Native Americans)?
- Can we identify signatures of natural selection on EDAR in Native Americans?
- Is selection targeting the same functional variant?

*CHALLENGES*
- Admixed population
- Low-depth sequencing data
- ...

*PLAN OF ACTION*

Goal day 1:

- Estimate allele frequencies for tested variant for African, European, East Asian and Native American samples from low-depth sequencing data

Goal day 2:

- Investigate population structure of American samples related to Europeans and Africans
- Estimate the site frequency spectrum for each population 

Goal day 3:

- Perfom a sliding windows scan based on allele frequency differentiation
- Assess statistical significance of selection signatures through simulations
- Test for extended haplotype homozygosity on high-depth sequencing data

Write the paper!

## Preparation

First create a directory where you will be working on.
```
mkdir Ex
cd Ex
```
Please set the path for all programs and data we will be using.
```
ANGSD=/truba/home/egitim/bin/angsd
NGSTOOLS=/truba/home/egitim/bin/ngsTools
NGSDIST=$NGSTOOLS/ngsDist
MS=/truba/home/egitim/bin/msHOT/ms
SS=/truba/home/egitim/bin/selscan-master/bin/linux
NGSADMIX=/truba/home/egitim/bin/NGSadmix
FASTME=/truba/home/egitim/bin/fastme-2.1.5/binaries/fastme-2.1.5-linux64
SAMTOOLS=/truba/home/egitim/bin/samtools-1.3.1/samtools
BGZIP=/truba/home/egitim/bin/htslib-1.3.1/bgzip
R=/truba/home/egitim/bin/R-3.2.2/bin/R
RSCRIPT=/truba/home/egitim/bin/R-3.2.2/bin/Rscript
```
You also need to provide the location of data and sequences:
```
DIR=/truba/home/egitim/Ankara
DATA=/truba/home/egitim/Data
REF=$DATA/ref.fa.gz
ANC=$DATA/anc.fa.gz
```
You can check that it worked by typing:
```
echo $REF
```
Finally, create a folder where you will put all the results and some temporary data.
```
mkdir Results
mkdir Data
```
Check if you are in the right location and all folders have been created:
```
pwd
ls
```

## Agenda

### Monday

#### Lecture

* Basics of data handling and filtering
* Maximum likelihood and Bayesian estimation
* Genotype likelihoods
* Allele frequencies, SNPs and genotypes calling
* Imputation

#### [Practical](Files/day1.md)

* Basic data filtering
* Estimation of allele frequencies and SNP calling
* Genotype calling
* Example: estimation of allele frequencies from low-depth sequencing data: the case of EDAR genetic variation in Native Americans

### Tuesday

* Population structure and site frequency spectrum
* Estimation summary statistics from low-depth sequencing data
* Brief notes on experimental design

#### [Practical](Files/day2.md)

* PCA from low-depth sequencing data
* SFS estimation and bias introduced by low-quality data
* Example: population structure and frequency variation in modern populations

### Wednesday

#### Lecture

* The effect of selection on the genome
* Methods to detect selection signals
* The problem of assessing significance
* Summary statistics from low-depth data

#### Practical

* Selection scan based on genetic [differentiation](Files/day3a.md) from low-depth data
* Assessing significance through [simulations](Files/day3b.md)
* Selection test based on [haplotype](Files/day3c.md) diversity
* Example: detection of natural selection from low-depth sequencing data and haplotype data: the case of EDAR genetic variation in Native Americans


## Credits

Some materials have been borrowed (and then adapted) from [Thorfinn Korneliussen](http://scholar.google.co.uk/citations?user=-YNWF4AAAAAJ&hl=en), [Anders Albrechtsen](http://popgen.dk/albrecht/web/WelcomePage.html), [Tyler Linderoth](http://scholar.google.com/citations?user=dTuxmzkAAAAJ&hl=en), [Filipe G. Vieira](http://scholar.google.com/citations?user=gvZmPNQAAAAJ&hl=en), [Dean Ousby](https://www.linkedin.com/in/deanousby), [Javier Mendoza](https://www.ucl.ac.uk/candela/candela-news/new-fellow-javiermendoza) and possibly many more. Many thanks to [Ryan Waples](http://www1.bio.ku.dk/english/Staff/?pure=en/persons/545443) for feedbacks.



