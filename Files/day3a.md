
Here you will learn how to perform a scan for selection by calculating PBS (population branch statistic) in windows from low-depth data.

As reference, these are the labelling for each population:

- LWK: Africans
- TSI: Europeans
- CHB: East Asians
- NAM: Native Americans

-----------------------------

### Allele frequency differentiation

The joint-SFS can be considered as a summary statistics for the level of genetic differentiation between populations.
We have seen how it can be used as prior information when computing FST (and related metrics) without relying on genotype calling.

Here we see how to compute the 2D-SFS, FST, PBS, and other summary statistics from low-depth data using ANGSD.
Our final goal is to detect signatures of selection in our data, with the specific example of EDAR gene in Native Americans.

To compute FST/PBS we first need to estimate the marginal and joint sites frequency spectra for our populations.
We have already seen how to estimate a single-population SFS.

---------------------------------------

Now, we need to estimate a **multi-dimensional SFS**, for instance the joint SFS between 2 populations (2D).
This can be used for making inferences on their divergence process (time, migration rate and so on).
However, here we are interested in estimating the 2D-SFS as prior information for our FST/PBS.

An important issue when doing this is to be sure that we are comparing the exactly same corresponding sites between populations.
ANGSD does that automatically and considers only a set of overlapping sites.

We are performing PBS assuming NAM (Native Americans) being the targeted population.
The 2D-SFS between all populations and NAM are computed with:
```
#!/bin/sh

# specify where the program is
ANGSD=/truba/home/egitim/bin/angsd

POP2=NAM
for POP in LWK TSI CHB
do
        echo $POP
        $ANGSD/misc/realSFS Results/$POP.saf.idx Results/$POP2.saf.idx 2> /dev/null > Results/$POP.$POP2.sfs
done
# we also need the comparison between LWK and TSI as we will see later 
$ANGSD/misc/realSFS Results/LWK.saf.idx Results/TSI.saf.idx 2> /dev/null > Results/LWK.TSI.sfs
```

The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] and so on.
```
less -S Results/LWK.NAM.sfs
```
You can plot it, but you need to define how many samples you have per population.
```
Rscript $DIR/Scripts/plot2DSFS.R Results/LWK.NAM.sfs 20 20
```
Transfer it to your local machine:
```
open Results/LWK.NAM.sfs.pdf
```

You can even estimate SFS with higher order of magnitude.
This command may take some time and you should skip it if not interested.
```
# $ANGSD/misc/realSFS Results/LWK.saf.idx Results/TSI.saf.idx Results/NAM.saf.idx 2> /dev/null > Results/LWK.TSI.NAM.sfs
```

------------------------------------

Here we are going to calculate **allele frequency differentiation** using the PBS (population branch statistic) metric.
Again, we can achieve this by avoid genotype calling using ANGSD.
From the sample allele frequencies likelihoods (.saf files) we can estimate PBS using the following pipeline.

Note that here we use the previously calculated SFS as prior information.
Also, PEL is our target population, while CHB and TSI are reference populations.
If not already done, you should calculate .saf.idx files for each population, as explained in the section above.

Therefore, we need to use the 2D-SFS between TSI and CHB and NAM (already done).

If we also assume CHB being the target population, as a possible separate analysis is (you don't have to run this if not interested):
```
#!/bin/sh

# specify where the program is
ANGSD=/truba/home/egitim/bin/angsd

POP2=CHB
for POP in LWK TSI
do
        echo $POP
        $ANGSD/misc/realSFS Results/$POP.saf.idx Results/$POP2.saf.idx 2> /dev/null > Results/$POP.$POP2.sfs
done
```

The 2D-SFS will be used as prior information for the joint allele frequency probabilities at each site.
From these probabilities we will calculate the population branch statistic (PBS) using the NAM (and/or CHB) as target population and LWK and TSI as reference populations.
Our goal is to detect selection in NAM (and/or CHB) in terms of allele frequency differentiation.

Specifically, we are computing a slinding windows scan, with windows of 50kbp and a step of 10kbp.
This can be achieved using the following commands.

1) This command will compute per-site FST indexes (please note the order of files):
```
#!/bin/sh

# specify where the program is
ANGSD=/truba/home/egitim/bin/angsd

# NAM
$ANGSD/misc/realSFS fst index Results/LWK.saf.idx Results/TSI.saf.idx Results/NAM.saf.idx -sfs Results/LWK.TSI.sfs -sfs Results/LWK.NAM.sfs -sfs Results/TSI.NAM.sfs -fstout Results/NAM.pbs &> /dev/null
# CHB
$ANGSD/misc/realSFS fst index Results/LWK.saf.idx Results/TSI.saf.idx Results/CHB.saf.idx -sfs Results/LWK.TSI.sfs -sfs Results/LWK.CHB.sfs -sfs Results/TSI.CHB.sfs -fstout Results/CHB.pbs &> /dev/null
```
and you can have a look at their values:
```
# NAM
$ANGSD/misc/realSFS fst print Results/NAM.pbs.fst.idx | less -S
# CHB
$ANGSD/misc/realSFS fst print Results/CHB.pbs.fst.idx | less -S
```
where columns are: chromosome, position, (a), (a+b) values for the three FST comparisons, where FST is defined as a/(a+b).
Note that FST on multiple SNPs is calculated as sum(a)/sum(a+b).

2) The next command will perform a sliding-window analysis:
```
#!/bin/sh

# specify where the program is
ANGSD=/truba/home/egitim/bin/angsd

# NAM
$ANGSD/misc/realSFS fst stats2 Results/NAM.pbs.fst.idx -win 50000 -step 10000 > Results/NAM.pbs.txt 2> /dev/null
# CHB
$ANGSD/misc/realSFS fst stats2 Results/CHB.pbs.fst.idx -win 50000 -step 10000 > Results/CHB.pbs.txt 2> /dev/null
```

Have a look at the output file:
```
# NAM
less -S Results/NAM.pbs.txt
# CHB
less -S Results/CHB.pbs.txt
```
The header is:
```
region  chr     midPos  Nsites  Fst01   Fst02   Fst12   PBS0    PBS1    PBS2
```
Where are interested in the column `PB2` which gives the PBS values assuming our population (coded here as 2) being the target population.
Note that negative PBS and FST values are equivalent to 0.

We are also provided with the individual FST values.
You can see that high values of PBS2 are indeed associated with high values of both Fst02 and Fst12 but not Fst01.
We can plot the results along with the gene annotation.
```
# NAM
Rscript $DIR/Scripts/plotPBS.R Results/NAM.pbs.txt Results/NAM.pbs.pdf
# CHB
Rscript $DIR/Scripts/plotPBS.R Results/CHB.pbs.txt Results/CHB.pbs.pdf
```

It will also print out the maximum PBS value observed as this value will be used in the next part.
This script will also plot the PBS variation in LWK as a control comparison.
```
# NAM
open Results/NAM.pbs.pdf
# CHB
open Results/CHB.pbs.pdf
```
Comment the results.

-------------------------

**OPTIONAL**

We are also interested in assessing whether an increase in allele frequency differentiation is also associated with a change of **nucleotide diversity** in CHB.
Again, we can achieve this using ANGSD by estimating levels of diversity without relying on called genotypes.

The procedure is similar to what done for PBS, and the SFS is again used as a prior to compute allele frequencies probabilities.
From these quantities, expectations of various diversity indexes are compute.
This can be achieved using the following pipeline.

First we compute the allele frequency posterior probabilities and associated statistics (-doThetas) using the SFS as prior information (-pest)
```
#!/bin/sh

# specify where the program is
ANGSD=/truba/home/egitim/bin/angsd

# specify where the data is
DATA=/truba/home/egitim/Data
REF=$DATA/ref.fa.gz
ANC=$DATA/anc.fa.gz

# specify the population label
POP=CHB

$ANGSD/angsd -b $DATA/$POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 100 -doCounts 1 \
	-GL 1 -doSaf 1 \
	-doThetas 1 -pest Results/$POP.sfs &> /dev/null
```
Then we need to index thess file and perform a sliding windows analysis using a window length of 50kbp and a step size of 10kbp.
```
#!/bin/sh

# specify where the program is
ANGSD=/truba/home/egitim/bin/angsd

# specify where the data is
DATA=/truba/home/egitim/Data
REF=$DATA/ref.fa.gz
ANC=$DATA/anc.fa.gz

POP=CHB

# index files
$ANGSD/misc/thetaStat make_bed Results/$POP.thetas.gz &> /dev/null
# perform a sliding-window analysis
$ANGSD/misc/thetaStat do_stat Results/$POP.thetas.gz -nChr 1 -win 50000 -step 10000 -outnames Results/$POP.thetas &> /dev/null
```

Look at the results:
```
less -S Results/CHB.thetas.pestPG
```
and plot the sliding windows scan for nucleotide diversity:
```
Rscript $DIR/Scripts/plotSS.R
```
```
open Results/CHB.ss.pdf
```

------------------------

[HOME](https://github.com/mfumagalli/Ankara)



