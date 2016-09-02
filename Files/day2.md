

#### Population structure

**IMPORTANT NOTE**: These commands are given as a mere example as, in practise, such analyses should be performed on larger genomic regions.

We now want to investigate population structure of our sample. 
We perform a principal component analyses (PCA) without relying on called genotypes, but rather by taking their uncertainty into account.
More specifically, the next program we are going to use (ngsTools) takes as input genotype probabilities in binary format, so we need to specify `-doGeno 32` in ANGSD.
Also, we are using a HWE-based prior with `-doPost 1`.

```
$ANGSD/angsd -b $DATA/ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 70 -setMaxDepth 200 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-2 \
        -doGeno 32 -doPost 1 &> /dev/null
```
Unzip the results but you cannot open it since it is in binary format:
```
gunzip Results/ALL.geno.gz
```
We are going to use `ngsCovar`, which estimates the covariance matrix between individuals based on genotype probabilities.
Then this matrix will be decomposed into principal componenets which will be investigated for population structure analyses.
If you type:
```
$NGSTOOLS/ngsPopGen/ngsCovar
```
you will see a list of possible options.

For instance, we need to define how many sites we have.
To retrieve such value, we can inspect the file with allele frequencies:
```
less -S Results/ALL.mafs.gz
```
How many sites do we have?
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```
Now we can perform a PCA by estimating the covariance matrix using:
```
$NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/ALL.geno -outfile Results/ALL.covar -nind 80 -nsites $N_SITES -call 0 -norm 0 &> /dev/null
```
with the options `-call 0` meaning that we do not perform genotype calling and `-norm 0` that we are not normalising by allele frequency.
The latter may give more weight to low frequency variants which are harder to estimate.

Look at the output file:
```
less -S Results/ALL.covar
```
which represents a matrix of NxN with N individuals giving the covariance.
Note that this matrix is symmetric.

Finally, we perform an eigenvector decomposition and plot the resulting map:
```
# create a cluster-like file defining the labelling (population) for each sample
Rscript -e 'write.table(cbind(rep(seq(1,20),4),rep(seq(1,20),4),c(rep("CHB",20),rep("LWK",20),rep("NAM",20),rep("TSI",20))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="Results/ALL.clst", quote=F)'
# run and plot
Rscript $DIR/Scripts/plotPCA_ngstools.R Results/ALL.covar Results/ALL.clst 1-2 Results/ALL.pca.pdf
```
where the parameter `1-2` specifies that we are plotting only the first and second component.
On the screen, you will see a series of numbers.
These are the percentage of explained variance for each component.

Finally, you can open the produced image:
```
evince Results/ALL.pca.pdf
```

You can inspect other components:
```
Rscript $DIR/Scripts/plotPCA_ngstools.R Results/ALL.covar Results/ALL.clst 2-3 Results/ALL.pca2.pdf
evince Results/ALL.pca2.pdf
```

Therefore, NAM samples appear very close to EUR but separated.
Indeed, among all Latin American populations present in the 1000G project, Peruvians are the least admixed population.
We can either use all these samples or compute admixture proportions in order to select a subset of putative Native American (unadmixed) samples.
However with such limited data set we cannot refine the structure between NAM and CHB.

-----------------------------------

**OPTIONAL**

#### Admixture

**IMPORTANT NOTE**: These commands are given as a mere example as, in practise, such analyses should be performed on larger genomic regions.

We want to select only a subset of PEL samples with high Native American ancestry, or check the level of admixture in our NAM samples.
We use ngsAdmix, which again works on genotype probabilities and not on individual calls.
This is suitable for low-depth data.

ngsAdmix requires genotype likelihoods in BEAGLE format as input.
We can compute these quantities with ANGSD with `-doGlf 2`.
```
$ANGSD/angsd -b $DATA/ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 50 -setMaxDepth 200 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-4 \
        -doGlf 2 &> /dev/null
```

We assume 4 ancestral populations making up the genetic diversity of our samples.
Therefore we compute admixture proportions with 4 ancestral components.
```
K=4
$NGSADMIX -likes Results/ALL.beagle.gz -K $K -outfiles Results/ALL.admix.K$K -P 4 -minMaf 0.02 &> /dev/null
```

We now combine samples IDs with admixture proportions and inspect the results.
```
paste $DATA/ALL.bamlist Results/ALL.admix.K$K.qopt > Results/ALL.admix.K$K.txt
less -S Results/ALL.admix.K$K.txt
```

From these quantities we can extract how many samples (and which ones) have a high proportion of Native American ancestry (e.g. >0.90). 
------------------

**OPTIONAL**

#### Genetic distances

**IMPORTANT NOTE**: These commands are given as a mere example as, in practise, such analyses should be performed on larger genomic regions.

Genotype probabilities can be used also to infer the structure of your population.
For instance, in our example, they can used to assess whether PEL samples are indeed admixed.

We can compute genetic distances as a basis for population clustering driectly from genotype probabilities, and not from assigned genotypes as we have seen how problematic these latters can be at low-depth.

First, we compute genotype posterior probabilities jointly for all samples:
```
# Assuming HWE, without filtering based on probabilities, with SNP calling
$ANGSD/angsd -b $DATA/ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 50 -setMaxDepth 200 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null
```

Next we record how many sites we retrieve.
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```

Then we create a file with labels indicating the population of interest for each sample.
```
Rscript -e 'cat(paste(rep(c("LWK","TSI","CHB","NAM"),each=20), rep(1:20, 4), sep="_"), sep="\n", file="Data/pops.label")'
cat Data/pops.label
```

With [ngsDist](https://github.com/fgvieira/ngsDist) we can compute pairwise genetic distances without relying on individual genotype calls.
```
$NGSDIST/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 80 -n_sites $N_SITES -labels Data/pops.label -o Results/ALL.dist -n_threads 4 &> /dev/null
less -S Results/ALL.dist
```

We can visualise the pairwise genetic distances in form of a tree.
```
$FASTME -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null
cat Results/ALL.tree
```

Finally, we plot the tree.
```
Rscript -e 'library(ape); library(phangorn); pdf(file="Results/ALL.tree.pdf"); plot(read.tree("Results/ALL.tree"), cex=0.5); dev.off();' &> /dev/null
evince Results/ALL.tree.pdf
```

One can also perform a PCA/MDS from such genetic distances to further explore the population structure.

[HOME](https://github.com/mfumagalli/Ankara)










