
## Test of genotype likelihoods

```
ssh life

ssh balloux

cd ~/Downloads

/data/data/Software/julia/julia 

include("/data//Software/ngsPoly/functions.jl")

mySite = Site("seq1", 272, 'A');

# """!"#\$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~""";

ascii2phred('5')-1

10^((33 - Int64('5')  )/10)

myReads = Reads("AAAT" , "5555");
```

```
myReads = Reads("AAAT" , "5555");
calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 4)
```

- AA:  -5.73393
- AT:  -2.79934
- TT: -17.1214

```
myReads = Reads("AAATT" , "55555");
calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 4)
```

- AA:  -11.4377 
- AT:  -3.49918
- TT: -17.1314

```
10^((33 - Int64('+')  )/10)
myReads = Reads("AAAT" , "555+");
calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 4)
```

- AA:  -3.43135
- AT:  -2.86165
- TT: -17.2167 

```
myReads = Reads("AAATAAC" , "555-43(");
calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 4)
```

- AA:  -6.63093
- AT:  -6.9514 
- TT: -30.6037 

--------------------

Allele frequencies

```
myReads = Reads("AAATAAC" , "555-43(");

haploid=calcGenoLogLike1(myReads, mySite)

(major, minor, minor2, minor3) = sortperm(haploid, rev=true)

freqsMLE = optimFreq_MajorMinor_GSS(myReads, major, minor, 1e-5)

maf = 0.0;
like_H0 = calcFreqLogLike1_MajorMinor(myReads, major, minor, maf)

like_H1 = freqsMLE[1]

lrtSnp = snpPval_MajorMinor(myReads, freqsMLE[1], major, minor)

```

- H0: -6.630932301040746
- H1: -5.420818481979177
- LRT: 2.420227638123137
- f_MLE: 0.14818398119374376


```
myReads = Reads("AAATAACAATT" , "555-43(8822");

haploid=calcGenoLogLike1(myReads, mySite)

(major, minor, minor2, minor3) = sortperm(haploid, rev=true)

freqsMLE = optimFreq_MajorMinor_GSS(myReads, major, minor, 1e-5)
maf = 0.0;
like_H0 = calcFreqLogLike1_MajorMinor(myReads, major, minor, maf)

like_H1 = freqsMLE[1]

lrtSnp = snpPval_MajorMinor(myReads, freqsMLE[1], major, minor)

```

- H0: -16.66699514233869
- H1: -8.899773700369137
- LRT: 15.534442883939104
- f_MLE: 0.29257488372333595






