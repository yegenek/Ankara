
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

```
myReads = Reads("AAATT" , "55555");
calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 4)
```

```
10^((33 - Int64('+')  )/10)
myReads = Reads("AAAT" , "555+");
calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 4)
```

```
myReads = Reads("AAATAAC" , "555-43(");
calcGenoLogLike2_MajorMinor(myReads, mySite, 1, 4)
```












