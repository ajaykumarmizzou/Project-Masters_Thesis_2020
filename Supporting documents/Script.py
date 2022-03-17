---
title: "Coursera_bioconductorforGDS"
author: "AJAY KUMAR"
date: "7/5/2020"
output: pdf_document
---

```{r}
print("hello")
```
IRANGES-BASIC USAGE
```{r}
library(IRanges)
ir1 = IRanges(start= c(1,3,5),end = c(2,4,7))
ir2 = IRanges(start = c(1,5,10),width = 3)
ir1
ir2
start(ir1) #Start column of ir1
end(ir2) #end column of ir2
width(ir2)=5 #set width of the ir2
ir2
dim(ir2) #It gives NULL, it means it has no dimension
length(ir2) #It gives no. of rows..
c(ir1,ir2) #IRanges data structure can be concatenated
#Plotting of IRange intervals
plotRanges = function(x,xlim=x,main = deparse(substitute(x)),col="black", sep = 0.5)
{height =1
if(is(xlim, "Ranges"))
  xlim = c(min(start(xlim)),max(end(xlim)))
  bins = disjointBins(IRanges(start(x),end(x)+1))
  plot.new()
  plot.window(xlim,c(0,max(bins)*(height+sep)))
  ybottom = bins*(sep+height)-height
  rect(start(x)-0.5,ybottom, end(x)+0.5,ybottom+height,col=col)
  title(main)
  axis(1)
}
par(mfrow=c(2,1))
ir = IRanges(c(1,3,7,9),end = c(4,4,8,10))
plotRanges(ir)
plotRanges(reduce(ir)) #Intersection of intervals
plotRanges((disjoin(ir))) #Undo intersection operation
#Union/Intersection/Resize of IRanges-intervals
ir1
ir2
union(ir1,ir2)
intersect(ir1,ir2)
resize(c(ir1,ir2),width = 1)
plotRanges(union(ir1,ir2))
#findOverlap - Real power of IRanges #memory efficient and faster
ir1
ir2
ov = findOverlaps(ir1,ir2)
ov
queryHits(ov)
unique(queryHits(ov))
subjectHits(ov)
```
Genomic Ranges - GRanges (Extension for IRanges)
```{r}
library(GenomicRanges)
gr = GRanges(seqnames = c("chr1","chr2","chr3"),strand = c("+","-","+"),ranges = IRanges(start = c(1,3,5),width = 3))
gr
flank(gr,5) #It is relevant to the direction of transcription
promoters(gr)
seqinfo(gr) #To find the sequence information
seqlengths(gr) =c("chr1"=10,"chr2"=10,"chr3"=10) #Enter/modify length of sequence
seqinfo(gr)
seqlevels(gr)
gaps(gr) # Gives all the range of chromosome that is either covered by a range or g range
sort(gr) # To sort sequence as per chromosome no.'s
genome(gr) = "h19" #assigned genome to the chromosome
seqinfo(gr)
```
Basic GenomicRanges - Basic GRanges Usage
```{r}
library(GenomicRanges)
ir = IRanges(start = 1:3,width = 2)
ir
df = DataFrame(ir = ir, score = rnorm(3))
df
df$ir #Return IRange data structure
df
gr = GRanges(seqnames = c("chr1","chr2","chr3"),strand = c("+","-","+"),ranges = IRanges(start = c(1,3,5),width = 3))
gr
values(gr) = DataFrame(score = rnorm(3))
gr
#To access a specific column
gr$score
mcols(gr)
values(gr)
gr2= GRanges(seqnames = c("chr1","chr2","chr1"),strand = c("*"),ranges = IRanges(start = c(1,3,5),width = 3))
gr2
plotRanges(gr)
findOverlaps(gr,gr2,ignore.strand=TRUE)
df = data.frame(chr="chr1",start=1:3,end=4:6, score = rnorm(3)) #We can also make dataframe then could convert into GRange
df
makeGRangesFromDataFrame(df,keep.extra.columns = TRUE) #as simple as its name
```

Genomic Ranges -seqinfo
```{r}
gr = GRanges(seqnames = c("chr1","chr2"),ranges = IRanges(start = 1:2,end = 4:5))
gr
#mapping the seqnames 1, chr1, I into chr1
newstyle = mapSeqlevels(seqlevels(gr),"NCBI")
newstyle
gr = renameSeqlevels(gr,newstyle)
gr
```

AnnotationHub

```{r}
library(AnnotationHub)
ah = AnnotationHub()
ah
ah[1] #Know about 1st species info.
unique(ah$dataprovider) 
length(unique(ah$species))
#Select subset of the data
ah = subset(ah,species=="Homo sapiens") #Finding specific species info.
ah
query(ah,c("Homo sapiens","H3K4me3","Gm12878")) #finding info. about certain characteristics
display(ah) #see data with GUI
# fetch the data as per cell type, chromosome, species et cetra.
```

Usecase: AnnotationHub and GRanges
This is chipseq-experiment
```{r}
library(AnnotationHub)
ahub = AnnotationHub()
ahub = subset(ahub,species == "Homo sapiens")
ahub = query(ahub,c("H3K4me3","Gm12878"))
ahub
gr1 = ahub[[2]] #we are choosing second element
gr2 = ahub[[4]]
gr2
summary(width(gr1)) #how small and big a sequence is!
summary(width(gr2)) #how small and big a sequence is!
table(width(gr1))
table(width(gr2))
peaks=gr2
#Retrieving Reference genome
ahub2 = AnnotationHub()
ahub2 = query(ahub2,c("RefSeq"))
ahub2
ahub2$genome
genes=ahub2[[1]]
head(genes)
table(table(genes$name))
prom = promoters(genes)
prom
table(width(prom))
#Now we have promoters and significant seqs. Now we have to ask do we have any promoter overlap? How often they overlap?
ov=findOverlaps(prom,peaks)
length(unique(queryHits(ov)))
length(unique(subjectHits(ov)))
length((subsetByOverlaps(peaks,prom,ignore.strand=TRUE)/length(peaks)))
length((subsetByOverlaps(prom,peaks,ignore.strand=TRUE)/length(prom)))
```

Quiz:1
```{r}
library(AnnotationHub)
ah1 = AnnotationHub()
ah1=ah1[["AH5086"]]
gr=keepStandardChromosomes(ah1,pruning.mode = "tidy")
gr
```



NOTE:
dropseqlevels() and keepSeqlevels()

Biostrings - 
This is a package that contains functionality for representing and
manipulating biological strings and biodata. 

```{r}
library(Biostrings)
dna1 = DNAString("ACGTagctatcgatcgatcgtacgG") #DNA string container
dna1
dnaset1 = DNAStringSet(c("ATCGCG","AGATGAC")) #constructing DNA seq set
dnaset1
IUPAC_CODE_MAP
dna1[2:3] #subsetting  of dna-seq
dnaset1[1] #omits first element of dnaset1
dnaset1[[1]] #omit first element of first element of the dnaset1
names(dnaset1) = paste("seq",1:2) #naming to dnaset elements dna-seqs
dnaset1
sort(dnaset1) #sort the string elements
rev(dna1) #reverse the order of the elements
reverse(reverseComplement(dna1)) #biological 3 prime to 5 prime string corresponding string for dna1
translate(dna1) #converted into amino acid chain
alphabetFrequency(dna1) #exact no. of times a term occurs
letterFrequency(dna1,letters = "GC") #finding how many times a letter is coming in the sequence
dinucleotideFrequency(dna1) #same as name suggests
consensusMatrix(dnaset1) #how many times a letter is coming at a specific position
```


BSgenome - package dealing with representing full genome in Bioconductor

```{r}
library(BSgenome)  #download all genome from bioconductor not from R
available.genomes() #how many genomes are available in Bioconductor
library("BSgenome.Scerevisiae.UCSC.sacCer1")
Scerevisiae #Yeast genome
seqlengths(Scerevisiae) #get the sequence lenght
Scerevisiae$chr1
#Now you have a sequence and can perform Biostrings operations
letterFrequency(Scerevisiae$chr1,letters = "GC") # it just finds the gc no.
#if you want to find the %of GC in this the use,
letterFrequency(Scerevisiae$chr1,letters = "GC",as.prob = TRUE)
```

#BSapply-package  apply a certain function on whole list elements
```{r}
param = new("BSParams",X=Scerevisiae,FUN = letterFrequency)
bsapply(param,"GC",as.prob= TRUE) #you can add even more parameters here to perform certain operations. #apply a certain function on whole genome
gc_per_chromo=unlist(bsapply(param,"GC"))  #unlist the ouput
(sum(gc_per_chromo)/sum(seqlengths(Scerevisiae)) )*100 #finding %gc overall
```

#Biostring Matching
##We can use it in,
1. Macthing a string to a string
2. Matching a set of string to a string
3. Matching a string to a set of string
4. Macthing a set of string to a set of string

```{r}
library("BSgenome.Scerevisiae.UCSC.sacCer1")
dnaseq = DNAString("ACGTACGT")
dnaseq
matchPattern(dnaseq,Scerevisiae$chr1) #one sequence against a sequence
countPattern(dnaseq,Scerevisiae$chr1) #counts the no. of times a sequnces comes
vmatchPattern(dnaseq,Scerevisiae) #checks a sequecne against a set of sequence
matchPWM() #It also, matches, precision weight matrix; It allows to find the transcription binding sites on a genome
PairwiseAlignments() #It allows to map , millions of reads against a short sequence such as a gene. This is underutilized tool.
trimLRPatterns() #trimming specific patterns; it has very rich use of facility
```

#BSgenome - Views object - It basically, sees the snapshot of the interval sequences and stores in Iranges intervals. Thus we do not need to store sequences; just to store intervals.

```{r}
library("BSgenome.Scerevisiae.UCSC.sacCer1")
dnaseq = DNAString("ACGTACGT")
dnaseq
vi = matchPattern(dnaseq,Scerevisiae$chr1)
vi
ranges(vi)
Scerevisiae$chr1[57933:57940]
shift(vi,10) #it shifts the view/snapshot in sequence
gr = vmatchPattern(dnaseq,Scerevisiae)
vi2 = Views(Scerevisiae,gr)
vi2
ranges(vi2)
```

#GenomicRanges - RLE (Run Lenght Encoding)
##Run length encoding is a way of representing very long vectors. where some elements are the same.
```{r}
```

