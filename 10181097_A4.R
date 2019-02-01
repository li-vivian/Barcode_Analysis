library(dplyr)
library(sangerseqR)
library(Biostrings)

detach("package:seqinr", unload = TRUE )


# 1. Extract and call the sequence from the .ab1 file provided 

rbcL<-read.abif("./Data/1Ipos_F_P1815443_064.ab1")

## extract
rbcLseq<-sangerseq(rbcL) 

## call
SeqX<-makeBaseCalls(rbcLseq) 


# 2. Use regular expressions to slice, paste and extract the 'Primary sequnece' only. 
primaryseq<-SeqX@primarySeq
paste(primaryseq)


# 3. Convert the sequence to a FASTA format, in which the first line starts with > and contains the file name and the second line contains the primary sequence

## first make primaryseq a data frame
primaryseq<-data.frame(primaryseq)

library(seqinr)
write.fasta(primaryseq, names = paste(" 1Ipos_F_P1815443_064"), file.out = "1Ipos_F_P1815443_064.fasta", open = "w", nbchar = 60)

## need to detach seqinr in order to use sangerseqR functions
detach("package:seqinr", unload = TRUE )


# 5. Save the BarcodePlateStats.csv to your Data folder and use it to exclude from your loop any sequences that do not pass the quality check 
BPS<-read.csv("./Data/BarcodePlateStats.csv")
true<-vector(mode = "character", length = sum(BPS$Ok == "TRUE"))

## filter through files that do not pass the quality check
q = 1

for (i in 1:length(BPS$Ok)){
  if (BPS$Ok[i] == "TRUE"){
    true[q]<- as.character(BPS$Chromatogram[i])
    q = q + 1
  }
}


# 4. Loop through every file in the Data folder, do steps 1-3 and save the output as a vector of strings in the FASTA format.

## empty vector to store good sequences and file names
goodfiles<-vector(mode = "character", length=length(true))
filename<-vector(mode = "character", length=length(true))

l = 1

for (i in true){
  ITS<-read.abif(paste0("./Data/", i))
  ITSseq<-sangerseq(ITS)
  SeqX<-makeBaseCalls(ITSseq)
  
  primseq<-primarySeq(SeqX, string = TRUE)
  
  filename[l]<-paste(i)
  
  goodfiles[l]<-paste(primseq)
  
  l = l + 1
}


## save output in FASTA format, goodfiles and filename

library(seqinr)

for (f in 1:length(filename)) {
  if (f == 1){
    write.fasta(goodfiles[f], file.out = "DNAOutput.fasta", as.string = TRUE, 
                names = paste(filename[f], "; length =", nchar(goodfiles[f]), sep = ""), open = "w")
  }
  else if (f > 1){
    write.fasta(goodfiles[f], file.out = "DNAOutput.fasta", as.string=TRUE,
                names = paste(filename[f], "; length =", nchar(goodfiles[f]), sep = ""), open = "a")
  }
}

