library("GenomicFeatures")
library("BSgenome.Mmusculus.UCSC.mm10")
library(Biostrings)

#load txDb of refGene.sqlite
refGene <- loadDb("refGene.sqlite")
refGene

tx_seqs1 <- extractTranscriptSeqs(Mmusculus, refGene)
cds_seqs <- extractTranscriptSeqs(Mmusculus, cdsBy(refGene, by="tx", use.names=TRUE))
protein <- translate(cds_seqs)

#write FASTA files
writeXStringSet(cds_seqs, "tmp/mm10_refGene.cds.fa")
writeXStringSet(protein, "tmp/mm10_refGene.protein.fa")

#grep all ATG
atg <- gregexpr2("ATG", cds_seqs)
names(atg) <- names(cds_seqs)

#some numbers regarding CDS missing ATG
cds_seqs[sapply(atg, function(x) { x[1] != 1 })] #none ATG starting CDS!!!!
length(grep("_.*_", names(cds_seqs[sapply(atg, function(x) { x[1] != 1 })])))
length(grep("_.*_", names(cds_seqs[sapply(atg, function(x) { x[1] != 1 })]), invert=TRUE))
#not in frame and no ATG
intersect(names(cds_seqs[(width(cds_seqs) %% 3 != 0 )]), names(cds_seqs[sapply(atg, function(x) { x[1] != 1 })]))

#extract in frame ATGs (resulting in a similar but shortend protein)
atg.frame <- lapply(atg,  function(x) {
  frame <- (x-1) %% 3
  x[frame != 0] <- 0
  return(x)  
})

#extract 2nd ATG which is in frame
atg.frame.2nd <- unlist(sapply(atg.frame, function(x) {
    #missing 1st ATG, hence look at first in list
    if (x[1] != 1)
    {
      #is it inframe else return NULL  
      if (x[1] != 0)
        {
          return(x[1])
        } else {
          return(NULL)
        }
    }
  #is 2nd ATG inframe else return NULL  
  if (length(x) > 1 && x[2] != 0)
    { 
      return(x[2])
    } else {
      return(NULL)
    }
}))

length(atg.frame.2nd)
#CDS is dividable by 3
atg.frame.2nd <- atg.frame.2nd[setdiff(names(atg.frame.2nd), names(cds_seqs[(width(cds_seqs) %% 3 != 0 )]))]

summary(atg.frame.2nd)
tail(sort(atg.frame.2nd))
write.table(cbind(refSeq=names(atg.frame.2nd), position2ndATG=atg.frame.2nd), file = "tmp/ATG2nd.txt", quote = FALSE, row.names = FALSE, sep = "\t")

