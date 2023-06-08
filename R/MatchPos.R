MatchPos <- function(data){
  loc <- read.table("/ix/ksoyeon/geneEXPLORE/Annotation_Probe_Location.txt",header=TRUE)
  loc2 <- loc[loc$ProbeID%in%rownames(data),]
  loc2$ProbeID<-as.character(loc2$ProbeID)
  match.loc = loc2[match(rownames(data), loc2$ProbeID),c(1,4,5)]
  meth.pos <- data.frame(probeid=match.loc$ProbeID, chr=substring(match.loc$CHROMOSOME,4), pos=match.loc$POSITION)
  meth.pos.range <- GRanges(seqnames = paste("chr",meth.pos$chr,sep=""), ranges = IRanges(start = meth.pos$pos, width = 1, names =rownames(meth.pos)), name = rownames(data))
  return(meth.pos.range)
}

