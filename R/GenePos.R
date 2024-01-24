
GenePos <- function(data){
  data("TSS")
  syn.name<- synonym.search(TSS$gene,rownames(data))
  org.name<-rownames(data)
  n <- sum(!org.name %in% syn.name)
  message(paste0("A total number of ",n," genes do not have TSS information and will be removed from the analysis"))
  rownames(data) <- syn.name

  genepos <- data.frame(geneid=TSS$gene, chr=TSS$chr, left=TSS$start, right=TSS$start)
  sub.gene.pose <- genepos[as.character(genepos$geneid)  %in% rownames(data),]
  sub.gene.pose2 <- sub.gene.pose[match(rownames(data), as.character(sub.gene.pose$geneid)),]
  gene.exp <- data[!is.na(sub.gene.pose2$geneid),] # 25911 * 466
  gene.pose <- sub.gene.pose2[!is.na(sub.gene.pose2$geneid),]
  return(list(gene.exp = gene.exp, gene.pose = gene.pose))
}
