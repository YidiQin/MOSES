#load("/ix/ksoyeon/eQTMs/GSE65205.data.Rdata") # observed expression for test

#write.table(mvalue[,1:56], "/ix/ksoyeon/YQ/data/Yang.meth.train.txt",quote=F,sep="\t",col.names = TRUE, row.names = TRUE)
#write.table(mvalue[,57:69], "/ix/ksoyeon/YQ/data/Yang.meth.test.txt",quote=F,sep="\t",col.names = TRUE, row.names = TRUE)
#write.table(nonna.gene.exp[unique(row.names(nonna.gene.exp))[1:50],1:56], "/ix/ksoyeon/YQ/data/Yang.exp.train.txt",quote=F,sep="\t",col.names = TRUE, row.names = TRUE)
#write.table(pheno[57:69,],"/ix/ksoyeon/YQ/data/Yang.pheno.txt",quote=F,sep="\t",col.names = TRUE, row.names = TRUE)

#train.meth.file <- "/ix/ksoyeon/YQ/code/MethylTWAS/data/train.meth.rda"
#test.meth.file <- "/ix/ksoyeon/YQ/data/Yang.meth.test.txt"
#train.exp.file <- "/ix/ksoyeon/YQ/data/Yang.exp.train.txt"
#output.file.path <- "/ix/ksoyeon/YQ/results/test/"
#pheno.file <-"/ix/ksoyeon/YQ/data/Yang.pheno.txt"

#train.meth.1 <- mvalue[1:200000,1:56]
#train.meth.2 <- mvalue[200001:423516,1:56]
#test.meth <- mvalue[,57:69]
#train.exp <- nonna.gene.exp[unique(row.names(nonna.gene.exp))[1:50],1:56]

#library(usethis)
#use_data(train.meth.1, overwrite = TRUE)
#use_data(train.meth.2, overwrite = TRUE)
#use_data(test.meth, overwrite = TRUE)
#use_data(train.exp, overwrite = TRUE)

MethylTWAS <- function(example, train.meth.file, train.exp.file, test.meth.file, pheno.file, output.file.path) {
  message("Importing data ...")
  if(example == TRUE){
    data(train.meth.1)
    data(train.meth.2)
    train.meth <- rbind(train.meth.1, train.meth.2)
    data(train.exp)
    data(test.meth)
  }
  else{
    load(file=train.meth.file)
    load(file=train.exp.file)
    load(file=test.meth.file)
  }
  #train.meth <- read.table(train.meth.file,sep="\t", header=TRUE)
  library(GenomicRanges)
  train.meth.pos.range <- MatchPos(train.meth)
  #train.exp <- read.table(train.exp.file,sep="\t", header=TRUE)
  #test.meth <- read.table(test.meth.file,sep="\t", header=TRUE)

  message("Matching probes ...")
  ##### find which methylation probes in testing data are also in annotation set as well as training data #####
  int_probe <- intersect(rownames(test.meth),rownames(train.meth))
  test.meth <- test.meth[rownames(test.meth) %in% int_probe,]
  train.meth <- train.meth[rownames(train.meth) %in% int_probe,]
  train.meth.pos.range <- train.meth.pos.range[train.meth.pos.range$name %in% rownames(train.meth),]

  ##### load promoter info #####
  #promoter<-read.delim("/ix/ksoyeon/eQTMs/hg19_promoter.txt")
  load("/ix/ksoyeon/YQ/code/MethylTWAS/data/promoter.rda")
  promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)

  ##### select genes with promoter info and in training data #####
  inter.gene.list <-promoter.range$gene.name[promoter.range$gene.name %in% rownames(train.exp)]
  name<- rownames(train.exp)
  dup.name<- name[duplicated(name)]

  ##### prediction parameters #####
  curi<-1
  nfolds<-10
  lambda.rule <- "lambda.min"
  enhancer.range=1e7
  seq.num <- 1:length(inter.gene.list)
  n <- ncol(test.meth)

  ###### prediction #####
  library(parallel)
  library(glmnet)
  message("Predicting expression ...")
  pred.gene.exp <- matrix(NA, ncol=n, nrow=length(inter.gene.list))
  rownames(pred.gene.exp) <- inter.gene.list
  colnames(pred.gene.exp) <- colnames(test.meth)
  k<-1
  while(length(seq.num) !=0 & k < 11) {
    print(paste(k,"th.running",sep=""))
    prediction(seq.num, k, inter.gene.list, promoter.range, enhancer.range,
               train.meth.pos.range, train.exp, train.meth, test.meth, lambda.rule, n, output.file.path)
    load(paste(output.file.path,".",k,"th.running.Rdata",sep=""))
    sub.exp <- sapply(pred.result, function(x) x[[1]])
    colnames(sub.exp) <- inter.gene.list[seq.num]
    pred.gene.exp[match(colnames(sub.exp), rownames(pred.gene.exp)),] <- t(sub.exp)
    k <- k+1
    seq.num <- which(rowSums(pred.gene.exp) == 0)
  }
  pred.gene.exp <- pred.gene.exp
  save(list=c('pred.gene.exp'), file=paste0(output.file.path,"prediction.Rdata"))

  ###### TWAS #####
  if(example == TRUE){
    data(pheno)
  }
  else{
    pheno <- read.table(pheno.file,sep="\t", header=TRUE)
  }
  library(limma)
  design <- model.matrix(~0+cc_new+gender+age, data=pheno)
  fit <- lmFit(pred.gene.exp, design)
  cont.matrix <- makeContrasts(AtopicvsControl=cc_newTRUE-cc_newFALSE, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  imputed.TWAS <- topTable(fit2, adjust="BH",number = Inf)
  write.table(imputed.TWAS, paste0(output.file.path,"TWAS.result.txt"),quote=F,sep="\t",col.names = TRUE, row.names = TRUE)
}

#MethylTWAS(train.meth.file = "/data/Yang.meth.train.txt", test.meth.file = "/ix/ksoyeon/YQ/data/Yang.meth.test.txt", train.exp.file = "/ix/ksoyeon/YQ/data/Yang.exp.train.txt",pheno.file="/ix/ksoyeon/YQ/data/Yang.pheno.txt",output.file.path = "/ix/ksoyeon/YQ/results/test/")
