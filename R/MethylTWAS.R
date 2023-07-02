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
#load("/ix/ksoyeon/YQ/code/MethylTWAS/data/promoter.rda")

#library(usethis)
#use_data(train.meth.1, overwrite = TRUE)
#use_data(train.meth.2, overwrite = TRUE)
#use_data(test.meth, overwrite = TRUE)
#use_data(train.exp, overwrite = TRUE)
#use_data(promoter, overwrite = TRUE)

MethylTWAS <- function(example = FALSE, train.meth.file, train.exp.file, test.meth.file, TWAS = TRUE, pheno.file, predictor, confounder, output.file.path) {
  message("Importing data ...")
  if(example == TRUE){
    data(train.meth.1)
    data(train.meth.2)
    train.meth <- rbind(train.meth.1, train.meth.2)
    data(train.exp)
    data(test.meth)
  }
  else{
    temp1 <- load(file=train.meth.file)
    train.meth <- get(temp1)
    rm(temp1)
    temp2 <- load(file=train.exp.file)
    train.exp <- get(temp2)
    rm(temp2)
    temp3 <- load(file=test.meth.file)
    test.meth <- get(temp3)
    rm(temp3)
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
  #load("/ix/ksoyeon/YQ/code/MethylTWAS/data/promoter.rda")
  data(promoter)
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
               train.meth.pos.range, train.exp, train.meth, test.meth, lambda.rule,
               n, output.file.path)
    load(paste(output.file.path,".",k,"th.running.Rdata",sep=""))
    sub.exp <- sapply(pred.result, function(x) x[[1]])
    colnames(sub.exp) <- inter.gene.list[seq.num]
    pred.gene.exp[match(colnames(sub.exp), rownames(pred.gene.exp)),] <- t(sub.exp)
    k <- k+1
    seq.num <- which(rowSums(pred.gene.exp) == 0)
  }
  pred.gene.exp <- pred.gene.exp
  save(list=c('pred.gene.exp'), file=paste0(output.file.path,"prediction.Rdata"))
  message("Saving predicted gene expression ...")

  ###### TWAS #####
  if(TWAS == FALSE){
    stop("\r skip TWAS ...")
  }
  if(example == TRUE){
    data(pheno)
  }
  else{
    temp4 <- load(file = pheno.file)
    pheno <- get(temp4)
    rm(temp4)
    #pheno <- read.table(pheno.file,sep="\t", header=TRUE)
  }
  library(limma)
  #print('pass000')
  confounder.var<- paste(unlist(strsplit(confounder, split = ",")),collapse="+")
  #print('pass00')
  #formula <- paste0("~0+",as.character(predictor),"+",confounder.var)
  #cmd0 <- paste("design <- model.matrix(as.formula(", formula, "), data=pheno)", sep = '')
  #eval(parse(text = cmd0))
  design <- model.matrix(as.formula(paste0("~0+",as.character(predictor),"+",confounder.var)), data=pheno)
  #print('pass0')
  fit <- lmFit(pred.gene.exp, design)
  #print("pass1")
  a <- paste0(as.character(predictor),"TRUE")
  b <- paste0(as.character(predictor),"FALSE")
  #cont.matrix <- makeContrasts(paste0("CasevsControl=",a,"-",b), levels=design)
  #print("pass2")
  contrast <- paste0("CasevsControl=",a,"-",b)
  cmd <- paste("cont.matrix <- makeContrasts(", contrast, ", levels = design)", sep ='')
  eval(parse(text = cmd))
  #print("pass3")
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  imputed.TWAS <- topTable(fit2, adjust="BH",number = Inf)
  write.table(imputed.TWAS, paste0(output.file.path,"TWAS.result.txt"),quote=F,sep="\t",col.names = TRUE, row.names = TRUE)
}

#MethylTWAS(train.meth.file = "/data/Yang.meth.train.txt", test.meth.file = "/ix/ksoyeon/YQ/data/Yang.meth.test.txt", train.exp.file = "/ix/ksoyeon/YQ/data/Yang.exp.train.txt",pheno.file="/ix/ksoyeon/YQ/data/Yang.pheno.txt",output.file.path = "/ix/ksoyeon/YQ/results/test/")
