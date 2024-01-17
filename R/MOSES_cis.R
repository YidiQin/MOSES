#' Train gene expression prediction model and conduct TWAS using customized training data
#'
#' @param example An indicator to input example datasets. Default is FALSE.
#' @param train.meth.file A file containing methylation matrix as training data.
#' @param train.exp.file A file containing gene expression matrix as training data.
#' @param test.meth.file A file containing methylation matrix as test data.
#' @param TWAS An indicator to perform TWAS analysis. Default is TRUE.
#' @param pheno.file A file containing phenotype information.
#' @param pheno The name of outcome of interest in TWAS.
#' @param confounder Names of variables to be adjusted in TWAS.
#' @param output.file.path The path/directory to store output files.
#' @param core.num The number of cores to be used. Default is 1.
#' @return Two output files in the folder specified by output.file.path: 1) prediction.Rdata is a matrix containing gene expression value predicted using train.meth.rda and train.exp.rda 2) TWAS.result.txt is a matrix containing TWAS results.
#'
#' @examples
#' MethylTWAS(example = F, train.meth.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/train.meth.rda",test.meth.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/test.meth.rda",train.exp.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/train.exp.rda",pheno.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/pheno.rda",output.file.path = "/ix/ksoyeon/YQ/results/test/",TWAS = T,phenotype = "cc_new",confounder = "gender, age")
#' MethylTWAS(example = T, output.file.path = "/ix/ksoyeon/YQ/results/test/",TWAS = T,phenotype = "cc_new",confounder = "gender, age")

MOSES_cis <- function(example = TRUE, train.meth.file, train.exp.file, test.meth.file, TWAS = TRUE, pheno.file, phenotype, confounder, output.file.path, core.num = 1) {
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
  library(GenomicRanges)
  train.meth.pos.range <- MatchPos(train.meth)
  message("Matching probes ...")

  ##### find which methylation probes in testing data are also in annotation set as well as training data #####
  int_probe <- intersect(rownames(test.meth),rownames(train.meth))
  test.meth <- test.meth[rownames(test.meth) %in% int_probe,]
  train.meth <- train.meth[rownames(train.meth) %in% int_probe,]
  train.meth.pos.range <- train.meth.pos.range[train.meth.pos.range$name %in% rownames(train.meth),]

  ##### load promoter info #####
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
  enhancer.range=1e6
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
               n, output.file.path, core.num)
    load(paste(output.file.path,".",k,"th.running.Rdata",sep=""))
    sub.exp <- sapply(pred.result, function(x) x[[1]])
    colnames(sub.exp) <- inter.gene.list[seq.num]
    pred.gene.exp[match(colnames(sub.exp), rownames(pred.gene.exp)),] <- t(sub.exp)
    k <- k+1
    seq.num <- which(rowSums(pred.gene.exp) == 0)
  }
  pred.gene.exp <- pred.gene.exp
  message("Saving predicted gene expression ...")
  save(list=c('pred.gene.exp'), file=paste0(output.file.path,"prediction.Rdata"))

  ###### TWAS #####
  if(TWAS == FALSE){
    message("Skip TWAS ...")
    stop_quietly()
  }
  if(example == TRUE){
    data(pheno)
  }
  else{
    temp4 <- load(file = pheno.file)
    pheno <- get(temp4)
    rm(temp4)
  }
  message("Running TWAS ...")
  library(limma)
  confounder.var<- paste(unlist(strsplit(confounder, split = ",")),collapse="+")
  design <- model.matrix(as.formula(paste0("~0+",as.character(phenotype),"+",confounder.var)), data=pheno)
  fit <- lmFit(pred.gene.exp, design)
  a <- paste0(as.character(phenotype),"TRUE")
  b <- paste0(as.character(phenotype),"FALSE")
  contrast <- paste0("CasevsControl=",a,"-",b)
  cmd <- paste("cont.matrix <- makeContrasts(", contrast, ", levels = design)", sep ='')
  eval(parse(text = cmd))
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  imputed.TWAS <- topTable(fit2, adjust="BH",number = Inf)
  message("Saving TWAS results ...")
  write.table(imputed.TWAS, paste0(output.file.path,"TWAS.result.txt"),quote=F,sep="\t",col.names = TRUE, row.names = TRUE)
}
