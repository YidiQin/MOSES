#' Train gene expression prediction model and conduct TWAS using estimated coefficients
#'
#' @param example An indicator to input example datasets. Default is FALSE.
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
#' MethylTWAS_BySummary(example = F, test.meth.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/test.meth.rda",pheno.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/pheno.rda",output.file.path = "/ix/ksoyeon/YQ/results/test/",TWAS = T,phenotype = "cc_new",confounder = "gender, age")
#' MethylTWAS_BySummary(example = T, output.file.path = "/ix/ksoyeon/YQ/results/test/",TWAS = T,phenotype = "cc_new",confounder = "gender, age")

MethylTWAS_BySummary <- function(example = TRUE, test.meth.file, TWAS = TRUE, pheno.file, phenotype, confounder, output.file.path, core.num = 1) {
  message("Importing data ...")
  if(example == TRUE){
    data(test.meth)
  }
  else{
    temp <- load(file=test.meth.file)
    test.meth <- get(temp)
    rm(temp)
  }

  ##### load coefficient info #####
  message("Loading coeffient ...")
  githubURL <- "https://github.com/YidiQin/MethylTWAS/blob/master/data/test.10M.EVA_PR.cvfit.Rdata?raw=true"
  download.file(githubURL, destfile= "./test.10M.EVA_PR.cvfit.Rdata", mode = "wb")
  load("./test.10M.EVA_PR.cvfit.Rdata")

  ##### select genes with promoter info and in test data #####
  data("valid_genes")
  all.test.cvfit <- all.test.cvfit[names(all.test.cvfit) %in% valid_genes$x]
  if(example == TRUE){
    all.test.cvfit <- all.test.cvfit[1:10]
  }
  inter.gene.list <- names(all.test.cvfit)

  ##### prediction parameters #####
  curi<-1
  lambda.rule <- "lambda.min"
  seq.num <- 1:length(inter.gene.list)
  n <- ncol(test.meth)

  ###### prediction #####
  library(parallel)
  message("Predicting expression ...")
  k <- 1
  print(paste(k,"th.running",sep=""))
  prediction_BySummary(all.test.cvfit, seq.num, k, inter.gene.list,
                       test.meth, lambda.rule, n, output.file.path, core.num)
  load(paste(output.file.path,".",k,"th.running.Rdata",sep=""))
  sub.exp <- sapply(pred.result, function(x) x[[1]])
  names(sub.exp) <- inter.gene.list[seq.num]
  var <- sapply(sub.exp, function(x) var(x))
  bad.pre <- which(var == 0 | is.na(var))
  sub.exp <- sub.exp[-bad.pre]
  pred.gene.exp <- as.data.frame(t(matrix(unlist(sub.exp), nrow=length(unlist(sub.exp[1])))))
  rownames(pred.gene.exp) <- names(sub.exp)
  colnames(pred.gene.exp) <- colnames(test.meth)
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
  message("Saving TWAS results...")
  write.table(imputed.TWAS, paste0(output.file.path,"TWAS.result.txt"),quote=F,sep="\t",col.names = TRUE, row.names = TRUE)
}
