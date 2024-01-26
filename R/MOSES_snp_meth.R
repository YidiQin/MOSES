#' Train gene expression prediction model and conduct TWAS using customized training data
#'
#' @param example An indicator to input example datasets. Default is FALSE.
#' @param train.meth.file A file containing methylation matrix as training data.
#' @param train.exp.file A file containing gene expression matrix as training data.
#' @param test.meth.file A file containing methylation matrix as test data.
#' @param genotype.file.list.file A txt file containing three columns named as "chr", "train.snp.file", and "test.snp.file". The "chr" column contains numeric chromosome number (e.g. 1, 2, ..., 22); the "train.snp.file" contains the directory of the training set genotype file from the corresponding chr; the "train.snp.file" contains the directory of the test set genotype file from the corresponding chr
#' @param TWAS An indicator to perform TWAS analysis. Default is TRUE.
#' @param pheno.file A file containing phenotype information.
#' @param pheno The name of outcome of interest in TWAS.
#' @param confounder Names of variables to be adjusted in TWAS.
#' @param output.file.path The path/directory to store output files.
#' @param core.num The number of cores to be used. Default is 1.
#' @return Two output files in the folder specified by output.file.path: 1) prediction.Rdata is a matrix containing gene expression value predicted using train.meth.rda and train.exp.rda 2) TWAS.result.txt is a matrix containing TWAS results.
#'
#' @examples
#' MOSES_snp_meth(example = F, train.meth.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/train.meth.rda",test.meth.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/test.meth.rda",train.exp.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/train.exp.rda",genotype.file.list="/ix/ksoyeon/YQ/code/MethylTWAS/data/genotype.file.list.txt",pheno.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/pheno.rda",output.file.path = "/ix/ksoyeon/YQ/results/test/",TWAS = T,phenotype = "cc_new",confounder = "gender, age")
#' MOSES_snp_meth(example = T, output.file.path = "/ix/ksoyeon/YQ/results/test/",TWAS = T,phenotype = "cc_new",confounder = "gender, age")

MOSES_snp_meth <- function(example = TRUE, train.meth.file, train.exp.file, test.meth.file, genotype.file.list.file, TWAS = TRUE, pheno.file, phenotype, confounder, output.file.path, core.num = 1) {
  message("Importing data ...")
  if(example == TRUE){
    data(train.meth.1)
    data(train.meth.2)
    train.meth <- rbind(train.meth.1, train.meth.2)
    data(train.exp)
    data(test.meth)
    data(genotype.file.list)
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
    genotype.file.list <- read.table(genotype.file.list.file)
  }

  train.exp <- GenePos(train.exp)$gene.exp
  train.gene.pose <- GenePos(train.exp)$gene.pose

  library(GenomicRanges)
  train.meth.pos.range <- MatchPos(train.meth)
  message("Matching probes ...")

  ##### find which methylation probes in testing data are also in annotation set as well as training data #####
  int_probe <- intersect(rownames(test.meth),rownames(train.meth))
  test.meth <- test.meth[rownames(test.meth) %in% int_probe,]
  train.meth <- train.meth[rownames(train.meth) %in% int_probe,]
  train.meth.pos.range <- train.meth.pos.range[train.meth.pos.range$name %in% rownames(train.meth),]

  ##### get genotype file information ####
  chr.list <- as.numeric(genotype.file.list$chr)

  ###### prediction by chr #####
  library(parallel)
  library(glmnet)
  library(data.table)
  message("Predicting expression ...")
  for(chr in chr.list) {

    train.id <- fread(genotype.file.list$train.snp.file[genotype.file.list$chr == chr],na.strings=c("",NA,"NULL"))[,1]
    train.SNPs <- fread(genotype.file.list$train.snp.file[genotype.file.list$chr == chr],na.strings=c("",NA,"NULL"))[,-1]
    test.id <- fread(genotype.file.list$test.snp.file[genotype.file.list$chr == chr],na.strings=c("",NA,"NULL"))[,1]
    test.SNPs <- fread(genotype.file.list$test.snp.file[genotype.file.list$chr == chr],na.strings=c("",NA,"NULL"))[,-1]

    # train snp set
    train.snp <- transpose(train.SNPs)
    train.snp <- as.matrix(train.snp)
    colnames(train.snp) <- t(train.id)
    rownames(train.snp) <- colnames(train.SNPs)
    train.sub1 <- sub("\\(?[0-9]+\\:", "", rownames(train.snp))
    train.pos <- sub("_.","", train.sub1)
    train.snps.pos <- data.frame(snpid=rownames(train.snp), chr=chr, pos=train.pos)
    train.snps.pos$pos <- as.numeric(as.character(train.snps.pos$pos))
    library(GenomicRanges)
    train.snps.pos.range <- GRanges(seqnames = paste("chr",train.snps.pos$chr,sep=""),
                                    ranges = IRanges(train.snps.pos$pos, width = 1, names =rownames(train.snp)),
                                    name = rownames(train.snps.pos))

    # test snp set
    test.snp <- transpose(test.SNPs)
    test.snp <- as.matrix(test.snp)
    colnames(test.snp) <- t(test.id)
    rownames(test.snp) <- colnames(test.SNPs)
    test.sub1 <- sub("\\(?[0-9]+\\:", "", rownames(test.snp))
    test.pos <- sub("_.","", test.sub1)
    test.snps.pos <- data.frame(snpid=rownames(test.snp), chr=chr, pos=test.pos)
    test.snps.pos$pos <- as.numeric(as.character(test.snps.pos$pos))
    test.snps.pos.range <- GRanges(seqnames = paste("chr",test.snps.pos$chr,sep=""),
                                    ranges = IRanges(test.snps.pos$pos, width = 1, names =rownames(test.snp)),
                                    name = rownames(test.snps.pos))

    # train and test snp sets with overlapped snps
    int_snp <- intersect(rownames(test.snp),rownames(train.snp))
    test.snp <- test.snp[rownames(test.snp) %in% int_snp,]
    train.snp <- train.snp[rownames(train.snp) %in% int_snp,]
    train.snps.pos.range <- train.snps.pos.range[train.snps.pos.range$name %in% rownames(train.snp),]

    # find genes on specific chr
    sub.train.exp <- train.exp[train.gene.pose$chr == paste0("chr",chr),]

    # load promoter info #####
    data(promoter)
    promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)

    # select genes with promoter info and in training data #####
    inter.gene.list <-promoter.range$gene.name[promoter.range$gene.name %in% rownames(sub.train.exp)]
    name<- rownames(sub.train.exp)
    dup.name<- name[duplicated(name)]

    # prediction parameters
    curi<-1
    nfolds<-10
    lambda.rule <- "lambda.min"
    enhancer.range=1e6
    seq.num <- 1:length(inter.gene.list)
    n <- ncol(test.meth)

    # prediction
    pred.gene.exp <- matrix(NA, ncol=n, nrow=length(inter.gene.list))
    rownames(pred.gene.exp) <- inter.gene.list
    colnames(pred.gene.exp) <- colnames(test.meth)
    k<-1
    while(length(seq.num) !=0 & k < 11) {
      print(paste(k,"th.running",sep=""))
      prediction2(seq.num, k, inter.gene.list, promoter.range, enhancer.range,
                 train.meth.pos.range, train.snps.pos.range, sub.train.exp,
                 train.meth, test.meth, train.snp, test.snp, lambda.rule,
                 n, output.file.path, chr, core.num)
      load(paste(output.file.path,"chr",chr,".",k,"th.running.Rdata",sep=""))
      sub.exp <- sapply(pred.result, function(x) x[[1]])
      colnames(sub.exp) <- inter.gene.list[seq.num]
      pred.gene.exp[match(colnames(sub.exp), rownames(pred.gene.exp)),] <- t(sub.exp)
      k <- k+1
      seq.num <- which(rowSums(pred.gene.exp) == 0)
    }
    sub.pred.gene.exp <- pred.gene.exp
    message(paste0("Saving predicted gene expression for chr",chr," ..."))
    save(list=c('sub.pred.gene.exp'), file=paste0(output.file.path,"chr",chr,".prediction.Rdata"))
  }

  ###### Merge predicted gene expression across all chr #####
  message("Merge predicted gene expression across all chr...")
  pred.gene.exp <- NULL
  for(chr in chr.list) {
    load(paste0(output.file.path,"chr",chr,".prediction.Rdata"))
    pred.gene.exp <- rbind(pred.gene.exp, sub.pred.gene.exp)
  }
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
