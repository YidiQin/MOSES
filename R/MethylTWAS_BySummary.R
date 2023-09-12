MethylTWAS_BySummary <- function(example, test.meth.file, TWAS = TRUE, pheno.file, phenotype, confounder, output.file.path, core.num = 1) {
  message("Importing data ...")
  if(example == TRUE){
    data(test.meth)
  }
  else{
    temp <- load(file=test.meth.file)
    test.meth <- get(temp)
    rm(temp)
  }
  #library(GenomicRanges)
  #test.meth.pos.range <- MatchPos(test.meth)

  #message("Matching probes ...")
  ##### find which methylation probes in testing data are also in annotation set as well as training data #####
  #int_probe <- intersect(rownames(test.meth),rownames(train.meth))
  #test.meth <- test.meth[rownames(test.meth) %in% int_probe,]
  #train.meth <- train.meth[rownames(train.meth) %in% int_probe,]
  #train.meth.pos.range <- train.meth.pos.range[train.meth.pos.range$name %in% rownames(train.meth),]

  ##### load promoter info #####
  #promoter<-read.delim("/ix/ksoyeon/eQTMs/hg19_promoter.txt")
  #load("/ix/ksoyeon/YQ/code/MethylTWAS/data/promoter.rda")
  #promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start=promoter$start, end=promoter$end), strand = promoter$strand, gene.name =promoter$gene.name)

  ##### load coefficient info #####
  message("Loading coeffient ...")
  githubURL <- "https://github.com/YidiQin/MethylTWAS/blob/master/data/test.10M.EVA_PR.cvfit.Rdata?raw=true"
  download.file(githubURL, destfile= "./test.10M.EVA_PR.cvfit.Rdata", mode = "wb")
  load("./test.10M.EVA_PR.cvfit.Rdata")

  ##### select genes with promoter info and in test data #####
  inter.gene.list <- names(all.test.cvfit)

  #inter.gene.list <-promoter.range$gene.name[promoter.range$gene.name %in% colnames(test.meth)]
  #name<- rownames(train.exp)
  #dup.name<- name[duplicated(name)]

  ##### prediction parameters #####
  curi<-1
  #nfolds<-10
  lambda.rule <- "lambda.min"
  #enhancer.range=1e7
  seq.num <- 1:length(inter.gene.list)
  n <- ncol(test.meth)

  ###### prediction #####
  library(parallel)
  message("Predicting expression ...")
  #pred.gene.exp <- matrix(NA, ncol=n, nrow=length(inter.gene.list))
  #rownames(pred.gene.exp) <- inter.gene.list
  #colnames(pred.gene.exp) <- colnames(test.meth)
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
    #pheno <- read.table(pheno.file,sep="\t", header=TRUE)
  }
  library(limma)
  #print('pass000')
  confounder.var<- paste(unlist(strsplit(confounder, split = ",")),collapse="+")
  #print('pass00')
  #formula <- paste0("~0+",as.character(phenotype),"+",confounder.var)
  #cmd0 <- paste("design <- model.matrix(as.formula(", formula, "), data=pheno)", sep = '')
  #eval(parse(text = cmd0))
  design <- model.matrix(as.formula(paste0("~0+",as.character(phenotype),"+",confounder.var)), data=pheno)
  #print('pass0')
  fit <- lmFit(pred.gene.exp, design)
  #print("pass1")
  a <- paste0(as.character(phenotype),"TRUE")
  b <- paste0(as.character(phenotype),"FALSE")
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
