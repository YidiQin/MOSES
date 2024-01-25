##### A function to predict gene exp for test data by model training (use meth and snp data as input) #####
library(glmnet)
prediction2 <- function(seq.num, k, inter.gene.list, promoter.range, enhancer.range, train.meth.pos.range, train.snps.pos.range, sub.train.exp, train.meth, test.meth, train.snp, test.snp, lambda.rule, n, output.file.path, chr, core.num) {
  sub.gene.list <- inter.gene.list[seq.num] # select a subset of gene
  pred.result <- mclapply(1:length(sub.gene.list),function(curi) {
    promoter.range1 <- promoter.range[promoter.range$gene.name == sub.gene.list[curi]] # find promoter loc of the tested gene
    candi.enhancer.range1<- GRanges(seqnames=seqnames(promoter.range1), IRanges(start=start(ranges(promoter.range1)) -enhancer.range,end=end(ranges(promoter.range1)) + enhancer.range),gene.name=promoter.range1$gene.name) # define enhancer region of the tested gene
    snp.id <- unique(queryHits(findOverlaps(train.snps.pos.range, candi.enhancer.range1))) # get loc of the snps overlapped with the enhancer region
    probe.id <- unique(queryHits(findOverlaps(train.meth.pos.range, candi.enhancer.range1))) # get loc of the methylation probes overlapped with the enhancer region
    np = length(c(snp.id, probe.id))

    if(length(c(snp.id, probe.id))!=0 ) {
      y=t(sub.train.exp[rownames(sub.train.exp) == sub.gene.list[curi],,drop=FALSE])
      x_snp=t(train.snp[snp.id,,drop=FALSE])
      x_meth=t(train.meth[probe.id,,drop=FALSE])
      testx_snp = t(test.snp[colnames(x_snp),])
      testx_meth = t(test.meth[colnames(x_meth),])
      if(length(snp.id)!=0 & length(probe.id)!=0) {
        x=cbind(x_snp,x_meth)
        testx=cbind(testx_snp,testx_meth)
      } else if(length(snp.id)!=0){
        x=x_snp
        testx=testx_snp
      } else {
        x=x_meth
        testx=testx_meth
      }

      set.seed(1234)

      x=apply(x,2,as.double)
      ny <- ncol(y)
      if(ny ==1) {
        y <- as.numeric(y)
        cv.fit <- cv.glmnet(x, y, keep=TRUE, alpha=0.5)
        pred <- predict(cv.fit, newx=testx,s=lambda.rule, type="response")
      } else {
        scvpearson <- rep(0, ny)
        cv.fit <-  vector("list", ny)
        for(i in 1:ny) {
          y.m <- as.numeric(y[,i])
          cv.fit[[i]] <- cv.glmnet(x, y.m, keep=TRUE, alpha=0.5)
          id<-which(cv.fit[[i]]$lambda == cv.fit[[i]]$lambda.min)

          if(var(cv.fit[[i]]$fit.preval[,id])!=0) {
            scvpearson[i] <- cor(cv.fit[[i]]$fit.preval[,id],y.m, method = "pearson")
          } else {
            scvpearson[i] <- 0
          }
        }
        ind <- which.max(abs(scvpearson))
        pred <- predict(cv.fit[[ind]], newx=testx,s=lambda.rule, type="response")
      }
    } else {
      pred <- rep(0, n)
    }
    list(pred=pred)
  },mc.cores=core.num)
  save(pred.result, file=paste(output.file.path,"chr",chr,".",k,"th.running.Rdata",sep=""))
}
