##### A function to predict gene exp for test data using coefficients #####
<<<<<<< HEAD
library(glmnet)
prediction_BySummary <- function(all.test.cvfit, seq.num, k, inter.gene.list,
                       test.meth, lambda.rule, n, output.file.path, core.num) {
  sub.gene.list <- inter.gene.list[seq.num]
  pred.result <- mclapply(1:length(sub.gene.list),function(curi) {
    print(paste(curi,"th gene",sep=""))
=======

library(glmnet)
prediction_BySummary <- function(all.test.cvfit, seq.num, k, inter.gene.list,
                       #promoter.range, enhancer.range, train.meth.pos.range, train.exp, train.meth,
                       test.meth, lambda.rule, n, output.file.path, core.num) {
  sub.gene.list <- inter.gene.list[seq.num] # select a subset of gene
  pred.result <- mclapply(1:length(sub.gene.list),function(curi) {

    print(paste(curi,"th running",sep=""))
    #promoter.range1 <- promoter.range[promoter.range$gene.name == sub.gene.list[curi]] # find promoter loc of the tested gene
    #candi.enhancer.range1<- GRanges(seqnames=seqnames(promoter.range1), IRanges(start=start(ranges(promoter.range1)) -enhancer.range,end=end(ranges(promoter.range1)) + enhancer.range),gene.name=promoter.range1$gene.name) # define enhancer region of the tested gene
    #probe.id <- unique(queryHits(findOverlaps(train.meth.pos.range, candi.enhancer.range1))) # get loc of the methylation probes overlapped with the enhancer region
    #probe.name.unique <- unique(c(promoter.probe.name, gene.body.probe.name))
    #if(length(probe.id)!=0 ) {
      #y=t(train.exp[rownames(train.exp) ==  sub.gene.list[curi],,drop=FALSE])
      #x=t(train.meth[probe.id,])
      #testx = t(test.meth[colnames(x),])
>>>>>>> 8544c7ca73e05fd4a361dc7d988f5b704acab973
      probe.id <- intersect(rownames(test.meth), names(all.test.cvfit[[curi]]))
      if(length(probe.id) == 0) {
        pred <- NA
      }
      else{
        testx = t(test.meth[probe.id,])
<<<<<<< HEAD
=======
        #testx = t(test.meth[rownames(all.test.cvfit[[curi]]$beta),])
>>>>>>> 8544c7ca73e05fd4a361dc7d988f5b704acab973
        set.seed(1234)
        testx <- cbind(1,testx)
        beta <- c(all.test.cvfit[[curi]]["intercept"], all.test.cvfit[[curi]][probe.id])
        pred  <- testx %*% beta
<<<<<<< HEAD
      }
=======
        #pred <- predict(all.test.cvfit[[curi]], newx=testx, s="lambda.min", type="response")
      }

      #x=apply(x,2,as.double)
      #ny <- ncol(y)
      #ntesty <- ncol(testy)
      #if(ny ==1) {
        #y <- as.numeric(y)
        #cv.fit <- cv.glmnet(x, y, keep=TRUE, alpha=0.5)
        #pred <- predict(cv.fit, newx=testx,s=lambda.rule, type="response")
      #} else {
        #scvpearson <- rep(0, ny)
        #cv.fit <-  vector("list", ny)
        #for(i in 1:ny) {
          #y.m <- as.numeric(y[,i])
          #cv.fit[[i]] <- cv.glmnet(x, y.m, keep=TRUE, alpha=0.5)
          #id<-which(cv.fit[[i]]$lambda == cv.fit[[i]]$lambda.min)

          #if(var(cv.fit[[i]]$fit.preval[,id])!=0) {
            #scvpearson[i] <- cor(cv.fit[[i]]$fit.preval[,id],y.m, method = "pearson")
          #} else {
            #scvpearson[i] <- 0
          #}
        #}
        #ind <- which.max(abs(scvpearson))
        #pred <- predict(cv.fit[[ind]], newx=testx,s=lambda.rule, type="response")
      #}
    #} else {
      #pred <- rep(0, n)
    #}
>>>>>>> 8544c7ca73e05fd4a361dc7d988f5b704acab973
    list(pred=pred)
  },mc.cores=core.num)
  save(pred.result, file=paste(output.file.path,".",k,"th.running.Rdata",sep=""))
}
