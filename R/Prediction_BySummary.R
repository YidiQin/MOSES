##### A function to predict gene exp for test data using coefficients #####
library(glmnet)
prediction_BySummary <- function(all.test.cvfit, seq.num, k, inter.gene.list,
                       test.meth, lambda.rule, n, output.file.path, core.num) {
  sub.gene.list <- inter.gene.list[seq.num]
  pred.result <- mclapply(1:length(sub.gene.list),function(curi) {
    print(paste(curi,"th gene",sep=""))
      probe.id <- intersect(rownames(test.meth), names(all.test.cvfit[[curi]]))
      if(length(probe.id) == 0) {
        pred <- NA
      }
      else{
        testx = t(test.meth[probe.id,])
        set.seed(1234)
        testx <- cbind(1,testx)
        beta <- c(all.test.cvfit[[curi]]["intercept"], all.test.cvfit[[curi]][probe.id])
        pred  <- testx %*% beta
      }
    list(pred=pred)
  },mc.cores=core.num)
  save(pred.result, file=paste(output.file.path,".",k,"th.running.Rdata",sep=""))
}
