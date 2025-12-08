counts <- function(out, alpha, method=c("aic","bic","cit",
                                        "par.cmst.joint.aic","par.cmst.aic","non.par.cmst.aic",
                                        "par.cmst.joint.bic","par.cmst.bic","non.par.cmst.bic"))
{
  ###
  get.counts.1 <- function(M, alpha)
  {
    M1 <- sum((M[,1] <= alpha) & (M[,2] > alpha) & (M[,3] > alpha) & (M[,4] > alpha))
    M2 <- sum((M[,1] > alpha) & (M[,2] <= alpha) & (M[,3] > alpha) & (M[,4] > alpha))
    M3 <- sum((M[,1] > alpha) & (M[,2] > alpha) & (M[,3] <= alpha) & (M[,4] > alpha))
    M4 <- sum((M[,1] > alpha) & (M[,2] > alpha) & (M[,3] > alpha) & (M[,4] <= alpha))
    no.call <- nrow(M) - M1 - M2 - M3 - M4
    output <- data.frame(M1,M2,M3,M4,no.call)
  }
  ###
  get.counts.2 <- function(M, alpha)
  {
    M1 <- sum((M[,1] <= alpha) & (M[,2] > alpha))
    M2 <- sum((M[,1] > alpha) & (M[,2] <= alpha))
    M3 <- sum((M[,1] > alpha) & (M[,2] > alpha))
    no.call <- nrow(M) - M1 - M2 - M3
    output <- data.frame(M1,M2,M3,no.call)
  }
  ###
  get.counts.3 <- function(M)
  {
    MM <- t(apply(M,1,rank))
    M1 <- sum(MM[,1]==1)
    M2 <- sum(MM[,2]==1)
    M3 <- sum(MM[,3]==1)
    M4 <- sum(MM[,4]==1)
    no.call <- nrow(MM) - M1 - M2 - M3 - M4
    output <- data.frame(M1,M2,M3,M4,no.call)
  }
  ###
  if(method=="par.cmst.joint.aic")
    output <- get.counts.1(out$pval.par.cmst.joint.AIC, alpha)
  if(method=="par.cmst.aic")
    output <- get.counts.1(out$pval.par.cmst.iu.AIC, alpha)
  if(method=="non.par.cmst.aic")
    output <- get.counts.1(out$pval.non.par.cmst.iu.AIC, alpha)
  if(method=="par.cmst.joint.bic")
    output <- get.counts.1(out$pval.par.cmst.joint.BIC, alpha)
  if(method=="par.cmst.bic")
    output <- get.counts.1(out$pval.par.cmst.iu.BIC, alpha)
  if(method=="non.par.cmst.bic")
    output <- get.counts.1(out$pval.non.par.cmst.iu.BIC, alpha)
  if(method=="cit")
    output <- get.counts.2(out$pval.cit, alpha)
  if(method=="aic")
    output <- get.counts.3(out$AICs)
  if(method=="bic")
    output <- get.counts.3(out$BICs)
  ###
  output
}
