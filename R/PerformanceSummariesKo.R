PerformanceSummariesKo <- function(alpha, nms, val.targets, all.orfs, 
                                   tests, cis.index)
PerformanceSummariesKo <- function(alpha, nms, val.targets, all.orfs,
                                   tests, cis.index)
{
  rnms <- c("aic", "bic", "j.bic", "p.bic", "np.bic", "j.aic", "p.aic", "np.aic", "cit")
  nt <- length(nms)
  TP <- FP <- TN <- FN <- NC <- data.frame(matrix(0, 9, nt))
  tar <- rep(NA, nt)
  names(TP) <- names(FP) <- names(TN) <- names(FN) <- names(NC) <- nms
  row.names(TP) <- row.names(FP) <- row.names(TN) <- row.names(FN) <- row.names(NC) <- rnms
  Causal <- NotCausal <- vector(mode = "list", length = 9)
  for (k in 1 : nt) {
    aux <- which(tests[[11]][,1] == nms[k])
    aux.nms <- tests[[11]][aux, 2]
    tar[k] <- length(aux)
    for (i in 2 : 3) {
      aux.rank <- apply(tests[[i]][aux, 1:4, drop = F], 1, rank)
      aux.best <- apply(aux.rank, 2, function(x) which.min(x))
      aux.index <- as.numeric(which(aux.best == 1))
      Causal[[i - 1]] <- aux.nms[aux.index]
      NotCausal[[i - 1]] <- aux.nms[-aux.index]
    }
    for (i in 4 : 9) {
      Causal[[i - 1]] <- aux.nms[which(tests[[i]][aux, 1] <= alpha)]
      NotCausal[[i - 1]] <-
        aux.nms[c(which(tests[[i]][aux, 2] <= alpha),
                  which(tests[[i]][aux, 3] <= alpha),
                  which(tests[[i]][aux, 4] <= alpha))]
    }
    Causal[[9]] <-
      aux.nms[which(tests[[10]][aux, 1] <= alpha & tests[[10]][aux, 2] > alpha)]
    NotCausal[[9]] <-
      aux.nms[c(which(tests[[10]][aux, 1] > alpha & tests[[10]][aux, 2] <= alpha),
                which(tests[[10]][aux, 1] >= alpha & tests[[10]][aux, 2] >= alpha))]
    val <- val.targets[[match(nms[k], names(val.targets))]]
    not.val <- all.orfs[-match(unique(c(nms[k], val)), all.orfs)]
    for (i in 1 : 9) {
      TP[i, k] <- length(!is.na(intersect(Causal[[i]], val)))
      FP[i, k] <- length(!is.na(intersect(Causal[[i]], not.val)))
      TN[i, k] <- length(!is.na(intersect(NotCausal[[i]], not.val)))
      FN[i, k] <- length(!is.na(intersect(NotCausal[[i]], val)))
    }
    for (i in 4 : 9) {
      NC[i - 1, k] <- length(c(which(tests[[i]][aux, 1] > alpha),
                               which(tests[[i]][aux, 2] > alpha),
                               which(tests[[i]][aux, 3] > alpha),
                               which(tests[[i]][aux, 4] > alpha)))
    }
    NC[9, k] <- length(c(which(tests[[10]][aux, 1] < alpha),
                         which(tests[[10]][aux, 2] < alpha)))
  }
  tp <- apply(TP, 1, sum)
  fp <- apply(FP, 1, sum)
  tn <- apply(TN, 1, sum)
  fn <- apply(FN, 1, sum)
  nc <- apply(NC, 1, sum)
  prec <- tp/(tp + fp)
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  overall.1 <- data.frame(prec, tp, fp, tpr, fpr, tn, fn, nc)
  tp <- apply(TP[, cis.index], 1, sum)
  fp <- apply(FP[, cis.index], 1, sum)
  tn <- apply(TN[, cis.index], 1, sum)
  fn <- apply(FN[, cis.index], 1, sum)
  nc <- apply(NC[, cis.index], 1, sum)
  prec <- tp/(tp + fp)
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  overall.2 <- data.frame(prec, tp, fp, tpr, fpr, tn, fn, nc)
  list(overall.1, overall.2, tar)
}

##############################################################################
performance.summaries.cmst <- function(out, model, alpha=0.05, method)
{
  ntests <- nrow(out[[1]])
  if((model=="A") || (model=="B")){
    ct <- counts(out, alpha, method)
    tp <- ct[1,1]
    fp <- ct[1,2]+ct[1,3]+ct[1,4]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)
  } 
  if(model=="C"){
    ct <- counts(out, alpha, method)
    tp <- ct[1,4]
    fp <- ct[1,1]+ct[1,2]+ct[1,3]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp) 
  }
  if(model=="D"){
    ct <- counts(out, alpha, method)
    tp <- ct[1,3]
    fp <- ct[1,1]+ct[1,2]+ct[1,4]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)  
  }
  if(model=="E"){
    ct <- counts(out, alpha, method)
    tp <- ct[1,4]
    fp <- ct[1,1]+ct[1,2]+ct[1,3]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)   
  } 
  data.frame(tp, fp, power, type1.err, prec)
}
##############################################################################
performance.summaries.cit <- function(out, model, alpha=0.05)
{
  ntests <- nrow(out[[1]])
  if((model=="A") || (model=="B")){
    ct <- counts(out, alpha, method="cit")
    tp <- ct[1,1]
    fp <- ct[1,2]+ct[1,3]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)
  } 
  if((model=="C") || (model=="D") || (model=="E")){
    ct <- counts(out, alpha, method="cit")
    tp <- ct[1,3]
    fp <- ct[1,1]+ct[1,2]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)   
  }
  data.frame(tp, fp, power, type1.err, prec)
}
##############################################################################
get.power.type1.prec.matrix <- function(out, models, alpha)
{
  n <- length(alpha)
  Power <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
                                          "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                        as.character(alpha)))
  Type1 <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
                                          "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                        as.character(alpha)))
  Prec <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
                                         "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                       as.character(alpha)))
  for(k in 1:n){
    outs <- array(NA, c(9,2,5), dimnames=list(c("aic","par.joint.aic","par.aic",
                                                "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                              c("TP","FP"),models))
    for(i in 1:5){
      outs[1,1:2,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="aic")[1:2])
      outs[2,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.aic")[1:2])
      outs[3,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.aic")[1:2])
      outs[4,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.aic")[1:2])
      outs[5,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="bic")[1:2])
      outs[6,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.bic")[1:2])
      outs[7,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.bic")[1:2])
      outs[8,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.bic")[1:2])
      outs[9,,i] <- as.numeric(performance.summaries.cit(out[[i]], model=models[i], alpha[k])[1:2])
    }
    all <- matrix(NA, 9, 2, dimnames=list(c("aic","par.joint.aic","par.aic",
                                            "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                          c("TP","FP")))
    for(i in 1:9)
      all[i,] <- apply(outs[i,1:2,],1,sum)
    Prec[,k] <- all[,1]/apply(all,1,sum)
    Power[,k] <- all[,1]/5000
    Type1[,k] <- all[,2]/5000
    print(k)
  }
  list(Power=Power, Type1=Type1, Prec=Prec)
}
##############################################################################
## without model C
get.power.type1.prec.matrix.2 <- function(out, models, alpha)
{
  n <- length(alpha)
  Power <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
                                          "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                        as.character(alpha)))
  Type1 <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
                                          "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                        as.character(alpha)))
  Prec <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
                                         "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                       as.character(alpha)))
  for(k in 1:n){
    outs <- array(NA, c(9,2,5), dimnames=list(c("aic","par.joint.aic","par.aic",
                                                "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                              c("TP","FP"),models))
    for(i in 1:5){
      outs[1,1:2,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="aic")[1:2])
      outs[2,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.aic")[1:2])
      outs[3,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.aic")[1:2])
      outs[4,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.aic")[1:2])
      outs[5,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="bic")[1:2])
      outs[6,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.bic")[1:2])
      outs[7,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.bic")[1:2])
      outs[8,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.bic")[1:2])
      outs[9,,i] <- as.numeric(performance.summaries.cit(out[[i]], model=models[i], alpha[k])[1:2])
    }
    all <- matrix(NA, 9, 2, dimnames=list(c("aic","par.joint.aic","par.aic",
                                            "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
                                          c("TP","FP")))
    for(i in 1:9)
      all[i,] <- apply(outs[i,1:2,-3],1,sum)
    Prec[,k] <- all[,1]/apply(all,1,sum)
    Power[,k] <- all[,1]/4000
    Type1[,k] <- all[,2]/4000
    print(k)
  }
  list(Power=Power, Type1=Type1, Prec=Prec)
}
