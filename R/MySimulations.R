#' @importFrom stats cor runif
#' @importFrom qtl calc.genoprob find.marker pull.geno sim.cross
mySimulations <- function(
  nSim,
  model,
  n.ind,
  mu.range=c(0,0),
  add.eff1.range=c(0,0),
  dom.eff1.range=c(0,0),
  add.eff2.range=c(0,0),
  dom.eff2.range=c(0,0),
  beta21.range=c(0,0),
  beta1h.range=c(0,0),
  beta2h.range=c(0,0),
  sig2.1.range=c(1,1),
  sig2.2.range=c(1,1),
  sig2.h.range=c(1,1),
  eq.spacing=FALSE,
  cross.type="f2",
  thr,
  peak.dist=2,
  normalize=FALSE)
{
  cor12 <- matrix(NA,nSim,1, 
    dimnames=list(c(1:nSim),c("cor.Y1.Y2")))
  R2s <- matrix(NA,nSim,2,
    dimnames=list(c(1:nSim),c("R2.Y1~Q","R2.Y2~Q")))
  BICs <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("BIC.1","BIC.2","BIC.3","BIC.4")))
  AICs <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("AIC.1","AIC.2","AIC.3","AIC.4")))
  z.scores <- matrix(NA,nSim,6,
    dimnames=list(c(1:nSim),c("z.12","z.13","z.14","z.23","z.24","z.34")))
  sig2s <- matrix(NA,nSim,6,
    dimnames=list(c(1:nSim),c("s.12.12","s.13.13","s.14.14","s.23.23","s.24.24","s.34.34")))
  pval.par.cmst.joint.BIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.par.cmst.iu.BIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.non.par.cmst.iu.BIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.par.cmst.joint.AIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.par.cmst.iu.AIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.non.par.cmst.iu.AIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.cit <- matrix(NA,nSim,2,
    dimnames=list(c(1:nSim),c("pval.1","pval.2")))

  if(model=="A"){
    k <- 1
    while(k <= nSim){
      mu <- stats::runif(1,mu.range[1],mu.range[2])
      beta21 <- stats::runif(1,beta21.range[1],beta21.range[2])
      add.eff1 <- stats::runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- stats::runif(1,dom.eff1.range[1],dom.eff1.range[2])
      sig2.1 <- stats::runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- stats::runif(1,sig2.2.range[1],sig2.2.range[2])
      Cross <- qtl::sim.cross(n.ind, mu, beta21, add.eff1, dom.eff1, sig2.1, 
        sig2.2, eq.spacing, cross.type, normalize)
      Cross <- qtl::calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- stats::cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(CMSTtests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(!inherits(aux, "try-error")) {
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- qtl::find.marker(Cross,cq[1,2],cq[1,3])
        LL <- qtl::pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(CitTests(LL, GG, TT),silent=TRUE)
        if(!inherits(aux2, "try-error"))
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  if(model=="B"){
    k <- 1
    while(k <= nSim){
      mu <- stats::runif(1,mu.range[1],mu.range[2])
      beta21 <- stats::runif(1,beta21.range[1],beta21.range[2])
      beta1h <- stats::runif(1,beta1h.range[1],beta1h.range[2])
      beta2h <- stats::runif(1,beta2h.range[1],beta2h.range[2])
      add.eff1 <- stats::runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- stats::runif(1,dom.eff1.range[1],dom.eff1.range[2])
      sig2.1 <- stats::runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- stats::runif(1,sig2.2.range[1],sig2.2.range[2])
      sig2.h <- stats::runif(1,sig2.h.range[1],sig2.h.range[2])
      Cross <- qtl::sim.cross(n.ind, mu, beta21, beta1h, beta2h, add.eff1, 
        dom.eff1, sig2.1, sig2.2, sig2.h, eq.spacing, cross.type, normalize)
      Cross <- qtl::calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- stats::cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(CMSTtests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(!inherits(aux, "try-error")) {
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- qtl::find.marker(Cross,cq[1,2],cq[1,3])
        LL <- qtl::pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(CitTests(LL, GG, TT),silent=TRUE)
        if(!inherits(aux2, "try-error"))
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  if(model=="C"){
    k <- 1
    while(k <= nSim){
      mu <- stats::runif(1,mu.range[1],mu.range[2])
      beta21 <- stats::runif(1,beta21.range[1],beta21.range[2])
      add.eff1 <- stats::runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- stats::runif(1,dom.eff1.range[1],dom.eff1.range[2])
      add.eff2 <- stats::runif(1,add.eff2.range[1],add.eff2.range[2])
      dom.eff2 <- stats::runif(1,dom.eff2.range[1],dom.eff2.range[2])
      sig2.1 <- stats::runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- stats::runif(1,sig2.2.range[1],sig2.2.range[2])
      Cross <- qtl::sim.cross(n.ind, mu, beta21, add.eff1, dom.eff1, add.eff2, 
        dom.eff2, sig2.1, sig2.2, eq.spacing, cross.type, normalize)
      Cross <- qtl::calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- stats::cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(CMSTtests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(!inherits(aux, "try-error")) {
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- qtl::find.marker(Cross,cq[1,2],cq[1,3])
        LL <- qtl::pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(CitTests(LL, GG, TT),silent=TRUE)
        if(!inherits(aux2, "try-error"))
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  if(model=="D"){
    k <- 1
    while(k <= nSim){
      mu <- stats::runif(1,mu.range[1],mu.range[2])
      add.eff1 <- stats::runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- stats::runif(1,dom.eff1.range[1],dom.eff1.range[2])
      add.eff2 <- stats::runif(1,add.eff2.range[1],add.eff2.range[2])
      dom.eff2 <- stats::runif(1,dom.eff2.range[1],dom.eff2.range[2])
      sig2.1 <- stats::runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- stats::runif(1,sig2.2.range[1],sig2.2.range[2])
      Cross <- qtl::sim.cross(n.ind, mu, add.eff1, dom.eff1, add.eff2, 
        dom.eff2, sig2.1, sig2.2, eq.spacing, cross.type, normalize)
      Cross <- qtl::calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- stats::cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(CMSTtests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(!inherits(aux, "try-error")) {
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- qtl::find.marker(Cross,cq[1,2],cq[1,3])
        LL <- qtl::pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(CitTests(LL, GG, TT),silent=TRUE)
        if(!inherits(aux2, "try-error"))
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  if(model=="E"){
    k <- 1
    while(k <= nSim){
      mu <- stats::runif(1,mu.range[1],mu.range[2])
      add.eff1 <- stats::runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- stats::runif(1,dom.eff1.range[1],dom.eff1.range[2])
      add.eff2 <- stats::runif(1,add.eff2.range[1],add.eff2.range[2])
      dom.eff2 <- stats::runif(1,dom.eff2.range[1],dom.eff2.range[2])
      beta1h <- stats::runif(1,beta1h.range[1],beta1h.range[2])
      beta2h <- stats::runif(1,beta2h.range[1],beta2h.range[2])
      sig2.1 <- stats::runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- stats::runif(1,sig2.2.range[1],sig2.2.range[2])
      sig2.h <- stats::runif(1,sig2.h.range[1],sig2.h.range[2])
      Cross <- qtl::sim.cross(n.ind, mu, add.eff1, dom.eff1, add.eff2, dom.eff2, 
        beta1h, beta2h, sig2.1, sig2.2, sig2.h, eq.spacing, cross.type,
        normalize)
      Cross <- qtl::calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- stats::cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(CMSTtests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(!inherits(aux, "try-error")) {
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- qtl::find.marker(Cross,cq[1,2],cq[1,3])
        LL <- qtl::pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(CitTests(LL, GG, TT),silent=TRUE)
        if(!inherits(aux2, "try-error"))
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  list(cor12=cor12,
       R2s=R2s,
       BICs=BICs,
       AICs=AICs,
       z.scores=z.scores,
       sig2s=sig2s,
       pval.par.cmst.joint.BIC=pval.par.cmst.joint.BIC,
       pval.par.cmst.iu.BIC=pval.par.cmst.iu.BIC,
       pval.non.par.cmst.iu.BIC=pval.non.par.cmst.iu.BIC,
       pval.par.cmst.joint.AIC=pval.par.cmst.joint.AIC,
       pval.par.cmst.iu.AIC=pval.par.cmst.iu.AIC,
       pval.non.par.cmst.iu.AIC=pval.non.par.cmst.iu.AIC,
       pval.cit=pval.cit)
}
