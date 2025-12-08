#' @export
#' @rdname CMSTtests
CMSTtestsList <- function(cross, 
                          pheno1, 
                          pheno2,
                          Q.chr,
                          Q.pos,
                          addcov1 = NULL, 
                          addcov2 = NULL, 
                          intcov1 = NULL, 
                          intcov2 = NULL, 
                          method = c("par", "non.par", "joint", "all"),
                          penalty = c("bic", "aic", "both"),
                          verbose = TRUE)
{
  if(length(pheno2) == 1)
    return(CMSTtestsList(cross, pheno1, pheno2, Q.chr, Q.pos,
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose))
  
  cross.type <- class(cross)[1]
  
  ntests <- length(pheno2)
  nms <- paste(pheno1, pheno2, sep = "_")
  pval.nms <- c("pval.1", "pval.2", "pval.3", "pval.4")
  if (penalty == "both") {
    AIC.nms <- c("AIC.1", "AIC.2", "AIC.3", "AIC.4", "z.12", "z.13", "z.14",
                 "z.23", "z.24", "z.34")
    BIC.nms <- c("BIC.1", "BIC.2", "BIC.3", "BIC.4", "z.12", "z.13", "z.14",
                 "z.23", "z.24", "z.34")
    if (method == "all") {
      out <- vector(mode = "list", length = 9)
      names(out) <- c("R2s", "AIC.stats", "BIC.stats", 
                      "pvals.j.BIC", "pvals.p.BIC", "pvals.np.BIC",
                      "pvals.j.AIC", "pvals.p.AIC", "pvals.np.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 4 : 9) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)]) 
        out[[3]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[4]][k,] <- aux$pvals.j.BIC
        out[[5]][k,] <- aux$pvals.p.BIC
        out[[6]][k,] <- aux$pvals.np.BIC
        out[[7]][k,] <- aux$pvals.j.AIC
        out[[8]][k,] <- aux$pvals.p.AIC
        out[[9]][k,] <- aux$pvals.np.AIC
        if(verbose)
          cat("pheno2 = ", k, "\n")   
      }
    }
    else if (method == "par") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "AIC.stats", "BIC.stats", "pvals.p.BIC", 
                      "pvals.p.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 4 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)]) 
        out[[3]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[4]][k,] <- aux$pvals.p.BIC
        out[[5]][k,] <- aux$pvals.p.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "non.par") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "AIC.stats", "BIC.stats", "pvals.np.BIC",
                      "pvals.np.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 4 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)]) 
        out[[3]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[4]][k,] <- aux$pvals.np.BIC
        out[[5]][k,] <- aux$pvals.np.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "joint") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "AIC.stats", "BIC.stats", "pvals.j.BIC",
                      "pvals.j.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 4 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)]) 
        out[[3]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[4]][k,] <- aux$pvals.j.BIC
        out[[5]][k,] <- aux$pvals.j.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
  }
  else if (penalty == "bic") {
    BIC.nms <- c("BIC.1", "BIC.2", "BIC.3", "BIC.4", "z.12", "z.13", "z.14",
                 "z.23", "z.24", "z.34")
    if (method == "all") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "BIC.stats", "pvals.j.BIC", "pvals.p.BIC",
                      "pvals.np.BIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 3 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[3]][k,] <- aux$pvals.j.BIC
        out[[4]][k,] <- aux$pvals.p.BIC
        out[[5]][k,] <- aux$pvals.np.BIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "joint") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "BIC.stats", "pvals.j.BIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[3]][k,] <- aux$pvals.j.BIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "par") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "BIC.stats", "pvals.p.BIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[3]][k,] <- aux$pvals.p.BIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "non.par") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "BIC.stats", "pvals.np.BIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[3]][k,] <- aux$pvals.np.BIC
        cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
  } 
  else if (penalty == "aic") {
    AIC.nms <- c("AIC.1", "AIC.2", "AIC.3", "AIC.4", "z.12", "z.13", "z.14",
                 "z.23", "z.24", "z.34")
    if (method == "all") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "AIC.stats", "pvals.j.AIC", "pvals.p.AIC",
                      "pvals.np.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      for (i in 3 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)])
        out[[3]][k,] <- aux$pvals.j.AIC
        out[[4]][k,] <- aux$pvals.p.AIC
        out[[5]][k,] <- aux$pvals.np.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "joint") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "AIC.stats", "pvals.j.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)])
        out[[3]][k,] <- aux$pvals.j.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "par") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "AIC.stats", "pvals.p.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)])
        out[[3]][k,] <- aux$pvals.p.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "non.par") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "AIC.stats", "pvals.np.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)])
        out[[3]][k,] <- aux$pvals.np.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
  } 
  out
}
