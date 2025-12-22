#' Perform CMST Tests on cross object
#' 
#' Performs 6 separate CMST tests (3 versions, 2 penalties).
#' 
#' Explain method and penalty here.
#' 
#' @param cross object of class `cross`
#' @param pheno1 first phenotype column number or character string name
#' @param pheno2 second phenotype column number or character string name; if
#' more than one, then all phenotypes will be tested against `pheno1`
#' @param Q.chr QTL chromosome (number or label)
#' @param Q.pos QTL position in cM
#' @param addcov1,addcov2 additive covariates for first and second phenotype,
#' respectively
#' @param intcov1,intcov2 interactive covariates for first and second
#' phenotype, respectively
#' @param method test method; see details
#' @param penalty type of penalty; see details
#' @param verbose verbose printout if `TRUE`
#' @param highobj High LOD object used for CMST tests.
#' @param cand.reg Candidate regions for QTL mapping.
#' @param lod.thr LOD threshold for filtering significant regions.
#' @seealso `CMSTCross`, `PrecTpFpMatrix`,
#' `FitAllTests`
#' @references Chaibub Neto E, Broman AT, Keller MP, Attie AD, Zhang B, Zhu J,
#' Yandell BS, Causal model selection hypothesis tests in systems genetics.
#' Genetics (in review).
#' @keywords utilities
#' @examples
#' \dontrun{
#' # Create CMSTCross object
#' example(SimCrossCausal)
#' nms <- names(CMSTCross$pheno)
#' out1 <- CMSTtests(CMSTCross, 
#'                   pheno1 = nms[1], 
#'                   pheno2 = nms[2],
#'                   Q.chr = 1,
#'                   Q.pos = 55,
#'                   addcov1 = NULL, 
#'                   addcov2 = NULL, 
#'                   intcov1 = NULL, 
#'                   intcov2 = NULL, 
#'                   method = "all",
#'                   penalty = "both")
#' out1[1:6]
#' out1[7]
#' out1[8:12]
#' out1[13:17]
#' ## list of phenotypes
#' out2 <- CMSTtests(CMSTCross, 
#'                   pheno1 = nms[1], 
#'                   pheno2 = nms[-1],
#'                   Q.chr = 1,
#'                   Q.pos = 55,
#'                   addcov1 = NULL, 
#'                   addcov2 = NULL, 
#'                   intcov1 = NULL, 
#'                   intcov2 = NULL, 
#'                   method = "par",
#'                   penalty = "bic")
#' out2
#' }
#' @export
#' @importFrom corpcor is.positive.definite make.positive.definite
#' @importFrom mnormt pmnorm
#' @importFrom qtl find.marker find.pseudomarker lodint makeqtl nind pull.geno
#'             scanone
CMSTtests <- function(cross, 
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
                      verbose = FALSE,
                      highobj = NULL,
                      cand.reg = NULL,
                      lod.thr = NULL)
{
  if (!any(class(cross) == "cross")) 
    stop("Input should have class \"cross\".")

  if(length(pheno2) > 1)
    return(CMSTtestsList(cross, pheno1, pheno2, Q.chr, Q.pos,
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose))
  
  cross.type <- class(cross)[1]

  ## Fit CMST tests
  HkDesignMatrix <- function(qtlo, cross.type) {
    # Create design matrix (Haley-Knott regression)
    nr <- nrow(qtlo$prob[[1]])
    ng <- length(qtlo$prob)
    if(cross.type == "f2"){
      tmp <- 
        unlist(lapply(qtlo$prob, function(x) cbind(x[,1] - x[,3], x[,2])))
      hkm <- matrix(tmp, nr, 2 * ng)
    }
    if(cross.type == "bc"){
      tmp <- unlist(lapply(qtlo$prob, function(x) x[,1] - x[,2]))
      hkm <- matrix(tmp, nr, ng)
    }
    cbind(rep(1, nr), hkm)
  }

  CreateDesignMatrix <- function(cross, chr, pos, addcov.nms, intcov.nms, 
                                 cross.type) {
    # Creates a design matrix
    le.mar <- length(chr)
    n <- qtl::nind(cross)
    if(le.mar > 0) {
      qtlo <- qtl::makeqtl(cross, chr, pos, what = "prob")
      XG <- HkDesignMatrix(qtlo, cross.type)
    }
    else {
      XG <- matrix(1,n,1)    
    }
    if(!is.null(intcov.nms)) {
      intcov.dat <- data.frame(cross$pheno[, intcov.nms])
      names(intcov.dat) <- intcov.nms
      int.sub.matrix <- stats::model.matrix( stats::as.formula(paste("~", intcov.nms)), 
        intcov.dat)[,-1]
      covs <- data.frame(cross$pheno[, unique(c(addcov.nms, intcov.nms))])
      names(covs) <- unique(c(addcov.nms, intcov.nms))
      form <- stats::as.formula(paste(" ~ ", paste(names(covs), collapse="+")))
      X <- as.matrix(stats::model.matrix(form, data=covs)[,-1])
      if(ncol(XG) > 1) {
        genobyintcov <- XG[,-1] * int.sub.matrix
        X <- cbind(XG, X, genobyintcov)
      }
      else {
        X <- cbind(XG,X)
      }
    }
    else {
      if(!is.null(addcov.nms)){
        covs <- data.frame(cross$pheno[, addcov.nms])
        names(covs) <- addcov.nms
        form <- stats::as.formula(paste(" ~ ", paste(names(covs), collapse = "+")))
        X <- stats::model.matrix(form, data = covs)[,-1]
        if(is.null(dim(X))) {
          X <- as.matrix(X)
        }
        X <- cbind(XG, X)
      }
      else {
        X <- XG
      }
    }
    X
  }

  ParametricIUCMST <- function(Z) {
    pv <- matrix(NA, 4, 4)
    for (i in 1 : 3) {
      for(j in (i + 1) : 4) {
        pv[i, j] <- stats::pnorm(Z[i, j], lower.tail = FALSE)
      }
    }
    pval.1 <- max(pv[1, 2], pv[1, 3], pv[1, 4])
    pval.2 <- max(1 - pv[1, 2], pv[2, 3], pv[2, 4])
    pval.3 <- max(1 - pv[1, 3], 1 - pv[2, 3], pv[3, 4])
    pval.4 <- max(1 - pv[1, 4], 1 - pv[2, 4], 1 - pv[3, 4])
    c(pval.1, pval.2, pval.3, pval.4)
  }

  NonparametricIUCMST <- function(penalty, n, k, vec.logLik) {
    # Computes the non-parametric CMST
    # k: vector of length 4 with the model dimensions
    # vec.logLik: matrix of loglik scores (n by 4)
    kM <- matrix(rep(k, each = n), n, 4)
    if (penalty == "bic") {
      vec.penal <- vec.logLik - 0.5 * kM * log(n)/n
    }
    if (penalty == "aic") {
      vec.penal <- vec.logLik - kM/n
    }
    pv <- matrix(NA, 4, 4)
    for (i in 1 : 3) {
      for (j in (i + 1) : 4) {
        vec.ratio <- vec.penal[, i] - vec.penal[, j]
        counts <- sum(vec.ratio > 0)
        nn <- n - sum(vec.ratio == 0)
        pv[i, j] <- stats::pbinom(counts - 1, nn, 0.5, lower.tail = FALSE)
        pv[j, i] <- stats::pbinom(counts, nn, 0.5, lower.tail = TRUE)
      }
    }
    pval.1 <- max(pv[1, 2], pv[1, 3], pv[1, 4])
    pval.2 <- max(pv[2, 1], pv[2, 3], pv[2, 4])
    pval.3 <- max(pv[3, 1], pv[3, 2], pv[3, 4])
    pval.4 <- max(pv[4, 1], pv[4, 2], pv[4, 3])
    c(pval.1, pval.2, pval.3, pval.4)
  }

  ParametricJointCMST <- function(Z, Cor.hat) {
    z <- min(Z[1, 2], Z[1, 3], Z[1, 4])
    pval.1 <- 1 - mnormt::pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[1]])
    z <- min(- Z[1, 2], Z[2, 3], Z[2, 4])
    pval.2 <- 1 - mnormt::pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[2]])
    z <- min(-Z[1, 3], -Z[2, 3], Z[3, 4])
    pval.3 <- 1 - mnormt::pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[3]])
    z <- min(-Z[1, 4], -Z[2, 4], -Z[3, 4])
    pval.4 <- 1 - mnormt::pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[4]])
    c(pval.1, pval.2, pval.3, pval.4)
  }

  GetLogLik <- function(cross, y, n, chr, pos, addcov.nms, intcov.nms, 
                        cross.type) {
    X <- CreateDesignMatrix(cross, chr, pos, addcov.nms, intcov.nms, 
                            cross.type)
    dX <- ncol(X)
    qrX <- qr(X)
    b <- qr.coef(qrX, y)
    RSS <- crossprod(y - X %*% b, y - X %*% b)
    log.lik <- as.vector(- (n/2) - (n/2) * log(2 * pi) - (n/2) * log(RSS/n))
    ss <- RSS/n
    vec.log.lik <- stats::dnorm(y, X %*% b, sqrt(ss), log = TRUE)
    list(log.lik = log.lik, vec.log.lik = vec.log.lik, d = dX, RSS = RSS)
  }

  method <- match.arg(method)
  penalty <- match.arg(penalty)

  to.drop <- DropMissing(cross, c(pheno1, pheno2, addcov1, addcov2, 
                         intcov1, intcov2))
  
  if (!is.null(to.drop)) {
    cross <- subset(cross, ind = -to.drop)
  }
  n <- qtl::nind(cross)
  y1 <- cross$pheno[, pheno1]
  y2 <- cross$pheno[, pheno2]

  # log.lik.1 #
  tmp.1 <- GetLogLik(cross, y1, n, Q.chr, Q.pos, addcov1, intcov1, 
                     cross.type)
  tmp.2 <- GetLogLik(cross, y2, n, Q.chr, Q.pos, addcov2, intcov2, 
                     cross.type)
  tmp.1g2 <- GetLogLik(cross, y1, n, NULL, NULL, c(addcov1, pheno2), 
                       intcov1, cross.type)
  tmp.2g1 <- GetLogLik(cross, y2, n, NULL, NULL, c(addcov2, pheno1), 
                       intcov2, cross.type)
  tmp.2g1.j <- GetLogLik(cross, y2, n, Q.chr, Q.pos, c(addcov2, pheno1), 
                         intcov2, cross.type)

  TSS1 <- sum((y1 - mean(y1))^2)
  TSS2 <- sum((y2 - mean(y2))^2)
  R2 <- c(1 - (tmp.1$RSS/TSS1), 1 - (tmp.2$RSS/TSS2))

  loglik <- c(tmp.1$log.lik + tmp.2g1$log.lik,
              tmp.2$log.lik + tmp.1g2$log.lik,
              tmp.2$log.lik + tmp.1$log.lik,
              tmp.1$log.lik + tmp.2g1.j$log.lik)  

  model.dim <- 2 + c(tmp.1$d + tmp.2g1$d,
                     tmp.2$d + tmp.1g2$d,
                     tmp.1$d + tmp.2$d,
                     tmp.1$d + tmp.2g1.j$d)

  BICs <- -2 * loglik + model.dim * log(n)
  AICs <- -2 * loglik + 2 * model.dim

  vec.logLik <- cbind(tmp.1$vec.log.lik + tmp.2g1$vec.log.lik,
                      tmp.2$vec.log.lik + tmp.1g2$vec.log.lik,
                      tmp.2$vec.log.lik + tmp.1$vec.log.lik,
                      tmp.1$vec.log.lik + tmp.2g1.j$vec.log.lik)

  vec.LR <- matrix(NA, n, 6, dimnames = list(NULL, 
                   c("12", "13", "14", "23", "24", "34")))

  vec.LR <- matrix(NA, n, 6)
  # i = 1, 12
  # i = 2, 13
  # i = 3, 14
  # i = 4, 23
  # i = 5, 24
  # i = 6, 34 
  ip <- matrix(c(1, 1, 1, 2, 2, 3, 2, 3, 4, 3, 4, 4), 6, 2)
  for (i in 1 : 6) {
    vec.LR[, i] <- vec.logLik[, ip[i, 1]] - vec.logLik[, ip[i, 2]]
  }
  S.hat <- (1 - 1/n) * stats::cov(vec.LR)
  dimnames(S.hat) <- list(dimnames(vec.LR)[[2]], dimnames(vec.LR)[[2]])

  if (method == "all") {
    Sig.hat <- vector(mode = "list", length = 4) 
    signs.2 <- matrix(c(1, -1, -1, -1, 1, 1, -1, 1, 1), 3, 3)
    signs.3 <- matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 1), 3, 3)
    # 1: 12.13.14
    Sig.hat[[1]] <- S.hat[c(1, 2, 3), c(1, 2, 3)]
    # 2: 21.23.24
    Sig.hat[[2]] <- S.hat[c(1, 4, 5), c(1, 4, 5)] * signs.2
    # 3: 31.32.34
    Sig.hat[[3]] <- S.hat[c(2, 4, 6), c(2, 4, 6)] * signs.3
    # 4: 41.42.43
    Sig.hat[[4]] <- S.hat[c(3, 5, 6), c(3, 5, 6)]
    Cor.hat <- vector(mode = "list", length = 4) 
    for (i in 1 : 4) {
      tmp <- corpcor::is.positive.definite(Sig.hat[[i]])
      if (!tmp) {
        Sig.hat[[i]] <- corpcor::make.positive.definite(Sig.hat[[i]])
      }
      Cor.hat[[i]] <- stats::cov2cor(Sig.hat[[i]])
      attr(Sig.hat[[i]], "is.positive.definite") <- tmp
    }
    if (penalty == "both") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.BIC <- ParametricIUCMST(Z.bic)
      pvals.np.BIC <- NonparametricIUCMST("bic", n, model.dim, vec.logLik)
      pvals.j.BIC <- ParametricJointCMST(Z.bic, Cor.hat)
      pvals.p.AIC <- ParametricIUCMST(Z.aic)
      pvals.np.AIC <- NonparametricIUCMST("aic", n, model.dim, vec.logLik)
      pvals.j.AIC <- ParametricJointCMST(Z.aic, Cor.hat)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.p.BIC = pvals.p.BIC,
                  pvals.np.BIC = pvals.np.BIC,
                  pvals.j.BIC = pvals.j.BIC,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.p.AIC = pvals.p.AIC,
                  pvals.np.AIC = pvals.np.AIC,
                  pvals.j.AIC = pvals.j.AIC)   
    }
    else if (penalty == "bic") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.BIC <- ParametricIUCMST(Z.bic)
      pvals.np.BIC <- NonparametricIUCMST("bic", n, model.dim, vec.logLik)
      pvals.j.BIC <- ParametricJointCMST(Z.bic, Cor.hat)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.p.BIC = pvals.p.BIC,
                  pvals.np.BIC = pvals.np.BIC,
                  pvals.j.BIC = pvals.j.BIC)         
    }
    else if (penalty == "aic") {
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.AIC <- ParametricIUCMST(Z.aic)
      pvals.np.AIC <- NonparametricIUCMST("aic", n, model.dim, vec.logLik)
      pvals.j.AIC <- ParametricJointCMST(Z.aic, Cor.hat) 
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.p.AIC = pvals.p.AIC,
                  pvals.np.AIC = pvals.np.AIC,
                  pvals.j.AIC = pvals.j.AIC)     
    }
  }
  else if (method == "joint") {
    Sig.hat <- vector(mode = "list", length = 4) 
    signs.2 <- matrix(c(1, -1, -1, -1, 1, 1, -1, 1, 1), 3, 3)
    signs.3 <- matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 1), 3, 3)
    # 1: 12.13.14
    Sig.hat[[1]] <- S.hat[c(1, 2, 3), c(1, 2, 3)]
    # 2: 21.23.24
    Sig.hat[[2]] <- S.hat[c(1, 4, 5), c(1, 4, 5)] * signs.2
    # 3: 31.32.34
    Sig.hat[[3]] <- S.hat[c(2, 4, 6), c(2, 4, 6)] * signs.3
    # 4: 41.42.43
    Sig.hat[[4]] <- S.hat[c(3, 5, 6), c(3, 5, 6)]
    Cor.hat <- vector(mode = "list", length = 4) 
    for (i in 1 : 4) {
      tmp <- corpcor::is.positive.definite(Sig.hat[[i]])
      if (!tmp) {
        Sig.hat[[i]] <- corpcor::make.positive.definite(Sig.hat[[i]])
      }
      Cor.hat[[i]] <- stats::cov2cor(Sig.hat[[i]])
      attr(Sig.hat[[i]], "is.positive.definite") <- tmp
    }
    if (penalty == "both") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.j.BIC <- ParametricJointCMST(Z.bic, Cor.hat)
      pvals.j.AIC <- ParametricJointCMST(Z.aic, Cor.hat)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.j.BIC = pvals.j.BIC,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.j.AIC = pvals.j.AIC)   
    }
    else if (penalty == "bic") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.j.BIC <- ParametricJointCMST(Z.bic, Cor.hat)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.j.BIC = pvals.j.BIC)         
    }
    else if (penalty == "aic") {
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.j.AIC <- ParametricJointCMST(Z.aic, Cor.hat) 
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.j.AIC = pvals.j.AIC)     
    }
  }
  else if (method == "par") {
    if (penalty == "both") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.BIC <- ParametricIUCMST(Z.bic)
      pvals.p.AIC <- ParametricIUCMST(Z.aic)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.p.BIC = pvals.p.BIC,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.p.AIC = pvals.p.AIC)
    }
    else if (penalty == "bic") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.BIC <- ParametricIUCMST(Z.bic)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.p.BIC = pvals.p.BIC)     
    }
    else if (penalty == "aic") {
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.AIC <- ParametricIUCMST(Z.aic)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.p.AIC = pvals.p.AIC)
    }
  }
  if (method == "non.par") {
    if (penalty == "both") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.np.BIC <- NonparametricIUCMST("bic", n, model.dim, vec.logLik)
      pvals.np.AIC <- NonparametricIUCMST("aic", n, model.dim, vec.logLik)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.np.BIC = pvals.np.BIC,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.np.AIC = pvals.np.AIC)   
    }
    else if (penalty == "bic") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.np.BIC <- NonparametricIUCMST("bic", n, model.dim, vec.logLik)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.np.BIC = pvals.np.BIC)         
    }
    else if (penalty == "aic") {
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.np.AIC <- NonparametricIUCMST("aic", n, model.dim, vec.logLik)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.np.AIC = pvals.np.AIC)     
    }
  }
  out
}
