#' Determine false positive and true positive rates for known targets.
#' 
#' Determine how well different tests do to predict candidates of regulation.
#' 
#' `FitAllTests()` invokes 7 tests. The hidden routine `CitTests()` is
#' invoked by call to `FitAllTests()`; this is hidden because we do not
#' recommend its use.
#' 
#' `JoinTestOutputs()` joins results of `FitAllTests()`, either
#' from a list `tests` or from a collection of files prefixed by
#' `file`. The joined tests from `JoinTestOutputs` are
#' summarized with `PrecTpFpMatrix()` using the biologically validated true
#' positives, false positives and precision, for the inferred causal relations.
#' We define a true positive as a statistically significant causal relation
#' between a gene and a putative target gene when the putative target gene
#' belongs to the known signature of the gene. Similarly, we define a false
#' positive as a statistically significant causal relation between a gene and a
#' putative target gene when the target gene does not belong to the signature.
#' (For the AIC and BIC methods that do not provide a p-value measuring the
#' significance of the causal call, we simply use the detected causal relations
#' in the computation of true and false positives). The validated precision is
#' computed as the ratio of true positives by the sum of true and false
#' positives. The `PrecTpFpMatrix` computes these measures to both all
#' genes, and to cis genes only. Simulations suggest only non-parametric tests
#' need to be adjusted using Benjamini-Hochberg via `p.adjust.np()`.
#' 
#' @param cross object of class `cross`
#' @param pheno1 first phenotype column number or character string name
#' @param pheno2 second phenotype column number or character string name; if
#' more than one, then all phenotypes will be tested against `pheno1`
#' @param Q.chr QTL chromosome (number or label)
#' @param Q.pos QTL position in cM
#' @param verbose verbose printout if `TRUE
#' @param comap list result of `GetComappingTraits()`
#' @param tests list object as list of `FitAllTests()` results, or of joined
#' output created by `JoinTestsOutputs()`
#' @param file prefix for file names when running `FitAllTests()` in
#' parallel and saving test results in separate files
#' \code{\link[stats]{p.adjust}}
#' @return List containing \item{Prec1,Prec2}{matrix of precision with rows for
#' significance level and columns for test; first is for all, second is for cis
#' candidates only} \item{Tp1,Tp2}{matrix of true positive rate with rows for
#' significance level and columns for test; first is for all, second is for cis
#' candidates only} \item{Fp1,Fp2}{matrix of false positive rate with rows for
#' significance level and columns for test; first is for all, second is for cis
#' candidates only}
#' @author Elias Chaibub Neto
#' @seealso \code{\link[stats]{p.adjust}}
#' @keywords utilities
#' @examples
#' \dontrun{
#' example(GetCandReg)
#' ## Suppose y1 is causal with targets y2 and y3.
#' targets <- list(y1 = c("y2","y3"))
#' 
#' tests <- list()
#' for(k in seq(names(comap.targets))) {
#'   tests[[k]] <- FitAllTests(CMSTCross, pheno1 = names(comap.targets)[k],
#'                       pheno2 = comap.targets[[k]],
#'                       Q.chr = cand.reg[k, 4],
#'                       Q.pos = cand.reg[k, 5])
#' }
#' names(tests) <- names(comap.targets)
#' tests <- JoinTestOutputs(comap.targets, tests)
#' 
#' PrecTpFpMatrix(alpha = seq(0.01, 0.10, by = 0.01),
#'   val.targets = targets, all.orfs = CMSThigh$names, tests = tests,
#'   cand.reg = cand.reg, cis.cand.reg = cis.cand.reg)
#' }
#' @export
FitAllTests <- function(cross, pheno1, pheno2, Q.chr, Q.pos, verbose = TRUE)
{
  out <- CMSTtests(cross, pheno1, pheno2, Q.chr, Q.pos, 
                   NULL, NULL, NULL, NULL, "all", "both", verbose)
  
  nms <- pheno2
  ntests <- length(pheno2)
  out$pvals.cit <- matrix(NA, ntests, 2, dimnames = list(nms, c("pval.1", "pval.2")))
  
  for(k in 1 : ntests) {
    cit.mar <- qtl::find.marker(cross, Q.chr, Q.pos)
    LL <- qtl::pull.geno(cross)[, cit.mar]
    GG <- cross$pheno[, pheno1]
    TT <- cross$pheno[, pheno2[k]]
    aux2 <- try(CitTests(LL, GG, TT), silent = TRUE)
    if(!inherits(aux2, "try-error")) {
      out$pvals.cit[k,] <- aux2
    }
    if(verbose)
      cat("CIT pheno2 = ", pheno2[k], "\n")   
  }
  out
}
#' @param LL,GG,TT Numeric vectors corresponding to genotype
#' @export
#' @rdname FitAllTests
CitTests <- function(LL, GG, TT)
{
  no.bootstrap <- 50
  ### remove missing values ###
  sel <- (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
  dat_f <- as.data.frame(cbind(LL, GG, TT), stringsAsFactors = FALSE)
  dat_f <- dat_f[sel,]
  names(dat_f) <- c("L", "G", "T")
  Lf <- as.factor(dat_f$L)
  dat_f$L <- as.integer(Lf) - 1
  llevels <- as.integer(unique(as.factor(dat_f$L)))
  dfL <- length(llevels) - 1
  pvec <- rep(NA, 4)
  
  if(dfL == 2){
    dat_f$L1 <- ifelse(dat_f$L == 1,1,0)
    dat_f$L2 <- ifelse(dat_f$L == 2,1,0)
    fit0 <- stats::lm(T ~ 1, data = dat_f)
    fit1 <- stats::lm(T ~ L1 + L2, data = dat_f)
    fit2 <- stats::lm(G ~ T, data = dat_f)
    fit3 <- stats::lm(T ~ G, data = dat_f)
    fit4 <- stats::lm(G ~ T + L1 + L2, data = dat_f)
    fit5 <- stats::lm(T ~ G + L1 + L2, data = dat_f)
    pvec[1] <- stats::anova(fit0, fit1)$"Pr(>F)"[2]
    pvec[2] <- stats::anova(fit2, fit4)$"Pr(>F)"[2]
    pvec[3] <- stats::anova(fit1, fit5)$"Pr(>F)"[2]
    f_ <- stats::anova(fit3, fit5)$F[2]
    fit1G <- stats::lm(G ~ L1 + L2, data = dat_f)
    alg <- summary(fit1G)$coefficients["(Intercept)", 1]
    blg1 <- summary(fit1G)$coefficients["L1", 1]
    blg2 <- summary(fit1G)$coefficients["L2", 1]
    alt <- summary(fit1)$coefficients["(Intercept)", 1]
    blt1 <- summary(fit1)$coefficients["L1", 1]
    blt2 <- summary(fit1)$coefficients["L2", 1]
    dat_f$eG <- stats::resid(fit1G)
    dat_f$eT <- stats::resid(fit1)
    ss <- dim(dat_f)[1]
    fvecr <- rep(NA, no.bootstrap)
    fvecr_r <- rep(NA, no.bootstrap)
    for (i in 1 : no.bootstrap) {
      nni <- trunc(1 + ss * stats::runif(ss, 0, 1)) ;
      dat_f$G_ <- alg + blg1 * dat_f$L1 + blg2 * dat_f$L2 + dat_f$eG[nni]
      fit_0 <- stats::lm(T ~ G_, data = dat_f)
      fit_1 <- stats::lm(T ~ G_ + L1 + L2, data = dat_f)
      fvecr[i] <- stats::anova(fit_0, fit_1)$F[2]
      dat_f$T_ <- alt + blt1 * dat_f$L1 + blt2 * dat_f$L2 + dat_f$eT[nni]
      fit_0 <- stats::lm(G ~ T_, data = dat_f)
      fit_1 <- stats::lm(G ~ T_ + L1 + L2, data = dat_f)
      fvecr_r[i] <- stats::anova(fit_0, fit_1)$F[2]
    }
  }#End dfL == 2
  
  if(dfL == 1){
    dat_f$L1 <- ifelse(dat_f$L == 1, 1, 0)
    fit0 <- stats::lm(T ~ 1, data = dat_f)
    fit1 <- stats::lm(T ~ L1, data = dat_f)
    fit2 <- stats::lm(G ~ T, data = dat_f)
    fit3 <- stats::lm(T ~ G, data = dat_f)
    fit4 <- stats::lm(G ~ T + L1, data = dat_f)
    fit5 <- stats::lm(T ~ G + L1, data = dat_f)
    pvec[1] <- stats::anova(fit0, fit1)$"Pr(>F)"[2]
    pvec[2] <- stats::anova(fit2, fit4)$"Pr(>F)"[2]
    pvec[3] <- stats::anova(fit1, fit5)$"Pr(>F)"[2]
    f_ <- stats::anova(fit3, fit5)$F[2]
    fit1G <- stats::lm(G ~ L1, data = dat_f)
    alt <- summary(fit1)$coefficients["(Intercept)", 1]
    blt1 <- summary(fit1)$coefficients["L1", 1]
    alg <- summary(fit1G)$coefficients["(Intercept)", 1]
    blg1 <- summary(fit1G)$coefficients["L1", 1]
    dat_f$eG <- stats::resid(fit1G)
    dat_f$eT <- stats::resid(fit1)
    ss <- dim(dat_f)[1]
    fvecr <- rep(NA, no.bootstrap)
    fvecr_r <- rep(NA, no.bootstrap)
    for (i in 1 : no.bootstrap) {
      nni <- trunc(1 + ss * stats::runif(ss, 0, 1)) 
      dat_f$G_ <- alg + blg1 * dat_f$L1 + dat_f$eG[nni]
      fit_0 <- stats::lm(T ~ G_, data = dat_f)
      fit_1 <- stats::lm(T ~ G_ + L1, data = dat_f)
      fvecr[i] <- stats::anova(fit_0, fit_1)$F[2]
      dat_f$T_ <- alt + blt1 * dat_f$L1 + dat_f$eT[nni]
      fit_0 <- stats::lm(G ~ T_, data = dat_f)
      fit_1 <- stats::lm(G ~ T_ + L1, data = dat_f)
      fvecr_r[i] <- stats::anova(fit_0, fit_1)$F[2]
    }
  } #End dfL == 1
  
  #####F Method
  fvecr <- fvecr[!is.na(fvecr)]
  df1 <- stats::anova(fit3, fit5)$Df[2]
  df2 <- stats::anova(fit3, fit5)$Res.Df[2]
  fncp <- mean(fvecr, na.rm = TRUE) * (df1/df2) * (df2 - df1) - df1
  if(fncp < 0) fncp <- 0
  ######### Transform F to normal
  npvals <- stats::pf(fvecr, df1, df2, ncp = fncp, lower.tail = TRUE)
  nfvecr <- stats::qnorm(npvals)
  npf <- stats::pf(f_, df1, df2, ncp = fncp, lower.tail = TRUE) #Transform observed F
  zf <- stats::qnorm(npf)
  pvec[4] <- stats::pnorm(zf, mean = 0, sd = stats::sd(nfvecr))
  pvalc <- max(pvec)  ###Causal p-value
  
  #### Reactive p-value
  fit0G <- stats::lm(G ~ 1, data = dat_f)
  pvec1 <- rep(NA, 4)
  pvec1[1] <- stats::anova(fit0G, fit1G)$"Pr(>F)"[2]
  pvec1[2] <- stats::anova(fit3, fit5)$"Pr(>F)"[2]
  pvec1[3] <- stats::anova(fit1G, fit4)$"Pr(>F)"[2]
  f_ <- stats::anova(fit2, fit4)$F[2]
  #####F Method
  fvecr_r <- fvecr_r[!is.na(fvecr_r)]
  df1 <- stats::anova(fit3, fit5)$Df[2]
  df2 <- stats::anova(fit3, fit5)$Res.Df[2]
  fncp <- mean(fvecr_r, na.rm = TRUE) * (df1/df2) * (df2 - df1) - df1
  if(fncp < 0) fncp <- 0
  ######### Transform F to normal
  npvals <- stats::pf(fvecr_r, df1, df2, ncp = fncp, lower.tail = TRUE)
  nfvecr <- stats::qnorm(npvals)
  npf <- stats::pf(f_, df1, df2, ncp = fncp, lower.tail = TRUE) #Transform observed F
  zf <- stats::qnorm(npf)
  pvec1[4] <- stats::pnorm(zf, mean = 0, sd = stats::sd(nfvecr))
  pvalr <- max(pvec1)  ###Reactive p-value
  ###
  c(pvalc, pvalr)
}
##############################################################################
##############################################################################
#' @export
#' @rdname FitAllTests
CombineTests <- function(comap, file)
{
  reg.nms <- names(comap)
  out <- NULL
  join.out <- list()
  for (k in seq_along(comap)) {
    load(paste(file, reg.nms[k], "Rdata", sep="."))
    join.out[[k]] <- out
  }
  names(join.out) <- reg.nms
  join.out
}
#' @export
#' @rdname FitAllTests
JoinTestOutputs <- function(comap, tests, file = NULL)
{
  if(!is.null(file) && missing(tests)) {
    tests <- CombineTests(comap, file)
  }
  
  reg.nms <- names(comap)
  join.out <- tests[[1]]
  ## Add extra element to join.out: phenos.
  join.out$phenos <- cbind(rep(reg.nms[1], length(comap[[1]])), comap[[1]])
  
  for (k in 2 : length(comap)) {
    out <- tests[[k]]
    if (length(out) != 10) { ## CMSTtests if length(pheno2) == 1
      tmp <- t(out$Z.aic)
      AIC.stats <- c(out$AICs, tmp[!is.na(tmp)])
      tmp <- t(out$Z.bic)
      BIC.stats <- c(out$BICs, tmp[!is.na(tmp)])      
      out <- list(R2s = out$R2, AIC.stats = AIC.stats, BIC.stats = BIC.stats, 
                  pvals.j.BIC = out$pvals.j.BIC, pvals.p.BIC = out$pvals.p.BIC, 
                  pvals.np.BIC = out$pvals.np.BIC,  pvals.j.AIC = out$pvals.j.AIC, 
                  pvals.p.AIC = out$pvals.p.AIC,  pvals.np.AIC = out$pvals.np.AIC,
                  pvals.cit = out$pvals.cit)
    }
    
    for (i in 1 : 10) {
      join.out[[i]] <- rbind(join.out[[i]], out[[i]])
    }
    join.out[[11]] <- 
      rbind(join.out[[11]], cbind(rep(reg.nms[k], length(comap[[k]])), comap[[k]]))
  }
  join.out
}
