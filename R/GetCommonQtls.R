#' Get common QTLs for phenotypes
#' 
#' Perform joint QTL mapping for phenotypes with marginal LOD peak positions
#' higher than LOD threshold and within set distance of each other
#' 
#' 
#' @param cross object of class `cross`
#' @param pheno1 first phenotype column number or character string name
#' @param pheno2 second phenotype column number or character string name; if
#' more than one, then all phenotypes will be tested against `pheno1`
#' @param thr LOD threshold
#' @param peak.dist maximal peak distance to be considered the same peak (in
#' cM)
#' @param addcov1,addcov2 additive covariates for first and second phenotype,
#' respectively
#' @param intcov1,intcov2 interactive covariates for first and second
#' phenotype, respectively
#' @references Chaibub Neto E, Broman AT, Keller MP, Attie AD, Zhang B, Zhu J,
#' Yandell BS, Causal model selection hypothesis tests in systems genetics.
#' Genetics (in review).
#' @keywords utilities
#' @examples
#' \dontrun{
#' # Create CMSTCross object
#' example(SimCrossCausal)
#' # data(CMSTCross) loaded lazily
#' commqtls <- GetCommonQtls(CMSTCross, 
#'                           pheno1 = "y1", 
#'                           pheno2 = "y3",
#'                           thr = 3,
#'                           peak.dist = 5,
#'                           addcov1 = NULL, 
#'                           addcov2 = NULL, 
#'                           intcov1 = NULL, 
#'                           intcov2 = NULL)
#' commqtls
#' }
#' @export
GetCommonQtls <- function(cross, 
                          pheno1, 
                          pheno2,
                          thr = 3,
                          peak.dist = 5,
                          addcov1 = NULL, 
                          addcov2 = NULL, 
                          intcov1 = NULL, 
                          intcov2 = NULL)
{
  if(length(pheno2) > 1 | length(pheno1) > 1)
    stop("pheno1 and pheno2 must have only one trait")
  
  CreateCovMatrix <- function(cross, cov.names) {
    # get covariates data
    if (!is.null(cov.names)) {
      myformula <- stats::formula(paste("~", paste(cov.names, collapse = "+")))
      out <- stats::model.matrix(myformula, cross$pheno)[, -1]
    }
    else {
      out <- NULL
    }
    out
  }
  FindCommonQtls <- function(cross, scanJ, scan1, scan2, thr, peak.dist) {
    # if multiple common QTLs, returns the strongest one from scanJ
    Q <- NA
    ssJ <- summary(scanJ)
    markers1 <- row.names(summary(scan1, thr))
    markers2 <- row.names(summary(scan2, thr))
    chr1 <- scan1[markers1, 1]
    chr2 <- scan2[markers2, 1]
    match.chr <- chr2[match(chr1, chr2, nomatch = 0)]
    if(length(match.chr) > 0) {
      aux1 <- match(match.chr, chr1)
      aux2 <- match(match.chr, chr2)
      markers1 <- markers1[aux1]
      markers2 <- markers2[aux2]
      peak1 <- scan1[markers1, 2]
      peak2 <- scan2[markers2, 2]
      aux3 <- abs(peak1 - peak2)
      aux4 <- which(aux3 <= peak.dist)
      if(length(aux4) > 0){
        cchr <- match.chr[aux4]
        aux5 <- match(cchr, ssJ[, 1])
        aux6 <- which.max(ssJ[aux5, 3])
        Q.chr <- as.numeric(ssJ[aux5[aux6], 1])
        Q.pos <- ssJ[aux5[aux6], 2]
        Q <- qtl::find.pseudomarker(cross, Q.chr, Q.pos, "prob")
        Q <- data.frame(Q, Q.chr, Q.pos, stringsAsFactors = FALSE)
      }
    }
    Q
  }
  
  to.drop <- DropMissing(cross, c(pheno1, pheno2, addcov1, addcov2,
                                  intcov1, intcov2))
  if (!is.null(to.drop)) {
    cross <- subset(cross, ind = -to.drop)
  }
  n <- qtl::nind(cross)
  y1 <- cross$pheno[, pheno1]
  y2 <- cross$pheno[, pheno2]
  addcov1.M <- CreateCovMatrix(cross, cov.names = addcov1)
  addcov2.M <- CreateCovMatrix(cross, cov.names = addcov2)
  intcov1.M <- CreateCovMatrix(cross, cov.names = intcov1)
  intcov2.M <- CreateCovMatrix(cross, cov.names = intcov2)
  scan1 <- qtl::scanone(cross, pheno.col = qtl::find.pheno(cross, pheno1), 
                        method = "hk", intcovar = intcov1.M,
                        addcovar = cbind(addcov1.M, intcov1.M))
  scan2 <- qtl::scanone(cross, pheno.col = qtl::find.pheno(cross, pheno2), 
                        method = "hk", intcovar = intcov2.M, 
                        addcovar = cbind(addcov2.M, intcov2.M))
  scan2g1 <- qtl::scanone(cross, pheno.col = qtl::find.pheno(cross, pheno2), 
                          method = "hk", intcovar = intcov2.M,
                          addcovar = cbind(addcov2.M, intcov2.M, y1))
  scanJ <- scan1
  scanJ[, 3] <- scan1[, 3] + scan2g1[, 3]
  FindCommonQtls(cross, scanJ, scan1, scan2, thr, peak.dist)
}
