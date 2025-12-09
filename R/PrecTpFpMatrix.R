#' Compute Precision, True Positives, and False Positives
#' 
#' The `PrecTpFpMatrix` function calculates precision, true positives, and false positives for various tests.
#' It is used to evaluate the performance of QTL mapping methods.
#' 
#' @param alpha A numeric vector of significance levels.
#' @param val.targets A vector of validated target genes.
#' @param all.orfs A vector of all open reading frames (ORFs).
#' @param tests A list of test results for different methods.
#' @param cand.reg Candidate regions for QTL mapping.
#' @param cis.cand.reg Cis-regulatory candidate regions.
#' 
#' @return A list of matrices containing precision, true positives, and false positives for each method and significance level.
#' @examples
#' \dontrun{
#' results <- PrecTpFpMatrix(alpha = seq(0.01, 0.1, by = 0.01),
#'                           val.targets = validated_genes,
#'                           all.orfs = all_genes,
#'                           tests = test_results,
#'                           cand.reg = candidate_regions,
#'                           cis.cand.reg = cis_regions)
#' }
#' @export
PrecTpFpMatrix <- function(alpha, val.targets, all.orfs, tests, cand.reg, cis.cand.reg)
{
  nms <- as.character(cand.reg[,1])
  cis.index <- attr(cis.cand.reg, "cis.index")
  
  le <- length(alpha)
  Prec1 <- Tp1 <- Fp1 <- matrix(NA, 9, le, 
                                dimnames=list(c("aic", "bic", "j.bic", "p.bic", "np.bic", "j.aic", 
                                                "p.aic", "np.aic", "cit"), as.character(alpha)))
  Prec2 <- Tp2 <- Fp2 <- Prec1
  for(i in 1:le){
    aux <- PerformanceSummariesKo(alpha = alpha[i], nms,
                                  val.targets = val.targets, 
                                  all.orfs = all.orfs, 
                                  tests = tests,
                                  cis.index = cis.index)
    Prec1[,i] <- round(aux[[1]][,1], 2) 
    Prec2[,i] <- round(aux[[2]][,1], 2)
    Tp1[,i] <- aux[[1]][, 2]
    Tp2[,i] <- aux[[2]][, 2]
    Fp1[,i] <- aux[[1]][, 3]
    Fp2[,i] <- aux[[2]][, 3]
  }
  list(Prec1 = Prec1,
       Prec2 = Prec2,
       Tp1 = Tp1,
       Tp2 = Tp2,
       Fp1 = Fp1,
       Fp2 = Fp2)
}
