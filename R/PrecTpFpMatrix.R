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
