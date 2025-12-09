#' Wrapper routine for simulations.
#' 
#' Simulate `nSim` realizations of cross object with `n.pheno`
#' phenotypes with correlation `latent.eff`. All simulations use the same
#' genotypes in the `cross` object.
#' 
#' @param nSim Number of simulated sets of phenotypes to create. See details.
#' @param cross Object of class `cross`. See
#' \code{\link[qtl]{read.cross}}.
#' @param n.pheno Number of traits, or phenotypes, to simulate for cross
#' object.
#' @param latent.eff Strength of latent effect, which is included in all
#' traits. See \code{\link{sim.null.cross}}.
#' @param res.var Residual variance for traits. Should not affect results.
#' @param n.quant maximum size of hotspots examined; ideally large enough to
#' exceed the largest Breitling alpha critical value.
#' @param n.perm Number of permutations to perform per realization. Good idea
#' to do 1000, but this takes time.
#' @param alpha.levels Vector of significance levels.
#' @param lod.thrs Vector of LOD thresholds, typically single-trait permutation
#' thresholds for various significance levels.
#' @param drop.lod Drop in LOD score examined. LODs below this drop from the
#' maximum for a chromosome will not be scored.
#' @param verbose verbose output if \code{TRUE}
#' @return `sim.null.cross` simulates an object of class `cross`.
#' `sim.null.pheno.data` simulates a data frame of phenotypes.
#' `sim.hotspot` uses these other routines to simulate a hotspot,
#' returning an list object.
#' @author Elias Chaibub Neto and Brian S. Yandell
#' @seealso \code{\link{sim.null.cross}}, \code{\link[qtl]{read.cross}}.
#' @keywords utilities
#' @examples
#' 
#' ncross1 <- sim.null.cross(chr.len = rep(100, 4),
#'                           n.mar = 51,
#'                           n.ind = 100,
#'                           type = "bc",
#'                           n.phe = 1000,
#'                           latent.eff = 3,
#'                           res.var = 1,
#'                           init.seed = 123457)
#' cross1 <- include.hotspots(cross = ncross1,
#'                            hchr = c(2, 3, 4),
#'                            hpos = c(25, 75, 50),
#'                            hsize = c(100, 50, 20),
#'                            Q.eff = 2,
#'                            latent.eff = 3,
#'                            lod.range.1 = c(2.5, 2.5),
#'                            lod.range.2 = c(5, 8),
#'                            lod.range.3 = c(10, 15),
#'                            res.var = 1,
#'                            n.phe = 1000,
#'                            init.seed = 12345)
#' 
#' @importFrom qtl calc.genoprob nchr nind pull.map sim.cross sim.map
#' @importFrom stats rnorm
sim.hotspot <- function(nSim, 
                        cross, 
                        n.pheno,
                        latent.eff,
                        res.var = 1,
                        n.quant,
                        n.perm,
                        alpha.levels,
                        lod.thrs,
                        drop.lod=1.5,
                        verbose = FALSE)
{
  s.quant <- seq(n.quant)
  
  nalpha <- length(alpha.levels)
  nlod <- length(lod.thrs)
  
  ## outputs count the number of times we detected
  ## a hotspot using the respective method
  outNL <- matrix(0, n.quant, nalpha)
  outN <- outWW <- matrix(0, nlod, nalpha)
  
  ## we are saving the thresholds of each simulation
  thrNL <- array(dim=c(n.quant, nalpha, nSim))
  thrN <- array(dim=c(nlod, nalpha, nSim))
  thrWW <- array(dim=c(nlod, nalpha, nSim))
  
  for(k in 1:nSim){
    mycat(k, verbose, TRUE)
    
    mycat("sim.null.pheno.data", verbose)
    ncross <- sim.null.pheno.data(cross, n.pheno, latent.eff, res.var)
    
    ## Simulate correlated phenotypes and create threshold summaries.
    out.sim <- filter.threshold(ncross, seq(n.pheno), latent.eff[k], res.var,
                                lod.thrs, drop.lod,
                                s.quant, n.perm, alpha.levels,
                                verbose)
    
    thrNL[,,k] <- out.sim$NL.thrs
    thrN[,,k] <- out.sim$N.thrs
    thrWW[,,k] <- out.sim$WW.thrs    
    outNL <- outNL + out.sim$NL
    outN <- outN + out.sim$N.counts
    outWW <- outWW + out.sim$WW.counts
  }
  
  
  NL.err <- outNL/nSim
  dimnames(NL.err) <- list(as.factor(s.quant), as.factor(alpha.levels))
  N.err <- outN / nSim
  dimnames(N.err) <- list(as.factor(lod.thrs), as.factor(alpha.levels))
  WW.err <- outWW / nSim
  dimnames(WW.err) <- list(as.factor(lod.thrs), as.factor(alpha.levels))
  list(nSim = nSim, NL.err=NL.err, N.err=N.err, WW.err=WW.err, thrNL=thrNL, thrN=thrN, 
       thrWW=thrWW)  
}
