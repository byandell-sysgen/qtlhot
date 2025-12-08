#' Generates a "null dataset" cross
#'
#' @param chr.len 
#' @param n.mar 
#' @param n.ind 
#' @param type 
#' @param n.pheno 
#' @param latent.eff 
#' @param res.var 
#' @param init.seed 
#'
#' @export
sim.null.cross <- function(chr.len = rep(400,16), n.mar=185, n.ind = 112, type = "bc",
                           n.pheno = 6000, latent.eff = 1.5, res.var = 1,
                           init.seed = 92387475)
{
  set.seed(init.seed)
  mymap <- qtl::sim.map(len = chr.len, n.mar = n.mar, include.x=FALSE, eq.spacing=TRUE)
  scross <- qtl::sim.cross(map = mymap, n.ind = n.ind, type = type)
  scross <- qtl::calc.genoprob(scross, step=0)

  sim.null.pheno.data(cross = scross, n.pheno = n.pheno, latent.eff = latent.eff, res.var = res.var)
}
sim.null.pheno.data <- function(cross, n.pheno, latent.eff, res.var)
{
  n <- qtl::nind(cross)
  latent <- stats::rnorm(n, 0, sqrt(res.var))
  ErrorM <- matrix(stats::rnorm(n * n.pheno, 0, sqrt(res.var)), n, n.pheno)
  pheno <- data.frame(latent*latent.eff + ErrorM)
  names(pheno) <- paste("P", 1:n.pheno, sep="")
  cross$pheno <- pheno
  
  cross
}

#' Generate hotspots for the simulated examples in the manuscript.
#'
#' @param cross 
#' @param hchr 
#' @param hpos 
#' @param hsize 
#' @param Q.eff 
#' @param latent.eff 
#' @param lod.range.1 
#' @param lod.range.2 
#' @param lod.range.3 
#' @param res.var 
#' @param n.pheno 
#' @param init.seed 
#' 
#' @export
include.hotspots <- function(cross,
                             hchr,
                             hpos,
                             hsize,
                             Q.eff,
                             latent.eff,
                             lod.range.1,
                             lod.range.2,
                             lod.range.3,
                             res.var=1,
                             n.pheno,
                             init.seed)
{
  get.closest.pos.nms <- function(pos, cross, chr)
  {
    ## map <- attributes(cross$geno[[chr]]$prob)$map
    ## This can be simplified.
    map <- attributes(cross$geno[[chr]]$prob)$map
    ## map <- qtl::pull.map(cross, chr)
    q.nms <- names(map)
    map.pos <- as.numeric(map)
    tmp <- which.min(abs(map.pos-pos))
    closest.pos <- map.pos[tmp]
    nms <- q.nms[tmp]
    list(nms,closest.pos)
  }
  pull.prob <- function(cross)
  {
    out <- vector(mode = "list", length = qtl::nchr(cross))
    names(out) <- names(cross$geno)
    for(i in names(out))
      out[[i]] <- cross$geno[[i]]$prob  
    out
  }
  get.qtl.eff <- function(n, lod, res.var, latent.eff, Q.eff)
  {
##    lod <- stats::runif(hsize, lod.range[1], lod.range[2])
##    r2 <- 1 - 10 ^ (-2 * lod / qtl::nind(cross))
##    beta <- sqrt(r2 * res.var * (1 + latent.eff ^ 2) / (Q.eff ^ 2 * (1 - r2) - res.var * r2))
    r2 <- 1 - 10^(-2*lod/n)
    sqrt(r2*(1 + latent.eff^2)/(Q.eff^2*(1 - r2) - r2))
  }
  update.pheno <- function(cross, hchr, hpos, hsize, Q.eff, latent.eff, lod.range,
                           res.var, index, hk.prob)
  {
    M.pos <- get.closest.pos.nms(hpos, cross, hchr)[[1]]
    M.dummy <- hk.prob[[hchr]][, M.pos, 1] - hk.prob[[hchr]][, M.pos, 2]
    M <- M.dummy * Q.eff + stats::rnorm(qtl::nind(cross), 0, sqrt(res.var))

    ## QTL effect
    beta <- get.qtl.eff(n = qtl::nind(cross),
                        lod = stats::runif(hsize, lod.range[1], 
                                           lod.range[2]),
                        res.var,
                        latent.eff,
                        Q.eff)

    for(j in seq(length(index)))
      cross$pheno[, index[j]] <- beta[j] * M + cross$pheno[, index[j]]

    cross
  }

  set.seed(init.seed)
  hk.prob <- pull.prob(cross)

  ## Why 50 for first, 500 for 2nd and 3rd?
  ## Why strange lod.range for 2nd?
  index1 <- sample(1:n.pheno, hsize[1], replace = FALSE)
  cross <- update.pheno(cross, hchr[1], hpos[1], hsize[1], Q.eff, latent.eff,
                        lod.range.1, res.var, index1, hk.prob)

  index2 <- sample((1:n.pheno)[-index1], hsize[2], replace = FALSE)
  cross <- update.pheno(cross, hchr[2], hpos[2], hsize[2], Q.eff, latent.eff,
                        lod.range.2, res.var, index2, hk.prob)

  index3 <- sample((1:n.pheno)[-c(index1, index2)], hsize[3], replace = FALSE)
  cross <- update.pheno(cross, hchr[3], hpos[3], hsize[3], Q.eff, latent.eff,
                        lod.range.3, res.var, index3, hk.prob)

  cross
}
#################################################################################
mySimulations <- function(...) sim.hotspot(...)
#################################################################################


#' Wrapper routine for simulations.
#' 
#' Wrapper routine for simulations
#' 
#' Simulate \code{nSim} realizations of cross object with \code{n.pheno}
#' phenotypes with correlation \code{latent.eff}. All simulations use the same
#' genotypes in the \code{cross} object.
#' 
#' @param nSim Number of simulated sets of phenotypes to create. See details.
#' @param cross Object of class \code{cross}. See
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
#' @param init.seed initial seed for pseudo-random number generation
#' @param chr.len vector of chromosome lengths
#' @param n.mar number of markers
#' @param n.ind number of individuals
#' @param type type of cross
#' @param hchr,hpos,hsize vectors for hotspot chromosomes, positions, and sizes
#' @param Q.eff QTL effect
#' @param lod.range.1,lod.range.2,lod.range.3 2-vectors of LOD ranges for
#' multiple purposes
#' @param verbose Verbose output if \code{TRUE}. More detailed output if
#' \code{2}.
#' @param \dots Arguments passed directly to \code{sim.hotspot}.
#' @return \code{sim.null.cross} simulates an object of class \code{cross}.
#' \code{sim.null.pheno.data} simulates a data frame of phenotypes.
#' \code{sim.hotspot} uses these other routines to simulate a hotspot,
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
