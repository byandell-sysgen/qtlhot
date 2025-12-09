#' Generate hotspots for the simulated examples in the manuscript.
#'
#' @param cross cross object
#' @param hchr character string for chromosome
#' @param hpos position on chromosome
#' @param hsize size
#' @param Q.eff Q effect
#' @param latent.eff latent effect
#' @param lod.range.1 lod range 1
#' @param lod.range.2 lod range 2
#' @param lod.range.3 lod range 3
#' @param res.var residual variance
#' @param n.pheno number of phenotypes
#' @param init.seed initial seed
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
