#' Generates a "null dataset" cross
#'
#' @param chr.len vector with length of chromosomes 
#' @param n.mar number of markers
#' @param n.ind number of individuals
#' @param type cross type as character
#' @param n.pheno number of phenotypes
#' @param latent.eff latent effect
#' @param res.var residual variance
#' @param init.seed initial seed
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
#################################################################################
mySimulations <- function(...) sim.hotspot(...)
#################################################################################


