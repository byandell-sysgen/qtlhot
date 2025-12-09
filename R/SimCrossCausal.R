#' Simulate Cross for Causal Tests
#' 
#' Creates cross with certain pattern of dependence across phenotypes.
#' 
#' @param n.ind number of individuals to simulate
#' @param len vector specifying the chromosome lengths (in cM)
#' @param n.mar vector specifying the number of markers per chromosome
#' @param beta causal effect (slope) of first phenotype on others
#' @param add.eff,add.eff.1,add.eff.h additive genetic effect
#' @param dom.eff,dom.eff.1,dom.eff.h dominance genetic effect
#' @param sig2.1 residual variance for first phenotype
#' @param sig2.2,sig2.h residual variance for all other phenotypes
#' @param eq.spacing if \code{TRUE}, markers will be equally spaced
#' @param cross.type type of cross (\code{bc} and \code{f2} for now)
#' @param normalize normalize values if \code{TRUE}
#' @references Chaibub Neto E, Broman AT, Keller MP, Attie AD, Zhang B, Zhu J,
#' Yandell BS, Causal model selection hypothesis tests in systems genetics.
#' Genetics (in review).
#' @keywords utilities
#' @examples
#' 
#' set.seed(987654321)
#' CMSTCross <- SimCrossCausal(n.ind = 100, 
#'   len = rep(100, 3), n.mar = 101,
#'   beta = rep(0.5, 2), add.eff = 1, dom.eff = 0, 
#'   sig2.1 = 0.4, sig2.2 = 0.1, eq.spacing = FALSE, 
#'   cross.type = "bc", normalize = TRUE)
#' CMSTCross <- calc.genoprob(CMSTCross, step = 1)
#' \dontrun{
#' save(CMSTCross, file = "CMSTCross.RData", compress = TRUE)
#' }
#' 
#' @export
#' @importFrom qtl calc.genoprob find.marker pull.geno sim.cross sim.map
SimCrossCausal <- function(n.ind, len, n.mar, beta, add.eff, dom.eff, 
                           sig2.1 = 1, sig2.2 = 1, eq.spacing = FALSE, 
                           cross.type = c("bc", "f2"), normalize = FALSE) {
  n.traits <- length(beta)
  beta <- matrix(rep(beta, each = n.ind), n.ind, n.traits)
  Map <- qtl::sim.map(len, n.mar, eq.spacing = eq.spacing, include.x = FALSE)
  Cross <- qtl::sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- qtl::pull.geno(Cross)
  q <- mygeno[, "D1M51"]
  
  cross.type <- match.arg(cross.type)
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- add.q * add.eff + stats::rnorm(n.ind, 0, sqrt(sig2.1))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- add.q * add.eff + dom.q * dom.eff + stats::rnorm(n.ind, 0, sqrt(sig2.1))
  }
  y <- beta * y1 + matrix(stats::rnorm(n.ind * n.traits, 0, sqrt(sig2.2)), n.ind, n.traits)
  y <- data.frame(y1, y)
  names(y) <- paste("y", 1 : (n.traits + 1), sep = "")
  if (normalize) {
    apply(y, 2, normal.trans)
  }
  Cross$pheno <- y
  Cross
}
################################################################################
SimCross1 <- function(n.ind, mu, beta21, add.eff1, dom.eff1, 
                      sig2.1 = 1, sig2.2 = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- qtl::sim.map(len = rep(100,3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- qtl::sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- qtl::pull.geno(Cross)
  q <- mygeno[, "D1M51"]
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + stats::rnorm(n.ind, 0, sqrt(sig2.1))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + 
          stats::rnorm(n.ind, 0, sqrt(sig2.1))
  }
  y2 <- mu + beta21 * y1 + stats::rnorm(n.ind, 0, sqrt(sig2.2))
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
################################################################################
SimCross2 <- function(n.ind, mu, beta21, beta1h, beta2h, add.eff1, dom.eff1, 
                      sig2.1 = 1, sig2.2 = 1, sig2.h = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- qtl::sim.map(len = rep(100,3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- qtl::sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- qtl::pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  h <- mu + stats::rnorm(n.ind, 0, sqrt(sig2.h))
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + beta1h * h + stats::rnorm(n.ind, 0, sqrt(sig2.1))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + beta1h * h + 
          stats::rnorm(n.ind, 0, sqrt(sig2.1))
  }
  y2 <- mu + beta21 * y1 + beta2h * h + stats::rnorm(n.ind, 0, sqrt(sig2.2))
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
##############################################################################
SimCross3 <- function(n.ind, mu, beta21, add.eff1, dom.eff1, add.eff2,
                      dom.eff2, sig2.1 = 1, sig2.2 = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- qtl::sim.map(len = rep(100, 3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- qtl::sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- qtl::pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + stats::rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + beta21 * y1 + stats::rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + 
          stats::rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + dom.q * dom.eff2 + beta21 * y1 + 
          stats::rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
##############################################################################
SimCross4 <- function(n.ind, mu, add.eff1, dom.eff1, add.eff2, dom.eff2, 
                      sig2.1 = 1, sig2.2 = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- qtl::sim.map(len = rep(100, 3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- qtl::sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- qtl::pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + stats::rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + stats::rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + 
          stats::rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + dom.q * dom.eff2 + 
          stats::rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
##############################################################################
SimCross5 <- function(n.ind, mu, add.eff1, dom.eff1, add.eff2, dom.eff2, 
                      beta1h, beta2h, sig2.1 = 1, sig2.2 = 1, sig2.h = 1, 
                      eq.spacing = FALSE, cross.type = "f2", 
                      normalize = FALSE) {
  Map <- qtl::sim.map(len = rep(100, 3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- qtl::sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- qtl::pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  h <- mu + stats::rnorm(n.ind, 0, sqrt(sig2.h))
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + h * beta1h + stats::rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + h * beta2h + stats::rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + h * beta1h + 
          stats::rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + dom.q * dom.eff2 + h * beta2h + 
          stats::rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
##############################################################################
SimCross6 <- function(n.ind, mu, add.eff, dom.eff, beta1h, beta2h, 
                      sig2.1 = 1, sig2.2 = 1, sig2.h = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- qtl::sim.map(len = rep(100, 3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- qtl::sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- qtl::pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  if (cross.type == "bc") {
    add.q <- q - 1.5
    h <- mu + add.q * add.eff + stats::rnorm(n.ind, 0, sqrt(sig2.h))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    h <- mu + add.q * add.eff + dom.q * dom.eff + 
         stats::rnorm(n.ind, 0, sqrt(sig2.h))
  }
  y1 <- mu + h * beta1h + stats::rnorm(n.ind, 0, sqrt(sig2.1))
  y2 <- mu + h * beta2h + stats::rnorm(n.ind, 0, sqrt(sig2.2))
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
