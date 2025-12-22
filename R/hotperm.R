#' Conduct NL and N permutation tests
#' 
#' This set of functions compute the permutation LOD thresholds for the NL-method
#' and the permutation hotspot size thresholds for the N-method. The output
#' is a list with two elements: the NL- and the N-method's threshold matrices.
#' The NL-method output is a nN (number of spurious hotspot sizes) by nalpha
#' (number of significance levels) threshold matrix. Note that for the
#' NL-method we have a single "alpha" since we use the same significance level
#' for QTL mapping and permutation significance.
#' The N-method output is a nlod (number of LOD thresholds) by nalpha (number
#' of significance levels) threshold matrix. Note that here we have two
#' "alphas", one for the QTL mapping (the LOD thresholds) and one for the 
#' permutation significance (alpha levels). 
#' 
#' @param cross object of class \code{cross}
#' @param n.quant maximum of \code{s.quant}
#' @param n.perm number of permutations
#' @param lod.thrs vector of LOD thresholds
#' @param alpha.levels vector of significance levels
#' @param quant.levels quantile levels, as number of traits, to show in
#' summary; default is 1, 2, 5, 10, \dots up to maximum recorded
#' @param drop.lod LOD drop amount for support intervals
#' @param window window size for smoothed hotspot size
#' @param verbose verbose output if \code{TRUE}
#' @param init.seed initial seed for pseudo-random number generation
#' @param x,object object of class \code{hotperm} or \code{summary.hotperm}
#' @param probs probability levels for quantiles (\code{1-probs} if all > 0.5);
#' default is \code{alpha.levels}
#' @param addcovar additive covariates as vector or matrix; see
#' \code{\link[qtl]{scanone}}
#' @param intcovar interactive covariates as vector or matrix; see
#' \code{\link[qtl]{scanone}}
#' @param level Significance level for hotspot detection.
#' @param \dots arguments passed along to \code{scanone}
#' @author Elias Chaibub Neto and Brian S Yandell
#' @keywords utilities
#' @examples
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
#' pt.scanone <- scanone(ncross1, method = "hk", n.perm = 1000)
#' alphas <- seq(0.01, 0.10, by=0.01)
#' lod.thrs <- summary(pt.scanone, alphas)
#' # This takes awhile, so we save the object.
#' \dontrun{
#' hotperm1 <- hotperm(cross = cross1,
#'                     n.quant = 300,
#'                     n.perm = 100,
#'                     lod.thrs = lod.thrs,
#'                     alpha.levels = alphas,
#'                     drop.lod = 1.5,
#'                     verbose = FALSE)
#' save(hotperm1, file = "hotperm1.RData", compress = TRUE)
#' # data(hotperm1) 
#' summary(hotperm1)
#' }
#' @export
#' @importFrom qtl nind nphe scanone 
hotperm <- function(cross, n.quant, n.perm, lod.thrs, alpha.levels, drop.lod = 1.5,
                    window = NULL, verbose = FALSE, init.seed = 0,
                    addcovar = NULL, intcovar = NULL, ...) 
{
  set.seed(init.seed)
  n.phe <- qtl::nphe(cross)
  pheno.col <- seq(n.phe)
  n.ind <- qtl::nind(cross)
  n.quant <- min(n.quant, n.phe)

  tmp <- table(sapply(cross$pheno, class))
  if(length(tmp) > 1 || names(tmp)[1] != "numeric")
    stop("all phenotypes in cross object must be numeric")
  
  n.lod <- length(lod.thrs)

  max.N <- matrix(0, n.perm, n.lod)
  dimnames(max.N) <- list(NULL, as.character(lod.thrs))
  if(is.null(window))
    max.N.window <- NULL
  else
    max.N.window <- max.N

  max.lod.quant <- matrix(0, n.perm, n.quant)
  dimnames(max.lod.quant) <- list(NULL, as.character(seq_len(n.quant)))

  for(i in 1:n.perm){
    if(verbose)
      cat("\n", i, "")
    
    ## permute rows of the phenotype data matrix
    perm.cross <- cross
    tmp <- sample(c(1:n.ind), n.ind, replace=FALSE)
    perm.cross$pheno <- cross$pheno[tmp,]

    ## perform mapping analysis in the permuted data
    mycat("nind...", verbose, last = "")
    ## NB: qtl::scanone groups phenos in batches based on missing data patterns.
    scanmat <- qtl::scanone(perm.cross, pheno.col = pheno.col, method = "hk", 
                        addcovar = addcovar, intcovar = intcovar, ...)

    ## Reduce to high LOD scores.
    mycat("highlod...", verbose, last = "")
    highs <- highlod(scanmat, min(lod.thrs), drop.lod, restrict.lod = TRUE)
    rm(scanmat)
    gc()

    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("max...", verbose, last = "")
    maxhi <- max(highs, lod.thr = lod.thrs, window = window)
    max.N[i, ] <- maxhi$max.N
    if(!is.null(window))
      max.N.window[i, ] <- maxhi$max.N.window
    
    ## get the maximum lod-quantile across the genome
    ## rows indexes the permutations
    ## columns indexes the s.quant quantiles
    mycat("quantile...", verbose, last = "")
    tmp <- quantile_highlod(highs, n.quant = n.quant)
    if(length(tmp))
      max.lod.quant[i, seq(tmp)] <- tmp
  }
  if(verbose) cat("\n")
  
  out <- list(max.N = max.N, max.N.window = max.N.window,
              max.lod.quant = max.lod.quant)
  class(out) <- c("hotperm", "list")
  attr(out, "lod.thrs") <- lod.thrs
  attr(out, "alpha.levels") <- alpha.levels
  
  out
}
#' @method print hotperm
#' @rdname hotperm
#' @export
print.hotperm <- function(x, ...) print(summary(x, ...))
#' @method summary hotperm
#' @rdname hotperm
#' @export
summary.hotperm <- function(object, quant.levels, ...)
{
  out <- quantile_hotperm(object, ...)

  attr(out, "window") <- attr(object, "window")

  alpha.levels <- attr(object, "alpha.levels")
  if(max(alpha.levels) < 0.5)
    alpha.levels <- 1 - alpha.levels
  n.quant <- ncol(object$max.lod.quant)
  if(missing(quant.levels)) {
    quant.levels <- log10(n.quant)
    quant.levels <- round(10 ^ c(outer(log10(c(1,2,5)), seq(0, floor(quant.levels)), "+")))
  }
  quant.levels <- quant.levels[quant.levels <= n.quant]
  if(max(quant.levels) < n.quant)
    quant.levels <- c(quant.levels, n.quant)
  out$max.lod.quant <- 
    t(apply(object$max.lod.quant[, quant.levels, drop = FALSE],
            2, stats::quantile, probs = alpha.levels, na.rm = TRUE))
  
  class(out) <- c("summary.hotperm", class(out))
  out
}
#' @method print summary.hotperm
#' @rdname hotperm
#' @export
print.summary.hotperm <- function(x, ...)
{
  cat("max.N: hotspot threshold by single-trait LOD threshold and significance level\n")
  print(ceiling(x$max.N))
  window <- attr(x, "window")
  if(!is.null(window)) {
    cat(paste("\nmax.N.window: smoothed hotspot threshold by single-trait LOD threshold and significance level ",
              "(window = ", window, ")\n", sep = ""))
    print(ceiling(x$max.N.window))
  }
  cat("\nmax.lod.quant: LOD threshold by hotspot size quantile and significance level\n")
  print(round(x$max.lod.quant, 2))
  invisible()
}
#' @method plot hotperm
#' @rdname hotperm
#' @export
plot.hotperm <- function(x, probs = seq(0.9, 0.99, by = 0.01), level = 0.95, ...)
{
  lod.thrs <- attr(x, "lod.thrs")
  alpha.levels <- attr(x, "alpha.levels")
  if(max(alpha.levels) <= 0.5)
    alpha.levels <- 1 - alpha.levels
  if(max(probs) <= 0.5)
    probs <- 1 - probs
  if(level < 0.5)
    level <- 1 - level
  
  lod.thr <- lod.thrs[which.min(abs(level - alpha.levels))[1]]
  out <- quantile_hotperm(x, probs, lod.thr = lod.thr, ...)

  wh.thr <- which.min(abs(lod.thr - lod.thrs))[1]
  wh.level <- which.min(abs(level - probs))[1]
  
  tmp.plot <- function(x.vals, quant, x.crit, probs, level, wh.thr,
                       is.quantile = FALSE, main = "", add.level = FALSE)
  {
    n.probs <- length(probs)
    wh.level <- which.min(abs(level - probs))[1]
    quant.thr <- quant[wh.thr, wh.level]

    xlabs <- "single trait LOD threshold"
    if(is.quantile)
      xlabs <- paste(xlabs, "quantile")
    
    graphics::plot(range(x.vals), c(0, max(quant)), type = "n", xlab = "", ylab = "")
    graphics::mtext(xlabs, 1, 2)
    graphics::mtext("hotspot size", 2, 2)
    graphics::abline(v = x.crit, col = "darkgray", lty = 2)
    graphics::abline(h = quant.thr, col = "darkgray", lty = 2)
    graphics::mtext(ceiling(quant.thr), 2, at = quant.thr, las = 2, cex = 1)
    for(i in seq(n.probs)) {
      graphics::lines(rev(sort(x.vals)), quant[,i],
            lwd = 1 + 2 * (round(probs[i] - level, 2) == 0))
    }

    graphics::text(x.crit, quant[wh.thr, n.probs] + 5, 1 - max(probs), adj = 0)
    graphics::text(x.crit, quant[wh.thr, 1] - 5, 1 - min(probs), adj = 1)
    graphics::text(graphics::par("usr")[1], quant.thr + 5, 1 - level, adj = 0)

    if(add.level)
      main <- paste(main, "\n hotspot size significance level =", 1 - max(probs), "to", 1 - min(probs))
    graphics::mtext(main, 3, 0.5)
  }

  tmpar <- graphics::par(mfrow = c(1 + !is.null(out$max.N.window),2), mar = c(3.1,3.1,2.1,0.1))
  if(!is.null(out$max.N.window)) {
    ## Jansen method, smoothing.
    tmp.plot(lod.thrs, out$max.N.window, lod.thr, probs, level, wh.thr, FALSE,
             "Jansen method 5cM window")
    
    tmp.plot(probs, out$max.N.window, level, probs, level, wh.thr, TRUE,
             "Jansen method 5cM window")
  }

  tmp.plot(lod.thrs, out$max.N, lod.thr, probs, level, wh.thr, FALSE,
           "Jansen method per locus")
  tmp.plot(probs, out$max.N, level, probs, level, wh.thr, TRUE,
           "Jansen method per locus")
  graphics::par(tmpar)

  if(!is.null(out$max.lod.quant)) {
    n.quant <- nrow(out$max.lod.quant)
    n.probs <- length(probs)    
    quant.thr <- max(which(out$max.lod.quant[,wh.level] >= lod.thr), na.rm = TRUE)
    
    ## Chaibub Neto method.
    graphics::plot(c(1,n.quant), range(out$max.lod.quant, na.rm = TRUE), type = "n",
         xlab = "significant hotspot size with given threshold",
         ylab = "hotspot LOD score threshold",
         log = "xy")
    graphics::abline(h = lod.thr, col = "darkgray", lty = 2)
    graphics::abline(v = quant.thr, col = "darkgray", lty = 2)
    graphics::mtext(ceiling(quant.thr), 1, at = quant.thr, las = 2)
    for(i in seq(n.probs)) {
      if(any(!is.na(out$max.lod.quant[,i])))
        graphics::lines(seq(n.quant), out$max.lod.quant[,i],
              lwd = 1 + 2 * (round(probs[i] - level, 2) == 0))
    }
    
    n.thr2 <- length(lod.thrs) / 2
    graphics::text(n.thr2 + 1, out$max.lod.quant[n.thr2, n.probs], 1 - max(probs), adj = 0)
    graphics::text(n.thr2 - 1, out$max.lod.quant[n.thr2, 1], 1 - min(probs), adj = 1)
    graphics::text(graphics::par("usr")[1], lod.thr, level, adj = 0)
    graphics::title(paste("hotspot LOD threshold by hotspot size\nsignificance level =",
                1 - max(probs), "to", 1 - min(probs)))
  }
}
