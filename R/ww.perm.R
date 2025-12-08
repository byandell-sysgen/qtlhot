## 
#' Conduct West-Wu (Q) permutation tests
#' 
#' This function computes the West/Wu permutation thresholds.
#' The output is a nlod (number of LOD thresholds) by nalpha (number of 
#' significance levels) matrix, where each entry shows the hotspot size 
#' significance threshold of the West/Wu approach. Note we have two "alphas" 
#' here, one for the QTL mapping (the LOD thresholds) and one for the 
#' permutation significance (alpha levels of lod.thrs).
#' 
#' Note that I separated the original ww.permutations() into a piece that do the 
#' actual permutations [ww.perm.matrix() function] and a piece that summarizes it
#' [the ww.summary() function] in the same way you did with the NL.N.permutations()
#' function.  
#' 
#' Perform permutation tests to assess the statistical significance of the
#' hotspots detected using the West-Wu \code{Q}-method permutations. The
#' \code{ww.perm} function implements the \code{Q}-method's permutation scheme
#' (see the Method's section of Chaibub Neto et a. 2012, for details). The
#' \code{n.perm} parameter specifies the number of simulations. Here we set it
#' to 100 in order to save time. In practice, we recommend at least 1,000
#' permutations. The function's output is a matrix with 100 rows representing
#' the permutations, and 10 columns representing the QTL mapping thresholds.
#' Each entry \code{ij}, represents the maximum number of significant linkages
#' across the entire genome detected at permutation \code{i}, using the LOD
#' threshold \code{j}. The \code{ww.summary} function computes the Q-method's
#' hotspot size permutation thresholds, that is, the \code{1-alpha} quantiles
#' for each one of the QTL mapping LOD thrsholds in \code{lod.thrs}. For
#' instance, the entry at row 10 and column 1 of the \code{Q.1.thr} matrix
#' tells us that the 99\% percentile of the permutation distribution of genome
#' wide maximum hotspot size based on a QTL mapping threshold of 2.11 is 27.00.
#' In other words, any hotspot greater than 27 is considered statistically
#' significant at a 0.01 significance level when QTL mapping is done using a
#' 2.11 LOD threshold.  In general, we are often interested in using the same
#' error rates for the QTL mapping and hotspot analysis. That is, if we adopt a
#' QTL mapping threshold that controls GWER at a 1\% level (in our case, 3.11)
#' we will also want to consider \code{alpha = 0.01} for the hotspot analysis,
#' leading to a hotspot threshold of 12.00. Therefore, we are usually more
#' interested in the diagonal of \code{Q.1.thr}. We adopted a GWER of 5\%, and
#' the corresponding \code{Q}-method's permutation threshold is 18. According
#' to this threshold, all hotspots are significant.
#' 
#' @param highobj object of class \code{\link{highlod}}
#' @param n.perm number of permutations
#' @param lod.thrs vector of LOD thresholds
#' @param alpha.levels vector of significance levels
#' @param x,object object of class \code{ww.perm}
#' @param \dots ignored
#' @param verbose verbose output if \code{TRUE}
#' @author Elias Chaibub Neto and Brian S Yandell
#' @keywords utilities
#' @examples
#' 
#' \dontrun{
#' ## All unspecified objects come from vignette qtlhot.
#' set.seed(12345)
#' Q.1 <- ww.perm(high1, n.perm = 100, lod.thrs, alphas)
#' Q.1.thr <- summary(Q.1, alphas)
#' Q.1.thr
#' diag(Q.1.thr)
#' 
#' set.seed(12345)
#' Q.2 <- ww.perm(high2, 100, lod.thrs, alphas)
#' Q.2.thr <- summary(Q.2, alphas)
#' }
#' 
ww.perm <- function(highobj, n.perm, lod.thrs, alpha.levels, verbose = FALSE)
  ww.perm.highlod(highobj, n.perm, lod.thrs, alpha.levels, verbose)
####################################################################################
#' @method print ww.perm
#' @rdname ww.perm
#' @export
print.ww.perm <- function(x, ...) print(summary(x, ...))
#' @method summary ww.perm
#' @rdname ww.perm
#' @export
summary.ww.perm <- function(object, alpha.levels = attr(object, "alpha.levels"), ...)
{
  nalpha <- length(alpha.levels)
  ww.thrs <- t(apply(object, 2, quantile_highlod, 1 - alpha.levels))
  dimnames(ww.thrs) <- list(dimnames(object)[[2]], as.character(alpha.levels))
  ww.thrs
}
####################################################################################
ww.perm.highlod <- function(highobj, n.perm, lod.thrs, alpha.levels, verbose = FALSE)
{
  ## These permutations differ from Elias's original (see inst/deprecated.R):
  ## 1. Only sampling on subset of phenotypes with peaks.
  ##    Thus the number of random sample calls is smaller.
  ## 2. chr.pos -> sample(chr.pos) rather than lod <- sample(lod).
  ##    Thus shuffle in inverted: s <- sample(1:10) vs. o <- order(s)
  ## 3. Original code shuffled scanmat repeatedly
  ##    rather than new shuffle each time. This is a minor bug rather than a feature.
  ## 4. New code is much faster.

  ## Both versions break up the support intervals with shuffling.
  ## This can lead to more peaks across genome.
  ## ww.perm.maxlod (below) resolves this by focusing only on peaks,
  ## but only makes sense if ww.perm is used with window = 0.
  
  highobj <- highlod.thr(highobj, min(lod.thrs))
  
  phenos <- sort(unique(highobj$highlod$pheno))
  n.chrpos <- nrow(highobj$chr.pos)

  nlod <- length(lod.thrs)
  max.ww <- matrix(NA, n.perm, nlod)
  dimnames(max.ww) <- list(NULL, as.character(lod.thrs))
  
  ## Want to replace row in highobj$highlod with permuted row.
  ## With separate permutation by pheno of seq(n.chrpos).

  ## Shuffle rows of chr.pos. Done separately for each pheno using tapply below.
  mysam <- function(row, n) sample(n)[row]

  ## Order highlod by pheno to make it easier
  highobj$highlod <- highobj$highlod[order(highobj$highlod$pheno), ]
  highs <- highobj
  
  for(i in 1:n.perm) {
    if(verbose)
      cat("\n", i, "")

    ## Permute columns of the scan object separately by trait.
    mycat("sample...", verbose, last = "")
    row.perm <- unlist(tapply(highobj$highlod$row, highobj$highlod$pheno, mysam, n.chrpos))
    highs$highlod$row <- row.perm
    
    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("max...", verbose, last = "")
    tmp <- max(highs, lod.thr = lod.thrs)
    if(length(tmp)) {
      m <- match(tmp$lod.thr, lod.thrs)
      max.ww[i, m] <- tmp$max.N
    }
  }
  if(verbose) cat("\n")
  
  class(max.ww) <- c("ww.perm", "list")
  attr(max.ww, "lod.thrs") <- lod.thrs
  attr(max.ww, "alpha.levels") <- alpha.levels
  
  max.ww
}
####################################################################################
ww.perm.maxlod <- function(highobj, n.perm, lod.thrs, alpha.levels, verbose = FALSE)
{
  highobj <- highlod.thr(highobj, min(lod.thrs))
  maxobj <- pull.max(highobj)
  
  phenos <- sort(unique(highobj$highlod$pheno))
  n.chrpos <- nrow(highobj$chr.pos)

  nlod <- length(lod.thrs)
  max.ww <- matrix(NA, n.perm, nlod)
  dimnames(max.ww) <- list(NULL, as.character(lod.thrs))
 
  ## Shuffle rows of chr.pos. Done separately for each pheno using tapply below.
  mysam <- function(row, n) sample(n)[row]
 
  ## Order highlod by pheno to make it easier
  maxobj$highlod <- maxobj$highlod[order(maxobj$highlod$pheno), ]
  highs <- maxobj
  
  for(i in 1:n.perm) {
    if(verbose)
      cat("\n", i)

    ## Permute columns of the scan object separately by trait.
    mycat("sample...", verbose, last = "")
    row.perm <- unlist(tapply(maxobj$highlod$row, maxobj$highlod$pheno, mysam, n.chrpos))
    highs$highlod$row <- row.perm
    
    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("max...", verbose, last = "")
    tmp <- max(highs, lod.thr = lod.thrs)
    if(length(tmp)) {
      m <- match(tmp$lod.thr, lod.thrs)
      max.ww[i, m] <- tmp$max.N
    }
  }
  class(max.ww) <- c("ww.perm", class(max.ww))
  attr(max.ww, "lod.thrs") <- lod.thrs
  attr(max.ww, "alpha.levels") <- alpha.levels
  
  max.ww
}
######################################################################
pull.max <- function(highobj)
{
  highobj$highlod <- highobj$highlod[order(highobj$highlod$pheno), ]
  chr.pheno <- ordered(paste(highobj$chr.pos$chr[highobj$highlod$row],
                             highobj$highlod$pheno, sep = "."))
  tmpfn <- function(x) sample(which.max(x), 1)
  tmp <- unlist(tapply(highobj$highlod$lod, chr.pheno, which.max))
  tmp2 <- table(chr.pheno)
  highobj$highlod <- highobj$highlod[tmp + c(0, cumsum(tmp2)[-length(tmp2)]), ]
  highobj
}




