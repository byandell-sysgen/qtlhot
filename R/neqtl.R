#' neqtl.R Ported from http://github.com/kbroman/neqtl
#'
#' @param sigpos.out position(s) of max
#' @param chr chromosome as character
#' @param pos position on chromosome
#' @param win window width in cM
#'
#' @export
#' @importFrom broman runningmean
neqtl <- function(sigpos.out,chr,pos,win=5) {
  smoothall(sigpos.out,chr,pos,window=win)
}

smoothall <- function(themax, thechr, thepos, window=5)
{
  thesmooth <- vector("list", length(themax))
  names(thesmooth) <- names(themax)
  for(i in names(themax))
      thesmooth[[i]] <- smoothchr(themax[[i]], thepos[thechr==i], window=window)
   out <- NULL
  for(i in 1:length(thesmooth))
    out <- rbind(out, data.frame(chr=rep(names(themax)[i], nrow(thesmooth[[i]])),
                    pos=thesmooth[[i]][,1], nqtl=thesmooth[[i]][,2]))
  class(out) <- c("scanone", "data.frame")

  ## This chokes right here!

  rownames(out) <- paste("c", out[,1], ".loc", 1:nrow(out), sep="")

  out
}

## Uses positions from thepos for smoothing: ATB 9/10/09 ##
## In theloc by=0.2 was outside the seq() function--moved it inside  ATB 12/15/09 ##
smoothchr <- function(themax, thepos, window=5)
{
  ## theloc <- sort(unique(c(thepos, seq(0, max(thepos), by=0.2))))
  theloc <- thepos

  temploc <- c(themax, theloc)
  tempval <- c(rep(1, length(themax)), rep(0, length(theloc)))
  o <- order(temploc)
  temploc <- temploc[o]
  tempval <- tempval[o]
  smoothed <- broman::runningmean(temploc, tempval, at=thepos, window=window, what="sum")
  ## NB: This differs from R/neqtl, where at=theloc.
  cbind(thepos, smoothed)
}
