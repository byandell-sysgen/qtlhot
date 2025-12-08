#' @export
#' @rdname CMSTtests
GetCisCandReg <- function(highobj, cand.reg, lod.thr = NULL)
{
  cand.names <- as.character(cand.reg[, 1])
  
  ## Restrict to being on same chromosome. This is fragile.
  cand.reg <- cand.reg[cand.reg[, 2] == cand.reg[, 4],]
  
  ## Restrict to LOD above lod.thr.
  highobj <- highlod.thr(highobj, lod.thr)
  
  chr.pos <- highobj$chr.pos
  
  ## Subset highlod to those phenos in cand.reg.
  pheno.cols <- unique(highobj$highlod$phenos)
  m.pheno <- match(as.character(cand.reg[,1]), highobj$names[pheno.cols])
  if(any(is.na(m.pheno)))
    stop("cannot find cand.reg traits in highobj$names")
  m.pheno <- pheno.cols[m.pheno]
  highlod <- highobj$highlod[highobj$highlod$phenos %in% m.pheno, ]
  
  ## Get start and end for each pheno. NB: may include multiple chr.
  h.index <- cumsum(table(highlod$phenos))
  h.index <- cbind(start = 1 + c(0, h.index[-length(h.index)]), end = h.index)
  ## Now get in right order.
  m.pheno <- order(m.pheno)
  h.index[m.pheno,] <- h.index
  
  ## Find lower and upper position around peak.
  tmpfn <- function(x, highlod, chr.pos) {
    h <- highlod[x[1]:x[2],, drop = FALSE]
    ## Only look at chr with peak LOD.
    wh <- which.max(h$lod)
    wh <- range(which(chr.pos$chr[h$row] == chr.pos$chr[h$row[wh]]))
    ## Could have non-contiguous regions. Don't sweat it for now.
    chr.pos$pos[h$row[wh]]
  }
  peak.pos <- t(apply(h.index, 1, tmpfn, highlod, chr.pos))
  dimnames(peak.pos)[[2]] <- c("peak.pos.lower", "peak.pos.upper")
  
  out <- data.frame(cand.reg, peak.pos)
  is.cis <- (out$phys.pos >= out$peak.pos.lower &
               out$phys.pos <= out$peak.pos.upper)
  ## Keep cis traits, but leave off peak.chr (since it == phys.chr).
  out <- out[is.cis, -4, drop = FALSE]
  if(nrow(out))
    attr(out, "cis.index") <- match(out[, 1], cand.names)
  out
}
GetCis <- function(x, window = 10) {
  xx <- x[x[, 2] == x[, 4],]
  xx <- xx[abs(xx[, 3] - xx[, 5]) <= window, ]
  index <- match(xx[, 1], x[, 1])
  list(cis.reg = xx, cis.index = index) 
}

