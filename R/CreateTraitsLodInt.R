CreateTraitsLodInt <- function(scan, annot, traits, lod.thr, drop = 1.5)
{
  traits <- unique(traits)
  n <- length(traits)
  out <- data.frame(matrix(NA, n, 7))
  names(out) <- c("gene", "chr", "phys.pos", "lower.pos", "peak.pos",
                  "upper.pos", "lod")
  nms <- names(scan)[-c(1, 2)]
  for (i in 1 : n) {
    ii <- match(traits[i], annot[, 1])
    out[i, 1:3] <- annot[ii, c(1, 3, 5)]
    trait.chr <- annot[ii, 3]
    trait.pos <- annot[ii, 5]
    if (!is.na(trait.pos)) {
      peak <- max(scan[scan[, 1] == trait.chr, traits[i]])
      if(peak >= lod.thr){
        trait.index <- match(traits[i], nms)
        sscan <- scan[, c(1, 2, trait.index + 2)]
        lod.interval <- qtl::lodint(sscan, chr = trait.chr, drop)
        lb <- lod.interval[1, 2]
        ub <- lod.interval[3, 2]
        out[i, 4] <- lb
        out[i, 5] <- lod.interval[2, 2]
        out[i, 6] <- ub
        out[i, 7] <- peak
      }
    }     
    cat(" ", i, "\n")
  }
  subset(out, !is.na(out[, 4]))
}
