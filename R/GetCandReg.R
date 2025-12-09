#' Get genetic information on candidate regulators and co-mapping traits.
#' 
#' Get chromosome (phys.chr) and physical position in cM (phys.pos), along with
#' the LOD score (peak.lod) at the peak position (peak.pos), and the chromosome
#' where the peak is located (peak.chr). Some candidates may map to the same
#' chromosome where they are physically located.
#' 
#' Traits that map to positions close to their physical locations are said to
#' map in cis (local linkages).  Traits that map to positions away from their
#' physical locations are said to map in trans (distal linkages). There is no
#' unambiguous way to determine how close a trait needs to map to its physical
#' location in order to be classified as cis. Our choice is to classify a trait
#' as cis if the 1.5-LOD support interval (Manichaikul et al. 2006) around the
#' LOD peak contains the trait's physical location, and if the LOD score at its
#' physical location is higher the the LOD threshold. The function
#' `GetCisCandReg` determines which of the candidate regulators map in
#' cis. The function `GetCoMappingTraits` returns a list with the putative
#' targets of each trait. A trait is included in the putative target list of a
#' trait when its LOD peak is greater than `lod.thr` and the `drop`
#' LOD support interval around the peak contains the location of the trait's
#' QTL. The function `JoinTestOutputs` currently relies on external files
#' that contain results of `FitAllTests`. It needs to be rewritten
#' to save space.
#' 
#' @param highobj data frame from `highlod`, which is sparse
#' summary of high LODs in large \code{\link[qtl]{scanone}} object
#' @param annot data frame with annotation information; must have first column
#' as unique identifier, third column as chromosome, and fifth column as
#' position in cM; typically column 2 has gene name, and column 4 has position
#' in Mb
#' @param traits names of traits to examine as candidate regulators; names must
#' correspond to phenotypes in `cross` object
#' @param cand.reg data frame with candidate regulator; see value section below
#' @param lod.thr LOD threshold for identifying candidate regions.
#' @return `GetCoMappingTraits` returns a list with each element being the
#' names of co-mapping traits for a particular name in `traits`.
#' `GetCandReg` returns a data frame while `GetCisCandReg` returns a
#' list with a similar candidate regulator data frame as the element
#' `cis.reg`, and the index of trait names as the element
#' `cis.index`. The elements of the candidate regulator data frame are as
#' follows (`peak.pos.lower` and `peak.pos.upper` only for
#' `GetCisCandReg`): \item{gene}{name of trait, which might be a gene
#' name} \item{phys.chr}{chromosome on which gene physically resides}
#' \item{phys.pos}{physical position (in cM)} \item{peak.chr}{chromosome where
#' peak LOD is located} \item{peak.pos}{position of peak (in cM)}
#' \item{peak.lod}{LOD value at peak}
#' \item{peak.pos.lower,peak.pos.upper}{lower and upper bounds of the 1.5-LOD
#' support interval around `peak.pos`}
#' @author Elias Chaibub Neto
#' @seealso `highlod`, `FitAllTests`,
#' \code{\link[qtl]{scanone}}
#' @references Manichaikul et al. (2006) Genetics
#' @keywords utilities
#' @examples
#' 
#' ## data(CMSTCross) is loaded lazily.
#' CMSTscan <- scanone(CMSTCross, pheno.col = 1:3, method = "hk")
#' CMSThigh <- highlod(CMSTscan)
#' traits <- names(CMSTCross$pheno)
#' annot <- data.frame(name = traits, traits = traits, chr = rep(1, 3),
#'  Mb.pos = c(55,10,100))
#' annot$cM.pos <- annot$Mb.pos
#' cand.reg <- GetCandReg(CMSThigh, annot, traits)
#' cis.cand.reg <- GetCisCandReg(CMSThigh, cand.reg)
#' comap.targets <- GetCoMappingTraits(CMSThigh, cand.reg)
#' 
#' @export
GetCandReg <- function(highobj, annot, traits)
{
  ## currently this only gets max over genome; want max by chr, yes?
  traits <- unique(traits)
  n <- length(traits)
  out <- data.frame(matrix(NA, n, 6))
  names(out) <- c("gene", "phys.chr", "phys.pos", "peak.chr", "peak.pos",
                  "peak.lod")
  
  ## Get annotation of trait name and physical chromosome and position (in cM).
  m <- match(traits, annot[,1])
  out[, 1:3] <- annot[m, c(1,3,5)]
  
  ## Get LOD peak information
  m <- !is.na(out[,3])
  
  pheno.cols <- match(traits[m], highobj$names)
  if(any(is.na(pheno.cols)))
    stop("some traits do not have scans")
  
  chr.pos <- highobj$chr.pos
  highlod <- highobj$highlod[highobj$highlod$phenos %in% pheno.cols,]
  tmp <- cumsum(table(highlod[,"phenos"]))
  tmp <- c(0, tmp[-length(tmp)])
  peak.index <- tmp + tapply(highlod$lod, highlod$phenos, which.max)
  
  ## now relate to lod, chr, pos, but get order right with traits
  m <- match(highobj$names[highlod[peak.index, "phenos"]], as.character(out[,1]))
  if(any(is.na(m)))
    stop("cannot match highlod with pheno names")
  
  out[m, 6] <- highlod[peak.index, "lod"]
  out[m, 4] <- chr.pos[highlod[peak.index, "row"], "chr"]
  out[m, 5] <- chr.pos[highlod[peak.index, "row"], "pos"]
  
  out[!is.na(out[,4]),, drop = FALSE]
}
##############################################################################
#' @export
#' @rdname GetCandReg
GetCoMappingTraits <- function(highobj, cand.reg)
{
  chr.pos <- highobj$chr.pos
  chrs <- unique(chr.pos$chr)
  chr <- ordered(cand.reg$peak.chr, chrs)
  phys <- cbind(phenos = match(as.character(cand.reg[,1]), highobj$names),
                chr = unclass(chr), pos = cand.reg$peak.pos)
  
  in.range <- function(pos, x.pos) {
    r <- range(pos)
    r[1] <= x.pos & r[2] >= x.pos
  }
  ## Find traits that are co-mapping.
  find.comap <- function(x, highlod, chr.pos, chrs, all.traits) {
    ## Traits must map to same chromosome.
    h <- highlod[chrs[x[2]] == chr.pos$chr[highlod$row],, drop = FALSE]
    ## And peaks must be in range.
    h <- tapply(chr.pos$pos[h$row], h$phenos, in.range, x[3])
    ## But have to remove trait from its list.
    h <- as.numeric(names(h[h]))
    h <- h[-match(x[1], h)]
    all.traits[h]
  }
  
  ## This list is too restrictive compared with earlier list of Elias.
  ## Try re-running deprecated code using scan.orf to compare.
  
  out <- apply(phys, 1, find.comap, highobj$highlod, chr.pos, chrs, highobj$names)
  if(is.matrix(out))
    out <- as.data.frame(out)
  names(out) <- as.character(cand.reg[,1])
  
  lapply(out, as.character)
}
