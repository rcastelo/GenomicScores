test_gscores <- function() {
  if (require(phastCons100way.UCSC.hg38)) {
    gr1 <- GRanges(seqnames="chr7", IRanges(start=117232380, width=5))

    gsco <- phastCons100way.UCSC.hg38
    res <- gscores(gsco, gr1)
    res$default <- round(res$default, digits=2)

    theoreticalres <- gr1
    seqlevels(theoreticalres) <- seqlevels(gsco)
    seqinfo(theoreticalres) <- seqinfo(gsco)
    theoreticalres$default <- 0.48

    checkIdentical(res, theoreticalres)
  }
}
