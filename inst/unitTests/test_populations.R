test_populations <- function() {
  if (require(phastCons100way.UCSC.hg19)) {
    gsco <- phastCons100way.UCSC.hg19
    pops <- populations(gsco)

    checkIdentical(pops, c("default", "DP2"))
  }
}
