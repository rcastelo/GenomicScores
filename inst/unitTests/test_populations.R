test_populations <- function() {
  if (require(phastCons100way.UCSC.hg38)) {
    gsco <- phastCons100way.UCSC.hg38

    pops <- populations(gsco)

    checkIdentical(pops, c("default", "DP2"))

    checkIdentical(unique(genome(gsco)), "Genome Reference Consortium GRCh38")
    checkIdentical(provider(gsco), "UCSC")
  }
}
