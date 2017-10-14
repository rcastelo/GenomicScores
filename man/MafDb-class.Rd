\name{MafDb-class}
\docType{class}
\alias{MafDb-class}
\alias{MafDb}

\alias{mafByOverlaps,MafDb-method}
\alias{mafByOverlaps}
\alias{mafById,MafDb-method}
\alias{mafById}
\alias{populations,MafDb-method}
\alias{populations}
\alias{organism,MafDb-method}
\alias{provider,MafDb-method}
\alias{providerVersion,MafDb-method}
\alias{referenceGenome,MafDb-method}
\alias{seqinfo,MafDb-method}
\alias{seqlengths,MafDb-method}
\alias{seqnames,MafDb-method}
\alias{seqlevelsStyle,MafDb-method}
\alias{show,MafDb-method}
\alias{$,MafDb-method}
\alias{citation,MafDb-method}

\title{MafDb class}

\description{
Class for annotation packages storing minor allele frequency data.
}
\usage{
\S4method{mafByOverlaps}{MafDb}(x, ranges, pop="AF", type=c("snvs", "nonsnvs"), caching=TRUE)
\S4method{mafById}{MafDb}(x, ids, pop, caching)
\S4method{populations}{MafDb}(x)
}
\arguments{
  \item{x}{A \code{MafDb} object.}
  \item{ranges}{Either a \code{GRanges} object, a \code{GPos} object or a character string vector with the format \code{"CHR:START[-END]"}.}
  \item{ids}{A character string vector with variant identifiers annotated by the MAF data source, typically dbSNP 'rs' identifiers. Note that the mapping of these identifiers to genomic positions and MAF values might be a subset of the most up to date dbSNP 'rs' identifier assignment to variants. To access the latter, please use the \code{snpsById()} method from the \code{BSgenome} package with the desired \code{SNPlocs.*} package.}
  \item{pop}{Character string vector with the populations for which we want to retrieve MAF values.}
  \item{type}{Character string setting the type of variant to seek, which can be either 'snvs' (default) when we seek single nucleotide variants or 'nonsnvs', otherwise.}
  \item{caching}{logical; TRUE (default) indicates that the function stores into main memeory the MAF data as it gets loaded from disk, improving performance; FALSE forces this function to load MAF data from disk each time, decreasing performance and memory requirements.}
}
\details{
The \code{MafDb} class is derived from the \code{\link{GScores}} class and it
serves the purpose of providing support to store and access minor allele frequency (MAF) data
from R and Bioconductor. Two annotation packages using the \code{MafDb} class are:

\tabular{ll}{
  \code{MafDb.1Kgenomes.phase1.hs37d5} \tab MAF values from the 1000 Genomes Project Phase 1. \cr
  \code{MafDb.1Kgenomes.phase3.hs37d5} \tab MAF values from the 1000 Genomes Project Phase 3. \cr
}

This object class tries to reduce the disk space required to store MAF values
for millions of SNPs by coding their double-precision values, which range
between 0 and 1, into a single-byte \code{raw} object type. To achieve this, the
original MAF values are rounded to one significant digit for AF < 0.1 and two
significant digits for AF >= 0.1. When a variant has multiple
alternate alleles, only the largest MAF value is stored.
}

\author{R. Castelo}
\examples{

\dontrun{
  ## lookup allele frequencies for rs1129038, a SNP associated to blue and brown eye colors
  ## as reported by Eiberg et al. Blue eye color in humans may be caused by a perfectly associated
  ## founder mutation in a regulatory element located within the HERC2 gene inhibiting OCA2 expression.
  ## Human Genetics, 123(2):177-87, 2008 [http://www.ncbi.nlm.nih.gov/pubmed/18172690]

  if (require(MafDb.1Kgenomes.phase1.hs37d5)) {
    mafdb <- MafDb.1Kgenomes.phase1.hs37d5
    mafdb

    ## specialized interface
    populations(mafdb)

    rng <- GRanges("15", IRanges(28356859, 28356859))
    mafByOverlaps(mafdb, rng)
    mafByOverlaps(mafdb, "15:28356859-28356859")
    mafByOverlaps(mafdb, "15:28356859")
    mafById(mafdb, "rs1129038")
  }
}

}

\keyword{datasets}
