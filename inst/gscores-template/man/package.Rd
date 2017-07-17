\name{@PKGNAME@}
\docType{package}

\alias{@PKGNAME@-package}
\alias{@PKGNAME@}

\title{@PKGTITLE@}

\description{
  This annotation package stores genomic scores for @ORGANISM@,
  genome version @GENOMEVERSION@.
}

\note{
  This data package was made from resources at @PROVIDER@
  based on version @PROVIDERVERSION@.
}

\author{@AUTHOR@}

\seealso{
  \code{\link[GenomicScores]{GenomicScores}}
}

\examples{
library(@PKGNAME@)

ls("package:@PKGNAME@")
@PKGNAME@
gsco <- @PKGNAME@
name(gsco)
type(gsco)
provider(gsco)
providerVersion(gsco)
organism(gsco)
referenceGenome(gsco)
seqlevelsStyle(gsco)
seqinfo(gsco)
head(seqnames(gsco))
head(seqlengths(gsco))
citation(gsco)
citation("@PKGNAME@")
}

\keyword{package}
\keyword{data}
