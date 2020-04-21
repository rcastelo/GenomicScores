\name{rgscores}

\alias{rgscores,GScores,missing-method}
\alias{rgscores,missing,GScores-method}
\alias{rgscores,numeric,GScores-method}
\alias{rgscores,integer,GScores-method}
\alias{rgscores}

\title{Sampling genomic gscores}

\description{
Function to sample genomic gscores uniformly at random from \code{GScores} objects.
}
\usage{
\S4method{rgscores}{GScores,missing}(n, object, ...)
\S4method{rgscores}{missing,GScores}(n, object, ...)
\S4method{rgscores}{numeric,GScores}(n, object, ...)
\S4method{rgscores}{integer,GScores}(n, object, ...)
}
\arguments{
  \item{n}{Number of scores to sample.}
  \item{object}{A \code{GScores} object.}
  \item{...}{Further arguments.}
}
\details{
The method \code{rgscores()} samples scores uniformly at random
from a \code{GScores} object.
}
\value{
A \code{GRanges} object with the sampled genomic positions and scores.
}
\references{
Puigdevall, P. and Castelo, R. GenomicScores: seamless access to genomewide
position-specific scores from R and Bioconductor.
\emph{Bioinformatics}, 18:3208-3210, 2018.
}
\author{R. Castelo}
\seealso{
  \code{\link[phastCons100way.UCSC.hg19]{phastCons100way.UCSC.hg19}}
  \code{\link[MafDb.1Kgenomes.phase1.hs37d5]{MafDb.1Kgenomes.phase1.hs37d5}}
}
\examples{
## accessing genomic gscores from an annotation package
if (require(phastCons100way.UCSC.hg19)) {
  library(GenomicRanges)

  gsco <- phastCons100way.UCSC.hg19
  gsco
}

if (require(MafDb.1Kgenomes.phase1.hs37d5)) {
  mafdb <- MafDb.1Kgenomes.phase1.hs37d5
  mafdb
}
}
\keyword{utilities}
