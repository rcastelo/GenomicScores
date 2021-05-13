\name{rgscores}

\alias{rgscores,GScores,missing-method}
\alias{rgscores,missing,GScores-method}
\alias{rgscores,numeric,GScores-method}
\alias{rgscores,integer,GScores-method}
\alias{rgscores}

\title{Sampling genomic gscores}

\description{
Function for randomly sampling genomic gscores from \code{GScores} objects.
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
  \item{...}{In the call to the \code{rgscores}() method one can additionally
                 set the following arguments:
                 \itemize{
                   \item{pop}{Character string vector of length one, specifying
                              the scores population from which we want to sample
                              scores. By default, its value is 'defaultPopulation(object)'.
                              Use \code{\link{populations}()} to find out the
                              available scores populations.}
                   \item{scores.only}{Flag setting whether only the scores should
                                      be returned as a numeric vector (\code{TRUE}),
                                      instead of returning them as a metadata column
                                      in a \code{GRanges} object (\code{FALSE}, default).}
                   \item{ranges}{Either a \code{GRanges} object or a character string
                                 vector of sequence names. Scores will be sampled from
                                 the given genomic regions, which by default correspond
                                 to the entire genomic space of \code{object}. Currently,
                                 only entire chromosomes are considered.}
                 }
            }
}
\details{
The method \code{rgscores()} samples scores randomly from a \code{GScores} object.
}
\value{
A \code{GRanges} object with the sampled genomic positions and scores. When
\code{scores.only=TRUE} then a numeric vector is returned with the sampled scores.
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

  phast <- phastCons100way.UCSC.hg19
  set.seed(123)
  rgscores(10L, phast, ranges=c("chr22", "chrY"))
}

if (require(MafDb.1Kgenomes.phase1.hs37d5)) {
  mafdb <- MafDb.1Kgenomes.phase1.hs37d5
  set.seed(123)
  rgscores(10L, mafdb, ranges=c("21", "22"))
}
}
\keyword{utilities}
