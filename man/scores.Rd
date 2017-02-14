\name{scores}

\alias{scores,GScores,GRanges-method}
\alias{scores}
\alias{availableGScores}
\alias{getGScores}

\title{Accessing genomic scores}

\description{
Functions to access genomic scores through \code{GScores} objects.
}
\usage{
availableGSscores()
getGScores(x)
\S4method{scores}{GScores,GRanges}(object, gpos, ...)
}
\arguments{
  \item{x}{A \code{character} vector of length 1 specifiying the genomic scores resource to
           fetch. The function \code{availableGScores()} shows the available genomic scores
           resources.}
  \item{object}{A \code{GScores} object.}
  \item{gpos}{A \code{GRanges} object with positions from where to retrieve
              genomic scores.}
  \item{...}{In the call to the \code{scores}() method one can additionally
             set the following arguments:
             \itemize{
               \item{summaryFun}{Function to summarize genomic scores when more than one
                                 position is retrieved. By default, this is set to the
                                 arithmetic mean, i.e., the \code{mean()} function.}
               \item{coercionFun}{Function to coerce the stored genomic scores, before the
                                  summary function is applied. By default genomic scores are
                                  coerced to real (\code{numeric}-class) values, i.e., to the
                                  \code{as.numeric()} function.}
               \item{caching}{Flag setting whether genomic scores per chromosome should be
                              kept cached in memory (\code{TRUE}, default) or not
                              (\code{FALSE}). The latter option minimizes the memory footprint
                              but slows down the performance when the \code{scores()} method
                              is called multiple times.}
             }}
}
\details{
The method \code{scores()} takes as first argument a \code{GScores-class} object
that can be loaded from an annotation package or from an \code{AnnotationHub} resource.
These two possibilities are illustrated in the examples below.
}

\author{R. Castelo}
\seealso{
  \code{\link[phastCons100way.UCSC.hg19]{phastCons100way.UCSC.hg19}}
  \code{\link[phastCons100way.UCSC.hg38]{phastCons100way.UCSC.hg38}}
}
\examples{
## accessing genomic scores from an annotation package
if (require(phastCons100way.UCSC.hg19)) {
  library(GenomicRanges)

  gsco <- phastCons100way.UCSC.hg19
  gsco
  scores(gsco, GRanges(seqnames="chr7", IRanges(start=117232380, width=5)))
}

## accessing genomic scores from AnnotationHub resources
\dontrun{
availableGScores()
gsco <- getGScores("phastCons100way.UCSC.hg19")
scores(gsco, GRanges(seqnames="chr7", IRanges(start=117232380, width=5)))
}
}
\keyword{datasets}
