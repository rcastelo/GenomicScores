\name{availableGScores}

\alias{availableGScores}
\alias{getGScores}

\title{Exploring genomic scores resources}

\description{
Functions to explore genomic scores resources.
}
\usage{
availableGScores(use.internet=FALSE)
getGScores(x)
}
\arguments{
  \item{x}{A\code{character} vector of length 1
           specifiying the genomic scores resource to fetch.}
  \item{use.internet}{A logical value specifying whether we want
                      to check through the internet whether the
                      annotation packages and \code{AnnotationHub}
                      resources are ready to be downloaded. By default,
                      i.e. \code{use.internet=FALSE}, this is not
                      checked to speed up the execution of this function.}
}
\details{
The function \code{availableGScores()} shows genomic score sets available as
\code{AnnotationHub} online resources.
}
\value{
The function \code{availableGScores()} returns a \code{data.frame} object
with a row for each available resource of genomic scores and the following
four columns:
  \itemize{
    \item{Name: }{Name of the Bioconductor package or \code{AnnotationHub}
                  resource.}
    \item{Organism: }{Organism on which the genomic scores are defined.}
    \item{Category: }{Category to which the genomic scores belong to.}
    \item{Installed: }{Whether the resource is installed as a package (\code{TRUE})
                       or not (\code{FALSE}).}
    \item{Cached: }{Whether the resource is available within the local cache
                    of the \code{AnnotationHub} (\code{TRUE}) or not (\code{FALSE}).}
    \item{BiocManagerInstall: }{Whether the resource can be installed as an
                                annotation package through \code{BiocManager::install()}
                                (\code{TRUE}) or not (\code{FALSE}).}
    \item{AnnotationHub: }{Whether the resource can be downloaded as a \code{GScores}
                           object through the \code{AnnotationHub}, using the function
                           \code{\link{getGScores}()} (\code{TRUE}), or not (\code{FALSE}).
                           (\code{TRUE}) or not (\code{FALSE}).}
  }
The function \code{getGScores()} returns a \code{GScores} object.
}
\references{
Puigdevall, P. and Castelo, R. GenomicScores: seamless access to genomewide
position-specific scores from R and Bioconductor.
\emph{Bioinformatics}, 18:3208-3210, 2018.
}
\author{R. Castelo}
\seealso{
  \code{\link{getGScores}()}
  \code{\link[phastCons100way.UCSC.hg19]{phastCons100way.UCSC.hg19}}
  \code{\link[MafDb.1Kgenomes.phase1.hs37d5]{MafDb.1Kgenomes.phase1.hs37d5}}
}
\examples{
availableGScores()

\dontrun{
gsco <- getGScores("cadd.v1.3.hg19")
}
}
\keyword{utilities}
