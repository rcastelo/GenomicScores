\name{availableGScores}

\alias{availableGScores}
\alias{getGScores}

\title{Exploring genomic gscore resources}

\description{
Functions to explore genomic gscore resources.
}
\usage{
availableGScores(installed=FALSE, ah.only=FALSE)
getGScores(x)
}
\arguments{
  \item{x}{A\code{character} vector of length 1
           specifiying the genomic scores resource to fetch.}
  \item{installed}{A logical value specifying whether we want
                   to see the installed packages only that contain
                   genomic score resources. By default (\code{installed=FALSE})
                   all available packages with genomic score resources are
                   shown.}
  \item{ah.only}{A logical value specifying whether we want to see
                 the packages the are available as \code{AnnotationHub}
                 resources. By default (\code{ah.only=FALSE}) all available
                 packages with genomic score resources are shown. If
                 \code{ah.only=TRUE} this function will attempt to check
                 whether those genomic score resources are ready to be
                 downloaded.}
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
    \item{Installed: }{Whether the resource is installed (\code{TRUE}) or
                       not (\code{FALSE}).}
    \item{AH: }{Whether the resource is available through the
                \code{AnnotationHub} (\code{TRUE}) or not (\code{FALSE}).}
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
