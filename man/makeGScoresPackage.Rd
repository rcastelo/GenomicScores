\name{makeGScoresPackage}

\alias{makeGScoresPackage}

\title{Building genomic gscores packages}

\description{
Build a genomic gscores packages from a \code{GScores} object.
}
\usage{
makeGScoresPackage(gsco, version, maintainer, author,
                   destDir=".", license="Artistic-2.0")
}
\arguments{
  \item{gsco}{\code{GScores}. The \code{GScores} object from which the package will be created.}
  \item{version}{\code{Character}. Version of the package.}
  \item{maintainer}{\code{Character}. Maintainer of the package, including email.}
  \item{author}{\code{Character}. Author of the package.}
  \item{destDir}{\code{Character}. The path to the directory where the package will be created.}
  \item{license}{\code{Character}. The license of the package.}
}
\details{
  This function allows one to create an R package from a \code{GScores} object. This may be useful if one wants to have a tar-ball package version of genomic scores available only through the AnnotationHub; see the vignette.
}
\value{
  It returns invisibily the package directory.
}
\references{
Puigdevall, P. and Castelo, R. GenomicScores: seamless access to genomewide
position-specific scores from R and Bioconductor.
\emph{Bioinformatics}, 18:3208-3210, 2018.
}
\author{R. Castelo}
\seealso{
  \code{\link{availableGScores}()}
  \code{\link{getGScores}()}
}
\examples{

## accessing genomic scores from AnnotationHub resources
## and building a package from them
\dontrun{
availableGScores()
gsco <- getGScores("fitCons.UCSC.hg19")
makeGScoresPackage(gsco, version="1.0", maintainer="me <me@example.com>", author="me")
}
}
\keyword{utilities}
