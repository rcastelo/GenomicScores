\name{GenomicScores-deprecated}
\alias{GenomicScores-deprecated}
\alias{scores,GScores,GenomicRanges-method}
\alias{scores}

\title{Deprecated methods in package \sQuote{GenomicScores}}

\description{
  These functions are provided for compatibility with older versions
  of \sQuote{GenomicScores} only, and will be defunct at the next release.
}

\details{
  The following classes are deprecated and will be made defunct; use
  the replacement indicated below:

  \itemize{
    \item{\code{MafDb}:  The \code{MafDb} class has been replaced by the
          \code{\link{GScores-class}} class.}
  }

  The following methods are deprecated and will be made defunct; use
  the replacement indicated below:

  \itemize{
    \item{\code{scores()}:  The \code{scores()} method has been replaced by the
          \code{\link{gscores}()} and \code{\link{score}()} methods.}
    \item{\code{mafByOverlaps()}:  The \code{mafByOverlaps()} method has been
          replaced by the \code{\link{gscores}()} and \code{\link{score}()}
          methods.}
    \item{\code{mafById()}:  The \code{mafById()} method has been
          replaced by the \code{\link{gscores}()} and \code{\link{score}()}
          methods.}
  }
}
\author{
  Robert Castelo <robert.castelo@upf.edu>
}
\keyword{internal}
\keyword{misc}
