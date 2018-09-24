\name{GenomicScores-defunct}
\alias{GenomicScores-defunct}
\alias{scores,GScores,GenomicRanges-method}
\alias{scores}

\title{Defunct methods in package \sQuote{GenomicScores}}

\description{
  Defunct classes and methods in \sQuote{GenomicScores}.
}

\details{
  The following classes are defunct; use the replacement indicated below:

  \itemize{
    \item{\code{MafDb}:  The \code{MafDb} class has been replaced by the
          \code{\link{GScores-class}} class.}
  }

  The following methods are defunct; use the replacement indicated below:

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
