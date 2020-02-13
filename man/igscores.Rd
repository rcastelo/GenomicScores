\name{igscores}
\alias{igscores}

\encoding{latin1}

\title{
GenomicScores shiny app
}
\description{
Starts an interactive GenomicScores shiny web app.
}
\usage{
igscores()
}
\details{
The goal of the \code{GenomicScores} package is to provide support to store
and retrieve genomic scores associated to physical nucleotide positions along
a genome.

The \code{igscores()} function starts an interactive shiny web app that allows the
user to query annotation packages storing genomic scores. Internally, it calls to
the function \code{\link{gscores}()}; see its manual page for a description of the
arguments and their default and alternative values.
}
\value{
  None.
}
\references{
  Puigdevall, P. and Castelo, R. GenomicScores: seamless access to genomewide 
  position-specific scores from R and Bioconductor. 
  \emph{Bioinformatics}, 18:3208-3210, 2018.
}
\author{P. \enc{Rodríguez}{Rodriguez} and R. Castelo}
\seealso{
  \code{\link{gscores}}
}
\examples{

\dontrun{
igscores() ## this will open your browser with the GenomicScores shiny web app
}

}
\keyword{GenomicScores}
\keyword{shiny}
