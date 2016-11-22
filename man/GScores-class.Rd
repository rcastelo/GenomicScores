\name{GScores-class}
\docType{class}
\alias{GScores-class}
\alias{GScores}
\alias{GenomicScores}

\alias{organism,GScores-method}
\alias{provider,GScores-method}
\alias{providerVersion,GScores-method}
\alias{referenceGenome,GScores-method}
\alias{scores,GScores,GRanges-method}
\alias{scores}
\alias{seqinfo,GScores-method}
\alias{seqlengths,GScores-method}
\alias{seqnames,GScores-method}
\alias{seqlevelsStyle,GScores-method}
\alias{show,GScores-method}
\alias{makePhastConsToAHM}

\title{The GenomicScores package and the GScores class}

\description{
Provide support to store and retrieve genomic scores associated
to nucleotide positions along a genome.
}
\usage{
GScores(provider, provider_version, download_url,
        download_date, reference_genome,
        data_pkgname, data_dirpath)
\S4method{scores}{GScores,GRanges}(object, gpos, ...)
}
\arguments{
  \item{provider}{Data provider.}
  \item{provider_version}{Version of the data, as given by the data provider.}
  \item{download_url}{URL from where data were downloaded.}
  \item{download_date}{Date when data were downloaded.}
  \item{reference_genome}{GenomeDescription object containing information about the reference genome to where data belongs to.}
  \item{data_pkgname}{Name of the package storing the data.}
  \item{data_dirpath}{Path to the directory where the actual data is stored.}
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
The \code{GScores} class and associated methods provide support to annotation
packages and resources that store genomic scores defined as score values each of
them associated to one physical nucleotide position along a genome, such as UCSC
phastCons conservation scores. Two such annotation packages are:

\describe{
  \item{\code{phastCons100way.UCSC.hg19}}{
    Nucleotide-level phastCons conservation scores from the UCSC Genome Browser
    calculated from multiple genome alignments from the human genome version hg19
    to 99 vertebrate species.
  }
  \item{\code{phastCons100way.UCSC.hg38}}{
    Nucleotide-level phastCons conservation scores from the UCSC Genome Browser
    calculated from multiple genome alignments from the human genome version hg38
    to 99 vertebrate species.
  }
}

The function \code{GScores()} is a constructor for a \code{GScores-class} object
and it is for internal use only. The end-user needs not to call this function.
}

\author{R. Castelo}
\seealso{
  \code{\link[phastCons100way.UCSC.hg19]{phastCons100way.UCSC.hg19}}
  \code{\link[phastCons100way.UCSC.hg38]{phastCons100way.UCSC.hg38}}
}
\examples{
if (require(phastCons100way.UCSC.hg19)) {
  library(GenomicRanges)

  phastCons100way.UCSC.hg19
  scores(phastCons100way.UCSC.hg19,
         GRanges(seqnames="chr7", IRanges(start=117232380, width=5)))
}
}
\source{
  \url{http://genome.ucsc.edu}
}
\keyword{datasets}
