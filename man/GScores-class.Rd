\name{GScores-class}
\docType{class}
\alias{GScores-class}
\alias{GScores}
\alias{GenomicScores}
\alias{name,GScores-method}
\alias{name}
\alias{type,GScores-method}
\alias{type}
\alias{organism,GScores-method}
\alias{provider,GScores-method}
\alias{providerVersion,GScores-method}
\alias{referenceGenome,GScores-method}
\alias{seqinfo,GScores-method}
\alias{seqlengths,GScores-method}
\alias{seqnames,GScores-method}
\alias{seqlevelsStyle,GScores-method}
%%\alias{qfun,GScores-method}
%%\alias{qfun}
%%\alias{dqfun,GScores-method}
%%\alias{dqfun}
\alias{show,GScores-method}

\title{The GenomicScores package and the GScores class}

\description{
Provide support to store and retrieve genomic scores associated
to nucleotide positions along a genome.
}
\usage{
GScores(provider, provider_version, download_url,
        download_date, reference_genome,
        data_pkgname, data_dirpath, data_serialized_objnames)
}
\arguments{
  \item{provider}{Data provider.}
  \item{provider_version}{Version of the data, as given by the data provider.}
  \item{download_url}{URL from where data were downloaded.}
  \item{download_date}{Date when data were downloaded.}
  \item{reference_genome}{GenomeDescription object containing information about the reference genome to where data belongs to.}
  \item{data_pkgname}{Name of the package storing the data.}
  \item{data_dirpath}{Path to the directory where the actual data is stored.}
  \item{data_serialized_objnames}{Serialized object names.}
}
\details{
The \code{GenomicScores} package defines the \code{GScores} class and associated
methods to provide support to annotation packages and resources that store genomic
scores. In the context of the \code{GenomicScores} package, genomic scores are
numeric values associated to physical nucleotide positions along a genome, such as
UCSC phastCons conservation scores for human. Two such annotation packages are:

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

}
\value{
The \code{GScores()} constructor function returns a \code{GScores-class} object.
In principle, the end-user needs not to call this function.
}
\author{R. Castelo}
\seealso{
  \code{\link[phastCons100way.UCSC.hg19]{phastCons100way.UCSC.hg19}}
  \code{\link[phastCons100way.UCSC.hg38]{phastCons100way.UCSC.hg38}}
}
\examples{
## supporting annotation packages with genomic scores
if (require(phastCons100way.UCSC.hg19)) {
  library(GenomicRanges)

  gsco <- phastCons100way.UCSC.hg19
  gsco
  scores(gsco, GRanges(seqnames="chr7", IRanges(start=117232380, width=5)))
}

## supporting AnnotationHub resources
\dontrun{
availableGScores()
gsco <- getGScores("phastCons100way.UCSC.hg19")
scores(gsco, GRanges(seqnames="chr7", IRanges(start=117232380, width=5)))
}

## meta information about a GScores object
name(gsco)
type(gsco)
provider(gsco)
providerVersion(gsco)
organism(gsco)
seqlevelsStyle(gsco)
seqinfo(gsco)
head(seqnames(gsco))
head(seqlengths(gsco))
}
\keyword{datasets}
