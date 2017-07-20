\name{GScores-class}
\docType{class}

% Class:
\alias{class:GScores}
\alias{GScores-class}
\alias{GScores}
\alias{GenomicScores}

% Accessors:
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
\alias{citation,missing-method}
\alias{citation,character-method}
\alias{citation,GScores-method}
\alias{citation}
\alias{qfun,GScores-method}
\alias{qfun}
\alias{dqfun,GScores-method}
\alias{dqfun}
\alias{show,GScores-method}
\alias{makeGScoresPackage}

\title{GScores objects}

\description{
  The goal of the \code{GenomicScores} package is to provide support to store
  and retrieve genomic scores associated to physical nucleotide positions along
  a genome. This is achieved through the \code{GScores} class of objects, which
  is a container for genomic score values.
}
 
\details{
  The \code{GScores} class attempts to provide a compact storage and efficient
  retrieval of genomic score values that have been typically processed and
  stored using some form of lossy compression. This class is currently based
  on a former version of the \code{SNPlocs} class defined in the
  \code{BSgenome} package, with the following slots:
  \describe{
    \item{\code{provider}}{(\code{character}), the data provider such as UCSC.}
    \item{\code{provider_version}}{(\code{character}), the version of the data
      as given by the data provider, typically a date in some compact format.}
    \item{\code{download_url}}{(\code{character}), the URL of the data provider
      from where the original data were downloaded.}
    \item{\code{download_date}}{(\code{character}), the date on which the data
      were downloaded.}
    \item{\code{reference_genome}}{(\code{GenomeDescription}), object with
      information about the reference genome whose physical positions have
      the genomic scores.}
    \item{\code{data_pkgname}}{(\code{character}), name given to the set
      of genomic scores associated to a particular genome. When the genomic
      scores are stored within an annotation package, then this corresponds to
      the name of that package.}
    \item{\code{data_dirpath}}{(\code{character}), absolute path to the local
      directory where the genomic scores are stored in one file per genome
      sequence.}
    \item{\code{data_serialized_objnames}}{(\code{character}), named vector of
      filenames pointing to files containing the genomic scores in one file per
      genome sequence. The names of this vector correspond to the genome sequence
      names.}
    \item{\code{.data_cache}}{(\code{environment}), data structure where objects
      storing genomic scores are cached into main memory.}
  }
  The goal of the design behind the \code{GScores} class is to load into main
  memory only the objects associated with the queried sequences to minimize the
  memory footprint, which may be advantageous in workflows that parallelize the
  access to genomic scores by genome sequence.

  \code{GScores} objects are created either from \code{AnnotationHub} resources
  or when loading specific annotation packages that store genomic score values.
  Two such annotation packages are:

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

\section{Constructor}{
  \describe{
    \item{}{
      \code{GScores(provider, provider_version, download_url,
        download_date, reference_genome, data_pkgname, data_dirpath,
        data_serialized_objnames)}:
      Creates a \code{GScores} object. In principle, the end-user needs not to call
      this function.
      \describe{
        \item{\code{provider}}{character, containing the data provider.}
        \item{\code{provider_version}}{character, containing the version of
          the data as given by the data provider.}
        \item{\code{download_url}}{character, containing the URL of the
          data provider from where the original data were downloaded.}
        \item{\code{reference_genome}}{GenomeDescription, storing the
          information about the associated reference genome.}
        \item{\code{data_pkgname}}{character, name given to the set of genomic
          scores stored through this object.}
        \item{\code{data_dirpath}}{character, absolute path to the local
          directory where the genomic scores are stored.}
        \item{\code{data_serialized_objname}}{character vector, containing
          filenames where the genomic scores are stored.}
      }
    }
  }
}

\section{Accessors}{
  \describe{
    \item{}{
      \code{name(x)}: get the name of the set of genomic scores.
    }
    \item{}{
      \code{type(x)}: get the substring of the name of the set of genomic scores
      comprised between the first character until the first period. This should
      typically match the type of genomic scores such as, \code{phastCons},
      \code{phyloP}, etc.
    }
    \item{}{
      \code{provider(x)}: get the data provider.
    }
    \item{}{
      \code{providerVersion(x)}: get the provider version.
    }
    \item{}{
      \code{organism(x)}: get the organism associated with the genomic scores.
    }
    \item{}{
      \code{referenceGenome(x)}: get the \code{GenomeDescription} object
      associated with the genome on which the genomic scores are defined.
    }
    \item{}{
      \code{seqlevelsStyle(x)}: get the genome sequence style.
    }
    \item{}{
      \code{seqinfo(x)}: get the genome sequence information.
    }
    \item{}{
      \code{seqnames(x)}: get the genome sequence names.
    }
    \item{}{
      \code{seqlengths(x)}: get the genome sequence lengths.
    }
    \item{}{
      \code{qfun(x)}: get the quantizer function.
    }
    \item{}{
      \code{dqfun(x)}: get the dequantizer function.
    }
    \item{}{
      \code{citation(x)}: get citation information for the genomic scores data
      in the form of a \code{bibentry} object.
    }
  }
}

\author{R. Castelo}
\seealso{
  \code{\link{scores}()}
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
gsco
scores(gsco, GRanges(seqnames="chr7", IRanges(start=117232380, width=5)))
}

## meta data from a GScores object
name(gsco)
type(gsco)
provider(gsco)
providerVersion(gsco)
organism(gsco)
referenceGenome(gsco)
seqlevelsStyle(gsco)
seqinfo(gsco)
head(seqnames(gsco))
head(seqlengths(gsco))
qfun(gsco)
dqfun(gsco)
citation(gsco)
}
\keyword{datasets}
