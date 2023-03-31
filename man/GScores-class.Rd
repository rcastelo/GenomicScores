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
\alias{seqinfo,GScores-method}
\alias{seqlengths,GScores-method}
\alias{seqnames,GScores-method}
\alias{seqlevelsStyle,GScores-method}
\alias{populations,GScores-method}
\alias{populations}
\alias{defaultPopulation,GScores-method}
\alias{defaultPopulation<-,GScores,character-method}
\alias{defaultPopulation<-}
\alias{defaultPopulation}
\alias{gscoresTag,GScores-method}
\alias{gscoresTag<-,GScores,character-method}
\alias{gscoresTag<-}
\alias{gscoresTag}
\alias{gscoresCategory,GScores-method}
\alias{gscoresCategory<-,GScores,character-method}
\alias{gscoresCategory<-}
\alias{gscoresCategory}
\alias{gscoresNonSNRs,GScores-method}
\alias{gscoresNonSNRs}
\alias{citation,missing-method}
\alias{citation,character-method}
\alias{citation,GScores-method}
\alias{citation}
\alias{qfun,GScores-method}
\alias{qfun}
\alias{dqfun,GScores-method}
\alias{dqfun}
\alias{nsites,GScores-method}
\alias{nsites}
\alias{show,GScores-method}

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
      genome sequence. The names of this vector correspond to the genome
      sequence names.}
    \item{\code{data_group}}{(\code{character}), name denoting a category of
      genomic scores to which the scores stored in the object belong to.
      Typical values are "Conservation", "MAF", "Pathogenicity", etc.}
    \item{\code{data_tag}}{(\code{character}), name identifying the genomic
      scores stored in the object and which can be used, for instance, to
      assign a column name storing these scores.}
    \item{\code{data_pops}}{(\code{character}), vector of character strings
      storing score population names. The term "default" is reserved to denote
      a score set that is not associated to a particular population name and
      is used by default.}
    \item{\code{data_nonsnrs}}{(\code{logical}), flag indicating whether the
      object stores genomic scores associated with non-single nucleotide ranges.}
    \item{\code{data_nsites}}{(\code{integer}), number of sites in the genome
      associated with the genomic scores stored in the object.}
    \item{\code{.data_cache}}{(\code{environment}), data structure where
      objects storing genomic scores are cached into main memory.}
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
      calculated from multiple genome alignments from the human genome version
      hg19 to 99 vertebrate species.
    }
    \item{\code{phastCons100way.UCSC.hg38}}{
      Nucleotide-level phastCons conservation scores from the UCSC Genome Browser
      calculated from multiple genome alignments from the human genome version
      hg38 to 99 vertebrate species.
    }
  }
}

\section{Constructor}{
  \describe{
    \item{\code{GScores(provider, provider_version, download_url,
                download_date, reference_genome, data_pkgname, data_dirpath,
                data_serialized_objnames, default_pop, data_tag)}:}{
      Creates a \code{GScores} object. In principle, the end-user needs not to
      call this function.
      \describe{
        \item{\code{provider}}{character string, containing the data provider.}
        \item{\code{provider_version}}{character string, containing the version
          of the data as given by the data provider.}
        \item{\code{download_url}}{character string, containing the URL of the
          data provider from where the original data were downloaded.}
        \item{\code{reference_genome}}{GenomeDescription, storing the
          information about the associated reference genome.}
        \item{\code{data_pkgname}}{character string, name given to the set of
          genomic scores stored through this object.}
        \item{\code{data_dirpath}}{character string, absolute path to the local
          directory where the genomic scores are stored.}
        \item{\code{data_serialized_objname}}{character string vector,
          containing filenames where the genomic scores are stored.}
        \item{\code{default_pop}}{character string, containing the name
          of the default scores population.}
        \item{\code{data_group}}{character string, containing a name that
          indicates a category of genomic scores to which the scores in the
          object belong to. Typical names could be "Conservation", "MAF", etc.}
        \item{\code{data_tag}}{character string, containing a tag that
          succintly labels genomic scores from a particular source. This can be
          used to automatically give, for instance, a name to a column storing
          genomic scores in data frame object. Its default value takes the
          prefix of the package name.}
      }
    }
  }
}

\section{Accessors}{
  \describe{
    \item{\code{name(x)}:}{ get the name of the set of genomic scores.
    }
    \item{\code{type(x)}:}{ get the substring of the name of the set of genomic
      scores comprised between the first character until the first period. This
      should typically match the type of genomic scores such as,
      \code{phastCons}, \code{phyloP}, etc.
    }
    \item{\code{provider(x)}:}{ get the data provider.
    }
    \item{\code{providerVersion(x)}:}{ get the provider version.
    }
    \item{\code{organism(x)}:}{ get the organism associated with the genomic scores.
    }
    \item{\code{seqlevelsStyle(x)}:}{ get the genome sequence style.
    }
    \item{\code{seqinfo(x)}:}{ get the genome sequence information.
    }
    \item{\code{seqnames(x)}:}{ get the genome sequence names.
    }
    \item{\code{seqlengths(x)}:}{ get the genome sequence lengths.
    }
    \item{\code{populations(x)}:}{ get the identifiers of the available scores
      populations. If only one scores population is available, then it shows
      only the term \code{default}.
    }
    \item{\code{defaultPopulation(x)}:}{ get or set the default population of scores.
    }
    \item{\code{gscoresCategory(x)}:}{ get or set the genomic scores category label.
    }
    \item{\code{gscoresTag(x)}:}{ get or set the genomic scores tag label.
    }
    \item{\code{gscoresNonSNRs(x)}:}{ get whether there are genomic scores associated
      with non-single nucleotide ranges.
    }
    \item{\code{nsites(x)}:}{ get the number of sites in the genome with genomic
      scores.
    }
    \item{\code{qfun(x)}:}{ get the quantizer function.
    }
    \item{\code{dqfun(x)}:}{ get the dequantizer function.
    }
    \item{\code{citation(x)}:}{ get citation information for the genomic scores data
      in the form of a \code{bibentry} object.
    }
  }
}

\references{
Puigdevall, P. and Castelo, R. GenomicScores: seamless access to genomewide
position-specific scores from R and Bioconductor.
\emph{Bioinformatics}, 18:3208-3210, 2018.
}
\author{R. Castelo}
\seealso{
  \code{\link{gscores}()}
  \code{\link{score}()}
  \code{\link[phastCons100way.UCSC.hg19]{phastCons100way.UCSC.hg19}}
}
\examples{
## one genomic range of width 5
gr1 <- GRanges(seqnames="chr7", IRanges(start=117232380, width=5))
gr1

## five genomic ranges of width 1
gr2 <- GRanges(seqnames="chr7", IRanges(start=117232380:117232384, width=1))
gr2

## supporting annotation packages with genomic scores
if (require(phastCons100way.UCSC.hg19)) {
  library(GenomicRanges)

  gsco <- phastCons100way.UCSC.hg19
  gsco
  gscores(gsco, gr1)
  score(gsco, gr1)
  gscores(gsco, gr2)
  populations(gsco)
  gscores(gsco, gr2, pop="DP2")
}

## supporting AnnotationHub resources
\dontrun{
availableGScores()
gsco <- getGScores("phastCons100way.UCSC.hg19")
gsco
gscores(gsco, gr1)
}

## metadata from a GScores object
name(gsco)
type(gsco)
provider(gsco)
providerVersion(gsco)
organism(gsco)
seqlevelsStyle(gsco)
seqinfo(gsco)
head(seqnames(gsco))
head(seqlengths(gsco))
gscoresTag(gsco)
populations(gsco)
defaultPopulation(gsco)
qfun(gsco)
dqfun(gsco)
citation(gsco)
}
\keyword{methods}
