\name{gscores}

\alias{gscores,GScores,GenomicRanges-method}
\alias{gscores,GScores,character-method}
\alias{gscores,MafDb,GenomicRanges-method}
\alias{gscores}
\alias{score,GScores-method}
\alias{score,MafDb-method}
\alias{availableGScores}
\alias{getGScores}

\title{Accessing genomic gscores}

\description{
Functions to access genomic gscores through \code{GScores} objects.
}
\usage{
availableGScores()
getGScores(x)
\S4method{gscores}{GScores,GenomicRanges}(x, ranges, ...)
\S4method{gscores}{GScores,character}(x, ranges, ...)
\S4method{gscores}{MafDb,GenomicRanges}(x, ranges, ...)
\S4method{score}{GScores}(x, ..., simplify=TRUE)
\S4method{score}{MafDb}(x, ..., simplify=TRUE)
}
\arguments{
  \item{x}{For \code{getGScores()}, a \code{character} vector of length 1
           specifiying the genomic scores resource to fetch. For
           \code{gscores()} and \code{score()}, a \code{GScores} object.}
  \item{ranges}{A \code{GenomicRanges} object with positions from where to
                retrieve genomic scores, or a character string vector with
                identifiers associated by the data producer to the genomic
                scores, e.g., dbSNP 'rs' identifiers.}
  \item{...}{In the call to the \code{gscores}() method one can additionally
             set the following arguments:
             \itemize{
               \item{pop}{Character string vector specifying the scores
                          populations to query, when there is more than one.
                          Use \code{\link{populations}()} to find out the
                          available scores populations.}
               \item{type}{Character string specifying the type of genomic
                           position being sought, which can be a single
                           nucleotide range (\code{snr}), by default, or a
                           \code{nonsnr} spanning multiple nucleotides. The
                           latter is the case of indel variants in minor
                           allele frequency data.}
               \item{summaryFun}{Function to summarize genomic scores when
                                 more than one position is retrieved. By
                                 default, this is set to the arithmetic mean,
                                 i.e., the \code{mean()} function.}
               \item{quantized}{Flag setting whether the genomic scores
                                should be returned quantized (\code{TRUE})
                                or dequantized (\code{FALSE}, default).}
               \item{ref}{Vector of reference alleles in the form of either
                          a character vector, a \code{DNAStringSet} object
                          or a \code{DNAStringSetList} object. This argument
                          is used only when there are multiple scores per
                          position.}
               \item{alt}{Vector of alternative alleles in the form of either
                          a character vector, a \code{DNAStringSet} object or
                          a \code{DNAStringSetList} object. This argument is
                          used only when there are multiple scores per
                          position.}
               \item{minoverlap}{Integer value passed internally to the
                                 function \code{findOverlaps()} from the
                                 \code{IRanges} package, when querying
                                 genomic positions associated with
                                 multiple-nucleotide ranges (nonSNRs). By
                                 default, \code{minoverlap=1L}, which assumes
                                 that the sought nonSNRs are stored as in
                                 \code{VCF} files, using the nucleotide
                                 composition of the reference sequence. This
                                 argument is only relevant for genomic scores
                                 associated with nonSNRs.}
               \item{caching}{Flag setting whether genomic scores per
                              chromosome should be kept cached in memory
                              (\code{TRUE}, default) or not (\code{FALSE}).
                              The latter option minimizes the memory
                              footprint but slows down the performance when
                              the \code{gscores()} method is called multiple
                              times.}
             }}
  \item{simplify}{Flag setting whether the result should be simplified to a
                  vector (\code{TRUE}, default) if possible. This happens
                  when scores from a single population are queried.}
}
\details{
The function \code{availableGScores()} shows genomic score sets available as
\code{AnnotationHub} online resources.

The method \code{gscores()} takes as first argument a \code{GScores-class}
object that can be loaded from an annotation package or from an
\code{AnnotationHub} resource. These two possibilities are illustrated in the
examples below.
}
\value{
The function \code{availableGScores()} returns a character vector with the
names of the \code{AnnotationHub} resources corresponding to different
available sets of genomic scores. The function \code{getGScores()} return a
\code{GScores} object. The method \code{gscores()} returns a \code{GRanges}
object with the genomic scores in a metadata column called \code{score}.
The method \code{score()} returns a numeric vector with the genomic scores.
}
\author{R. Castelo}
\seealso{
  \code{\link[phastCons100way.UCSC.hg19]{phastCons100way.UCSC.hg19}}
  \code{\link[MafDb.1Kgenomes.phase1.hs37d5]{MafDb.1Kgenomes.phase1.hs37d5}}
}
\examples{
## one genomic range of width 5
gr1 <- GRanges(seqnames="chr7", IRanges(start=117232380, width=5))
gr1

## five genomic ranges of width 1
gr2 <- GRanges(seqnames="chr7", IRanges(start=117232380:117232384, width=1))
gr2

## accessing genomic gscores from an annotation package
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

if (require(MafDb.1Kgenomes.phase1.hs37d5)) {
  mafdb <- MafDb.1Kgenomes.phase1.hs37d5
  mafdb
  populations(mafdb)

  ## lookup allele frequencies for SNP rs1129038, located at 15:28356859, a
  ## SNP associated to blue and brown eye colors as reported by Eiberg et al.
  ## Blue eye color in humans may be caused by a perfectly associated founder
  ## mutation in a regulatory element located within the HERC2 gene
  ## inhibiting OCA2 expression. Human Genetics, 123(2):177-87, 2008
  ## [http://www.ncbi.nlm.nih.gov/pubmed/18172690]
  gscores(mafdb, GRanges("15:28356859"), pop=populations(mafdb))
  gscores(mafdb, "rs1129038", pop=populations(mafdb))
}

## accessing genomic scores from AnnotationHub resources
\dontrun{
availableGScores()
gsco <- getGScores("phastCons100way.UCSC.hg19")
gscores(gsco, gr1)
}
}
\keyword{utilities}
