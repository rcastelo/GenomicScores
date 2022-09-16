\name{gscores}

\alias{gscores,GScores,GenomicRanges-method}
\alias{gscores,GScores,character-method}
\alias{gscores}
\alias{score,GScores-method}

\title{Accessing genomic gscores}

\description{
Functions to access genomic gscores through \code{GScores} objects.
}
\usage{
\S4method{gscores}{GScores,GenomicRanges}(x, ranges, ...)
\S4method{gscores}{GScores,character}(x, ranges, ...)
\S4method{score}{GScores}(x, ..., simplify=TRUE)
}
\arguments{
  \item{x}{A \code{GScores} object.}
  \item{ranges}{A \code{GenomicRanges} object with positions from where to
                retrieve genomic scores, or a character string vector with
                identifiers associated by the data producer with the genomic
                scores, e.g., dbSNP 'rs' identifiers in the case of some
                \code{MafDb.*} packages.}
  \item{...}{In the call to the \code{gscores}() method one can additionally
             specify the following optional parameters:
             \itemize{
               \item{pop: }{Character string vector specifying the scores
                            populations to query, when there is more than one.
                            By default, its value is 'defaultPopulation(x)'.
                            Use \code{\link{populations}()} to find out the
                            available scores populations.}
               \item{type: }{Character string specifying the type of genomic
                             position being sought, which can be a single
                             nucleotide range (\code{snr}), by default, or a
                             \code{nonsnr} spanning multiple nucleotides. The
                             latter is the case of indel variants in minor
                             allele frequency data.}
               \item{scores.only: }{Flag setting whether only the scores should
                                    be returned as a numeric vector (\code{TRUE}),
                                    instead of returning them as a metadata column
                                    in a \code{GRanges} object (\code{FALSE}, default).}
               \item{summaryFun: }{Function to summarize genomic scores when
                                   more than one position is retrieved. By
                                   default, this is set to the arithmetic mean,
                                   i.e., the \code{mean()} function.}
               \item{quantized: }{Flag setting whether the genomic scores
                                  should be returned quantized (\code{TRUE})
                                  or dequantized (\code{FALSE}, default).}
               \item{ref: }{Vector of reference alleles in the form of either
                            a character vector, a \code{DNAStringSet} object
                            or a \code{DNAStringSetList} object. This argument
                            is used only when either there are multiple scores per
                            position or \code{x} is a \code{MafDb.*} package.}
               \item{alt: }{Vector of alternative alleles in the form of either
                            a character vector, a \code{DNAStringSet} object or
                            a \code{DNAStringSetList} object. This argument is
                            used only when either there are multiple scores per
                            position or \code{x} is a \code{MafDb.*} package.}
               \item{minoverlap: }{Integer value passed internally to the
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
               \item{caching: }{Flag setting whether genomic scores per
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
The method \code{gscores()} takes as first argument a \code{GScores}
object, previouly loaded from either an annotation package or an
\code{AnnotationHub} resource; see \code{\link{getGScores}()}.

The arguments \code{ref} and \code{alt} serve two purposes. One, when there
are multiple scores per position, such as with CADD or M-CAP, and we want to
select a score matching a specific combination of reference and alternate
alleles. The other purpose is when the \code{GScores} object \code{x} is a
\code{MafDb.*} package, then by providing \code{ref} and \code{alt} alelles
we will get separate frequencies for reference and alternate alleles. The
current lossy compression of these values yields a correct assignment for
biallelic variants in the corresponding \code{MafDb.*} package and an
approximation for multiallelic ones.
}
\value{
The method \code{gscores()} returns a \code{GRanges} object with the genomic
scores in a metadata column called \code{score}. The method \code{score()}
returns a numeric vector with the genomic scores.
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
}
\keyword{utilities}
