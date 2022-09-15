\name{wgscores}

\alias{wgscores,GScores,GenomicRanges-method}
\alias{wgscores}

\title{Discovering genomic gscores}

\description{
Functions to discover which genomic scores are present in given genomic ranges through \code{GScores} objects.
}
\usage{
\S4method{wgscores}{GScores,GenomicRanges}(x, ranges, ...)
}
\arguments{
  \item{x}{A \code{GScores} object.}
  \item{ranges}{A \code{GenomicRanges} object with genomic ranges where to
                search for genomic scores.}
  \item{...}{In the call to the \code{wgscores}() method one can additionally
             set the following arguments:
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
               \item{caching: }{Flag setting whether genomic scores per
                                chromosome should be kept cached in memory
                                (\code{TRUE}, default) or not (\code{FALSE}).
                                The latter option minimizes the memory
                                footprint but slows down the performance when
                                the \code{wgscores()} method is called multiple
                                times.}
             }}
}
\details{
The method \code{wgscores()} takes as first argument a \code{GScores}
object, previouly loaded from either an annotation package or an
\code{AnnotationHub} resource; see \code{\link{getGScores}()} and a
\code{GenomicRanges} object as a second argument. It will search for
which genomic scores fall within the provided genomic ranges and return
them in an ordered \code{GenomicRanges} object with the scores in the
metadata columns.
}
\value{
The method \code{wgscores()} returns a \code{GRanges} object with the
genomic scores in metadata columns named after the corresponding score
population name.
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
if (require(MafDb.1Kgenomes.phase1.hs37d5)) {
  mafdb <- MafDb.1Kgenomes.phase1.hs37d5
  mafdb
  populations(mafdb)

  ## lookup allele frequencies for SNP rs1129038, located at 15:28356859 in
  ## GRCh37, a SNP associated to blue and brown eye colors as reported by
  ## Eiberg et al. (Human Genetics, 2008; http://www.ncbi.nlm.nih.gov/pubmed/18172690).
  ## Blue eye color in humans may be caused by a perfectly associated founder
  ## mutation in a regulatory element located within the HERC2 gene inhibiting
  ## OCA2 expression.
  ##
  ## for the sake of illustrating this functionality let's create a
  ## GenomicRanges object with the SNP rs1129038 and enlarge its range
  ## by 200nt at each flank.
  snp <- GRanges("15:28356859")
  rngsnp <- flank(snp, width=100, both=TRUE)
  width(rngsnp)

  ## now use this genomic range to search for the rs1129038 SNP
  wgscores(mafdb, rngsnp)

  ## let's illustrate this same functionality for INDELs, starting
  ## from the specific INDEL rs113993960 that leads to the loss of
  ## phenylalanine at amino acid position 508 of the CFTR protein,
  ## commonly referred to as F508del in the CFTR gene, which is
  ## concretely a deletion of four nucleotides at position
  ## 7:117199644 in GRCh37 and enlarge its range by 20nt on each flank.
  indel <- GRanges(seqnames="chr7",
                   ranges=IRanges(start=117199644, width=4))
  rngindel <- flank(indel, width=20, both=TRUE)
  width(rngindel)

  ## now use this genomic range to search for the rs113993960 INDEL
  wgscores(mafdb, rngindel, type="nonsnrs")
}
}
\keyword{utilities}
