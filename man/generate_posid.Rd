\name{generate_posid}
\alias{generate_posid}
\title{Generate position IDs}
\description{
Given a set of genomic integration positions, the function generates a character
vector of integrated positions in the following format: chr(+/-)pos, where chr
is the chromosome, (+/-) designates the orientation of the integration, and pos
is the position of the integrated element.
}
\usage{
generate_posid(sites = NULL)

generate_posid(seqnames = NULL, strand = NULL, start = NULL, end = NULL)
}
\arguments{
  \item{sites}{: a GRanges object where each row represents one integrated element.}
  \item{seqnames}{: a character vector of chromosomes for integrated elements, i.e. chr1, chr3, chrX.}
  \item{strand}{: orientation or strand on which the integrated element is on (+, -, *)}
  \item{start}{: the lower numerical position of the integrated element, considered the start for "+" integrated elements.}
  \item{end}{: the greater numerical position of the integrated element, considered the end for "+" integrated elements.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Christopher Nobles, Ph.D.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
chr <- c("chr1", "chr3", "chrX")
strands <- c("+", "-", "+")
starts <- c(900231, 13254892, 603292)
ends <- c(900431, 13255292, 603592)
ranges <- IRanges(start = starts, end = ends)
gr <- GRanges(seqnames = chr, ranges = ranges, strand = strands)

generate_posid(sites = gr)

generate_posid(seqnames = chr, strand = strands, start = starts, end = ends)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
