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
  \item{sites}{
%%     ~~Describe \code{x} here~~
}
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
