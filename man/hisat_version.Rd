% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hisat.R
\name{hisat_version}
\alias{hisat_version}
\title{Print version information of hisat-0.1.6-beta}
\usage{
hisat_version()
}
\value{
An invisible \code{Integer} of call status.
The value is 0 when there is not any mistakes
Otherwise the value is non-zero.
}
\description{
Print version information of hisat-0.1.6-beta
}
\examples{
cmdout<-hisat_version();cmdout
}
\references{
Kim, D., Langmead, B. & Salzberg, S. HISAT: a fast 
spliced aligner with low memory requirements. Nat Methods 12, 
357-360 (2015).
}
\author{
Qixiu Du
}
