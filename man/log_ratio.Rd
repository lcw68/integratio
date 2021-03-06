% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{log_ratio}
\alias{log_ratio}
\title{Two version ratio calculation}
\usage{
log_ratio(species_data, version = 1, epsilon = 1)
}
\arguments{
\item{species_data}{p*2 matrix, p represents the samples, DNA and RNA columns}

\item{version}{which ratio version should be calculated}

\item{epsilon}{the constant we add in two versions ratio}
}
\value{
A data frame contaning ratio and which type it belongs to for every sample
}
\description{
Two version ratio calculation
}
\examples{

log_ratio(otu[1,,],version = 1, epsilon = 1)




}
