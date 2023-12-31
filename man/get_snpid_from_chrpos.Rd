% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.r
\name{get_snpid_from_chrpos}
\alias{get_snpid_from_chrpos}
\title{Get variant information based on chr/pos lookup}
\usage{
get_snpid_from_chrpos(
  chrpos,
  build = c("19", "38")[1],
  map_files = get_mapfiles()
)
}
\arguments{
\item{chrpos}{vector of variants by chr/pos e.g. \code{c("10:60494", "10:60515", "10:60523")}}

\item{build}{build of chrpos being queries - 19 or 38}

\item{map_files}{path to the \code{Metadata/SNP RSID maps} directory}
}
\value{
data frame
}
\description{
Provide a path to the \code{Metadata/SNP RSID maps} directory and lookup the snp IDs by chromosome position
}
