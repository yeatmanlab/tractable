% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\name{read.afq.hbn}
\alias{read.afq.hbn}
\title{Load tract profiles from the Healthy Brain Network dataset}
\usage{
read.afq.hbn(truncate = FALSE, ...)
}
\arguments{
\item{truncate}{if TRUE, truncate the data to 49 rows. default = FALSE}

\item{...}{arguments to be passed to read.csv}
}
\value{
A merged dataframe with data from HBN
}
\description{
Load tract profiles from the Healthy Brain Network dataset
}
\examples{
\dontrun{
  df_hbn <- read.afq.hbn()
}
}
