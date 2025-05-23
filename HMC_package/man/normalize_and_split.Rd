% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/31_pre_processing_functions.R
\name{normalize_and_split}
\alias{normalize_and_split}
\title{Normalize and split two datasets using pooled mean and standard deviation}
\usage{
normalize_and_split(df1, df2)
}
\arguments{
\item{df1}{A data frame or matrix. Typically group 1.}

\item{df2}{A data frame or matrix. Typically group 2.}
}
\value{
A list with elements:
\describe{
\item{df1}{Normalized version of \code{df1}.}
\item{df2}{Normalized version of \code{df2}.}
}
}
\description{
Combines two datasets, normalizes features using pooled mean and standard deviation,
and returns the normalized datasets separately.
}
\examples{
set.seed(123)
df1 <- matrix(rnorm(20), nrow = 5)
df2 <- matrix(rnorm(20), nrow = 5)
normalize_and_split(df1, df2)

}
