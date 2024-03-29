% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hreal.R
\name{hreal}
\alias{hreal}
\alias{print.hreal}
\alias{summary.hreal}
\title{Realization of Hawkes process}
\usage{
\method{print}{hreal}(x, n = 20, ...)

\method{summary}{hreal}(object, n = 20, ...)
}
\arguments{
\item{x}{S3-object of \code{hreal}.}

\item{n}{Number of rows to display.}

\item{...}{Further arguments passed to or from other methods.}

\item{object}{S3-object of \code{hreal}.}
}
\description{
\code{hreal} is the list of the following:
\itemize{
\item \code{hspec} : S4 object \code{\link{hspec-class}} that specifies the parameter values.
\item \code{inter_arrival} : the time between two consecutive events.
\item \code{arrival} : cumulative sum of \code{inter_arrival}.
\item \code{type} : integer, the type of event.
\item \code{mark} : the size of mark, an additional information associated with event.
\item \code{N} : counting process that counts the number of events.
\item \code{Nc} : counting process that counts the number of events weighted by mark.
\item \code{lambda} : intensity process, left-continuous version.
\item \code{lambda_component} : the component of intensity process with \code{mu} not included.
\item \code{rambda} : intensity process, right-continuous version.
\item \code{rambda_component} : the right-continuous version of \code{lambda_component}.
}

Print functions for \code{hreal} are provided.
}
