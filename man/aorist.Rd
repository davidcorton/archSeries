% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aorist.R
\name{aorist}
\alias{aorist}
\title{Calculate aoristic sum.}
\usage{
aorist(data, weight = 1, start.date = NULL, end.date = NULL,
  bin.width = 100)
}
\arguments{
\item{data}{Data table with (minimally) two numeric columns called Start and End.}

\item{weight}{Numeric vector. The weight to be applied to each row in `data`, or a constant weight to be applied to all.
Defaults to 1.}

\item{start.date}{Numeric. Start of time period to be considered. Defaults to lowest value in data$Start.}

\item{end.date}{Numeric. End of time period to be considered. Defaults to highest value in data$End.}

\item{bin.width}{Numeric. The resolution of the analysis, in units of time. Defaults to 100.}
}
\value{
data table with two named columns: 'bin', a character vector specifying the date range represented by each chronological bin;
     'aorist', a numeric vector giving the total probability mass assigned to each bin.
}
\description{
Calculates aoristic sum from a table of entities with defined date ranges, based on assumption of uniform probability between limits.
}
\examples{
date.ranges <- data.table(Start=c(450, 450, 600), End=c(700, 800, 650), frag.count=c(3, 6, 25))
x <- aorist(date.ranges, weight=date.ranges$frag.count, 500, 1500, bin.width=50)
}

