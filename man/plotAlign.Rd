% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotAlign}
\alias{plotAlign}
\alias{plotAlign,MetaboSet-method}
\title{Plotting chromatographic profile with and without alignment}
\usage{
plotAlign(Experiment,AlignId, per.class = T, xlim = NULL)

\S4method{plotAlign}{MetaboSet}(Experiment, AlignId, per.class = T, xlim = NULL)
}
\arguments{
\item{Experiment}{A 'MetaboSet' S4 object containing the experiment after being deconolved, aligned and (optionally) identified.}

\item{AlignId}{the Id identificator for the compound to be shown.}

\item{per.class}{logical. if TRUE the profiles are shown one color per class, if FALSE one color per sample.}

\item{xlim}{x axsis (retention time) limits (see \code{\link{plot.default}}).}
}
\description{
Plots the chromatophic profiles of the compounds found by eRah. Similarly to plotProfile, but with two sub-windows, showing the chromatophic profiles before and after alignment.
}
\seealso{
\code{\link{plotSpectra}} \code{\link{plotProfile}}
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}
