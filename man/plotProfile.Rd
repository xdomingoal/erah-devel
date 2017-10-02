\name{plotProfile}
\alias{plotProfile}
\title{Plotting chromatophic profile}
\description{Plots the chromatophic profiles of the compounds found by eRah.}
\usage{
plotProfile(Experiment,AlignId, per.class=T, xlim=NULL)
}
\arguments{
 \item{Experiment}{A 'MetaboSet' S4 object containing the experiment after being deconolved, aligned and (optionally) identified.}
   \item{AlignId}{the Id identificator for the compound to be shown.}
   \item{per.class}{logical. if TRUE the profiles are shown one color per class, if FALSE one color per sample.}
   \item{xlim}{x axsis (retention time) limits (see \code{\link{plot.default}}).}
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{plotSpectra}}, \code{\link{plotAlign}}
}
