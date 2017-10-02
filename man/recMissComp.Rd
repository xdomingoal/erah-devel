\name{recMissComp}
\alias{recMissComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Missing compound recovery}
\description{Missing compounds recovery: fits a general model (all the compounds above a certain minimum number of samples) to all the samples.}
\usage{
recMissComp(Experiment, min.samples, free.model=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data previously created by newExp, deconvolved by deconvolveComp and aligned by alignComp.}
  \item{min.samples}{The minimum number of samples in which a compound has to appear to be considered for searching into the rest of the samples where this compound missing.}
  \item{free.model}{If TRUE, the spectra found in the samples where the compound is missing is used to get the final average spectra. (See details)}

}
\details{
WARNING: If compounds were previously identified, they have to be identified again after applying the "recMissComp" function. This means that "identifyComp" function has to be executed always after "recMissComp" for identification of compounds, even if "identifyComp" has been previously applied before.

The free.model parameter is recomended to be always FALSE (except for carbon tracking applications). This is because the spectra of the samples where the compound is missing is usually affected by noise, and this could decrease the matching score for a certain compound.

}
\value{
The function returns an updated S4 'MetaboSet' class, where the GC-MS samples have been now aligned.
}

\references{
[1] Domingo-Almenara X, et al. Compound deconvolution in GC-MS-based metabolomics by blind source separation. Journal of Chromatography A (2015). Vol. 1409: 226-233. DOI: 10.1016/j.chroma.2015.07.044
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{newExp}}, \code{\link{alignComp}}, \code{\link{setAlPar}}, \code{\link{setDecPar}}
}




