\name{deconvolveComp}
\alias{deconvolveComp}
\title{Deconvolution of compounds in samples}
\description{Deconvolution of GC-MS data}
\usage{
deconvolveComp(Experiment, decParameters, 
samples.to.process=NULL, down.sample=FALSE,
virtual.scans.ps=NULL)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment data previously created by newExp.}
  \item{decParameters}{The software deconvolution parameters object previously created by setDecPar}
  \item{samples.to.process}{Vector indicating which samples are to be processed.} 
  \item{down.sample}{If TRUE, chromatograms are down sampled to define one peak with 10 scan points (according to the minimum peak width). This is to process longer chromatograms with wider peak widths (more than 20 seconds peak width and small scans per second values). See details. }
   \item{virtual.scans.ps}{Manually correction of scans per second. When chromatograms are downsampled (too few scans per second, or too many), a virtual scans per second can be defined, and data interpolation will correct the data. However, reanalysis of experimental data is advised.}  
}
\details{
See eRah vignette for more details. To open the vignette, execute the following code in R:
vignette("eRahManual", package="erah")

eRah uses multivariate methods which run-time performance depend on the amount of data to be analyzed. When peaks are wide and the scans per second used for acquisition is too large, the number of points (scans) that define a peak might be too many, leading eRah to a poor run-time performance. To solve that, use down.sample=TRUE to allow eRah to define a peak with 10 seconds and analyze the data more efficiently
}
\value{
The function returns an updated S4 'MetaboSet' class, where the GC-MS samples have been now deconvolved.
}
\references{
[1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\seealso{
\code{\link{newExp}}, \code{\link{setAlPar}}
}
\examples{
\donttest{
# Deconvolve data from a created experiment by newExp().
ex <- newExp(instrumental="path")

# The following will set eRah for analyzing the chromatograms
# from minutes 5 to 15, and withouth taking into account the masses
# 35:69,73:75,147:149, with a minimum peak widht of 0.7 seconds.

ex.dec.par <- setDecPar(min.peak.width=0.7, min.peak.height=5000, 
noise.threshold=500, avoid.processing.mz=c(35:69,73:75,147:149), 
analysis.time=c(5,15))

# An now deconvolve the compounds in the samples:
ex <- deconvolveComp(ex, decParameters=ex.dec.par)
}
}


