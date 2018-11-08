\name{newExp}
\alias{newExp}
\title{New Experiment}
\description{Sets a new experiment for eRah}
\usage{
newExp(instrumental, phenotype=NULL, info=character())
}
\arguments{
  \item{instrumental}{The path where the instrumental .csv file is located.}
  \item{phenotype}{(optional) The path where the phenotypic .csv file is located.}
  \item{info}{Experiment description}
}
\details{
See eRah vignette for more details. To open the vignette, execute the following code in R:
vignette("eRahManual", package="erah")

}
\value{
\code{newExp} returns an S4 object of the class 'MetaboSet'.
}
\references{
[1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
}
\author{
Xavier Domingo-Almenara. xavier.domingo@urv.cat
}

\examples{
\donttest{
# Store all the raw data files in one different folder per class,
# and all the class-folders in one folder, which is the experiment
# folder. Then execute

createdt(path)

# where path is the experiment folder path.
# The experiment can be now started by:

ex <- newExp(instrumental="path/DEMO_inst.csv", 
phenotype="path/DEMO_pheno.csv", info="DEMO Experiment")


}}

\seealso{
\code{\link{createdt}}, \code{\link{setDecPar}}, \code{\link{setAlPar}} 
}

