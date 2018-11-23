\name{computeRIerror}
\alias{computeRIerror}
\title{Retention Index Error Computation}
\description{Computes the retention indexes (RI) errors given an internal or external calibration curve}
\usage{
computeRIerror(Experiment, id.database=mslib, reference.list, 
ri.error.type=c('relative','absolute'), plot.results=TRUE)
}
\arguments{
  	\item{Experiment}{A 'MetaboSet' S4 object containing the experiment data. The experiment has to be previously deconvolved, aligned and identified.}
  	\item{id.database}{The mass-spectra library to be compared with the empirical spectra. By default, the MassBank - Mass Bank of North America (MoNa) database are employed (mslib object).}
 	\item{reference.list}{The list containing the reference data to create the calibration curve (Internal or external calibration). Please see the Vignette and Examples.}
 	\item{ri.error.type}{Specify wether absolute or relative (default) RI error is to be computed.}
	\item{plot.results}{Shows the RI/RT graphic.}
}

\value{
Returns an 'Experiment' object. To visualize the output use: idList() function.	

}
\references{
[1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
}
\author{
Xavier Domingo-Almenara. xdomingo@scripps.edu
}

\examples{
\donttest{

#Option A: (The RT and RI of an external calibration is provided)

ex <- computeRIerror(ex, mslib, reference.list=list(RT=c(4.4683, 7.4402, 8.8121, 11.5103), 
RI=c(1081.68, 1251.31, 1346.8, 1456.8)))
id.list <- idList(ex)

#Option B: (The AlignID of internally identified metabolites by
# eRah are provided, and used as a reference. )

ex <- computeRIerror(ex, mslib, reference.list=list(AlignID = c(45,67,92,120)))
id.list <- idList(ex) 

#Please, see Vignette for extended details, usage and advices!
}}

\seealso{
\code{\link{idList}}, \code{\link{dataList}}
}