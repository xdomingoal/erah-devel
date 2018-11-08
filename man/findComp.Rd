\name{findComp}
\alias{findComp}
\title{Find a Compound}
\description{Finds compounds in the MS library by Name, CAS or chemical formula.}
\usage{
findComp(name = NULL, id.database = mslib, CAS = NULL, 
chem.form = NULL) 
}
\arguments{
  \item{name}{The name of the compound to be found.}
    \item{id.database}{The mass-spectra library to be compared with the empirical spectra. By default, the MassBank - Mass Bank of North America (MoNa) database are employed (mslib object).}
  \item{CAS}{The CAS number of the compound to be found.}
  \item{chem.form}{The chemical formula of the compound to be found.}
}

\value{
\code{findComp} returns an S3 object:
	\item{DB.Id}{The identification number of the library. Each metbolite in the reference library has a different DB.Id number.}
	\item{Compound Name}{Compound Name.}
	\item{CAS}{CAS number}
	\item{Formula}{Chemical Formula.}
}

\examples{

# finding proline

findComp("proline")

# be careful, exact matches are not supported, 
# as well as special cases of partial searches, like these cases:

findComp("L-proline (2TMS)")

findComp("proline 2")

# which will not report any results despite being in the database
}

\seealso{
\code{\link{compInfo}}

}