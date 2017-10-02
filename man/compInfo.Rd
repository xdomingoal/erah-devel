\name{compInfo}
\alias{compInfo}
\title{Information of a Compound}
\description{Displays basic information of a compound in the MS library.}
\usage{
compInfo(comp.id, id.database=mslib)
}
\arguments{
  \item{comp.id}{The DB.Id number of the compound.}
   \item{id.database}{The mass-spectra library to be compared with the empirical spectra. By default, the MassBank - Mass Bank of North America (MoNa) database are employed (mslib object).}
}
\details{
Returns details on a given compound such as the synonyms, CAS, KEGG, retention index, among others.
}

\examples{

# finding proline
findComp("proline")

#we see that proline 2TMS has the DB.Id number 42, then:
compInfo(42)

}

\seealso{
\code{\link{findComp}}

}