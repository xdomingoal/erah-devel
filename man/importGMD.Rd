\name{importGMD}
\alias{importGMD}
\title{Import MSP files from GMD to R}
\description{Import the Golm Metabolome Database.}
\usage{
importGMD(filename, DB.name, DB.version, DB.info, 
type=c("VAR5.ALK","VAR5.FAME","MDN35.ALK", "MDN35.FAME"))
}
\arguments{
  \item{filename}{The filepath containing the GMD database file.}
    \item{DB.name}{The name of the database (each user may chose its own name)}
    \item{DB.version}{The version of the database (each user may chose its own version)}
  \item{DB.info}{Some info about the database for further reference}
  \item{type}{The type of RI to be imported from the database}
}

\details{
	For more details, please see the eRah manual	
}


