\name{importMSP}
\alias{importMSP}
\title{Import MSP files to R}
\description{Import MS libraries in MSP format to eRah DB format.}
\usage{
importMSP(filename, DB.name, DB.version, DB.info)
}
\arguments{
  \item{filename}{The filepath containing the MSP library file.}
    \item{DB.name}{The name of the database (each user may chose its own name)}
    \item{DB.version}{The version of the database (each user may chose its own version)}
  \item{DB.info}{Some info about the database for further reference}
}

\details{
	
	The MSP input file should look like:

-----
	
Name: Metabolite_name 

Formula: H2O 

MW: 666

ExactMass: 666.266106 

CAS#: 11-22-3 

DB#: 1 

Comments: Metabolite_name reference standard

Num Peaks: XX

53 1; 54 2; 55 5; 56 2; 57 2; 

58 14; 59 18; 60 1000; 61 2; 67 1; 

Name: Metabolite_name_2

Formula: H2O2

MW: 999

ExactMass: 999.266106

CAS#: 22-33-4

DB#: 2

Comments: Metabolite_name_"" reference standard

Num Peaks: XX

66 10; 67 1000; 155 560; 156 800; 157 2; 

158 14; 159 1; 160 100; 161 2; 167 1; 

-------	

OR

-----
	
Name: Metabolite_name 

Formula: H2O 

MW: 666

ExactMass: 666.266106 

CASNO: 11-22-3 

DB#: 1 

Comment: Metabolite_name reference standard

Num peaks: XX

53 1

54 2

55 5 

Name: Metabolite_name_2

Formula: H2O2

MW: 999

ExactMass: 999.266106

CASNO: 22-33-4

DB#: 2

Comment: Metabolite_name_"" reference standard

Num Peaks: XX

66 10

67 1000

155 560

-------	

Or combinations of both.

	
	For more details, please see the eRah manual.
	

		
}


