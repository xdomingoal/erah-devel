\name{createdt}
\alias{createdt}
\title{Creating Experiment Tables}
\description{eRah requieres of a instrumental and (optionally) phenotype .csv file for starting/creating a new eRah project/experiment. This function automatically creates the Phenoytpe and Instrumental data .csv files.}
\usage{
createdt(path)
}
\arguments{
  \item{path}{the path where the experiment-folder is (where the experiment samples are stored).}
}
\details{
The experiment has to been organized as follows: all the samples related to each class have to be stored in the same folder (one folder = one class), and all the class-folders in one folder, which is the experiment folder. 

Two things have to be considered at this step: .csv files are different when created by American and European computers, so errors may raise due to that fact. Also, the folder containing the samples, must contain only folders. If the folder contains files (for example, already created .csv files), eRah will prompt an error.

See eRah vignette for more details. To open the vignette, execute the following code in R:
vignette("eRahManual", package="erah")

}

\examples{

# Store all the raw data files in one different folder per class,
# and all the class-folders in one folder, which is the experiment
# folder. Then execute

# createdt(path)

# where path is the experiment folder path.
# The experiment can be now startd by:

# ex <- newExp(instrumental="path/DEMO_inst.csv", 
# phenotype="path/DEMO_pheno.csv", info="DEMO Experiment")


}

\seealso{
\code{\link{newExp}}

}