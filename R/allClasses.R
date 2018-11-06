## MetaboSet Class Definition:

setClass(Class = "eRahSoftParameters", 
         representation = representation(algorithm = "character", 
                                         min.peak.width = "numeric", 
                                         min.peak.height = "numeric", 
                                         noise.threshold = "numeric", 
                                         avoid.processing.mz = "vector",
                                         compression.coef = "numeric", 
                                         analysis.time = "vector"))

setClass(Class = "eRahAlParameters", 
         representation = representation(algorithm = "character", 
                                         min.spectra.cor = "numeric", 
                                         max.time.dist = "numeric", 
                                         mz.range = "vector", 
                                         method = "character")
)	

setClass(Class = "eRahIdParameters", 
         representation = representation(algorithm = "character", 
                                         database = "character", 
                                         compare.only.mz = "vector")
)

setClass(Class = "MetaData", 
         representation = representation(Instrumental = "data.frame", 
                                         Phenotype = "data.frame", 
                                         DataDirectory = "character")
)

setClass(Class = "Statistics", 
         representation = representation(Univariate = "data.frame", 
                                         Multivariate = "data.frame")
)	

setClass(Class = "MSResultsParameters", 
         representation = representation(Alignment = "list", 
                                         Identification = "list")
)

setClass(Class = "Data", 
         representation = representation(FeatureList = "list",
                                         FactorList = "list", 
                                         Parameters = "list")
)

setClass(Class = "Results", 
         representation = representation(Parameters = "MSResultsParameters", 
                                         Alignment = "data.frame", 
                                         Identification = "data.frame", 
                                         Statistics = "Statistics"))

#' @name MetaboSet-class
#' @docType class
#' @aliases MetaboSet-class
#' @title Class \code{"MetaboSet"}
#' @description The MetaboSet class is a single generic class valid for all sorts of metabolomic studies regardless of the experimental platform, the statistical processing and the annotation stage. It is the core operation class of eRah.
#' @slot Info Slot Info stores the general information of the experiment and the experimental platform used in the analysis of the biological samples.
#' @slot Data Slot Data contains either the raw data or the path of the files. It also contains the list of the selected features (deconvolved compounds). In the subslot Parameters it is saved the information regarding the feature selector algorithm (type, parameters, version...) and the experimental platform used.
#' @slot MetaData Slot MetaData has two slots. In the Instrumental slot it is saved a data frame with some mandatory fields (filename, date, time, sampleID) and optional fields related to the experimental platform (Column ID, Column Type, Ioniser,...). Slot Phenotypic contains a data frame with the sample and experimental information (phenotypes, longitudinal data,...).
#' @slot Results In the Results slot it is saved the information related to the statistical and identification results. The slot Parameters contains all the values of the parameters used in the identification and statistical functions. Slot Identification has the results of the identification process as well as the identification or/and annotation steps. The results of the statistical functions are saved in the Statistics slot.
#' @author Xavier Domingo-Almenara, Arnald Alonso and Francesc Fernandez-Albert.

setClass(Class = "MetaboSet",
         representation = representation(Info = "character", 
                                         Data = "Data", 
                                         MetaData = "MetaData", 
                                         Results = "Results")
)

## Intern Classes for eRah:

#' @name eRah_DB-class
#' @docType class
#' @aliases eRah_DB-class
#' @title Class \code{"eRah_DB"}
#' @description The eRah_DB class contains the slots for storing and accessing a MS library.
#' @slot name The name of the stored library
#' @slot version The version of the stored library (and which is the database identifier, should be unique and used to check if is the database used in other experiments)
#' @slot info Character vector containing complementary information about the library.
#' @slot database A list of S3 objects, which each object contains the information on a different compound.
#' @author Xavier Domingo-Almenara.

setClass(Class = "eRah_DB", 
         representation = representation(name = "character", 
                                         version = "character", 
                                         info = "character", 
                                         database = "list")
)

#'  @name RawDataParameters-class
#'  @docType class
#'  @aliases RawDataParameters-class
#'  @title Class \code{"RawDataParameteres"}
#'  @description The RawDataParameters class contains the slots for storing and accessing into a MS sample, and the essential parameters for performing its processing (deconvolution).
#'  @slot data The data matrix of the sample to be processed
#'  @slot min.mz The minimum adquired mz number
#'  @slot max.mz The maximum adquired mz number
#'  @slot start.time Starting time of adquisition
#'  @slot mz.resolution Mz resolution
#'  @slot scans.per.second Scans per second
#'  @slot avoid.processing.mz Which mz do not have to be processed
#'  @slot min.peak.width Minimum peak width (stored in scans)
#'  @slot min.peak.height Minimum peak height
#'  @slot noise.threshold The noise threshold
#'  @slot compression.coef Compression coefficient (parameter for Orthogonal Signal Deconvolution)
#'  @author Xavier Domingo-Almenara.

setClass(Class = "RawDataParameters", 
         representation = representation(data = "matrix", 
                                         min.mz = "numeric", 
                                         max.mz = "numeric", 
                                         start.time = "numeric", 
                                         mz.resolution = "numeric", 
                                         scans.per.second = "numeric", 
                                         avoid.processing.mz = "vector", 
                                         min.peak.width = "numeric", 
                                         min.peak.height = "numeric", 
                                         noise.threshold = "numeric", 
                                         compression.coef = "numeric")
)

# setClass(Class = "expClasses",
#          representation = representation(classes.type = "character",
#                                         classes.summary = "data.frame")
# )