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

setClass(Class = "MetaboSet",
         representation = representation(Info = "character", 
                                         Data = "Data", 
                                         MetaData = "MetaData", 
                                         Results = "Results")
         )

## Intern Classes for eRah:

setClass(Class = "eRah_DB", 
         representation = representation(name = "character", 
                                         version = "character", 
                                         info = "character", 
                                         database = "list")
         )

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
