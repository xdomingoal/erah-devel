\name{RawDataParameters-class}
\docType{class}
\alias{RawDataParameters-class}
\title{Class \code{"RawDataParameters"}}
\description{
The RawDataParameters class contains the slots for storing and accessing into a MS sample, and the essential parameters for performing its processing (deconvolution).
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{The data matrix of the sample to be processed}
    \item{\code{min.mz}:}{The minimum adquired mz number}
	\item{\code{max.mz}:}{The maximum adquired mz number}
    \item{\code{start.time}:}{Starting time of adquisition}
    \item{\code{mz.resolution}:}{Mz resolution}
    \item{\code{scans.per.second}:}{Scans per second}
    \item{\code{avoid.processing.mz}:}{Which mz do not have to be processed}
    \item{\code{min.peak.width}:}{Minimum peak width (stored in scans)}
    \item{\code{min.peak.height}:}{Minimum peak height}
    \item{\code{noise.threshold}:}{The noise threshold}
    \item{\code{compression.coef}:}{Compression coefficient (parameter for Orthogonal Signal Deconvolution)}
  }
}

\author{
Xavier Domingo-Almenara.
}