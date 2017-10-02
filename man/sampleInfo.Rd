\name{sampleInfo}
\alias{sampleInfo}
\title{Information of the samples}
\description{Returns basic information on the samples.}
\usage{
sampleInfo(Experiment, N.sample=1)
}
\arguments{
  \item{Experiment}{A 'MetaboSet' S4 object containing the experiment.}
    \item{N.sample}{Integer. The number of the sample to query.}
}
\details{
Returns details on a given sample of the experiment, such as name, start time, end time, minium and maximum adquired m/z and scans per second.
}

\seealso{
\code{\link{plotChr}}

}