
globalVariables("mslib")

#' @useDynLib erah 

## Main Software functions:

#' @name setDecPar
#' @aliases setDecPar
#' @title Set Software Parameters
#' @description Sets Software Parameters for eRah.
#' @param min.peak.width Minimum compound peak width (in seconds). This is a critical parameter that conditions the efficiency of eRah. Typically, this should be the half of the mean compound width.
#' @param min.peak.height Minimum compound peak height
#' @param noise.threshold Data above this threshold will be considered as noise
#' @param avoid.processing.mz The masses that do not want to be considered for processing. Typically, in GC-MS those masses are 73,74,75,147,148 and 149, since they are they are ubiquitous mass fragments typically generated from compounds carrying a trimethylsilyl moiety.
#' @param compression.coef Data is compressed when using the orthogonal signal deconvolution (OSD) algorithm according to this value. A level 2 of compression is recomended.
#' @param analysis.time The chromatographic retention time window to process. If 0, all the chromatogram is processed.
#' @details See eRah vignette for more details. To open the vignette, execute the following code in R:
#' vignette("eRahManual", package="erah")
#' @references [1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{newExp}} \code{\link{deconvolveComp}} \code{\link{alignComp}} \code{\link{setAlPar}}
#' @examples \dontrun{
#' # The following will set eRah for analyzing the chromatograms
#' #from minutes 5 to 15, and withouth taking into account the masses
#' #35:69,73:75,147:149, widht a minimum peak widht of 0.7 seconds.
#' ex.dec.par <- setDecPar(min.peak.width = 0.7, 
#'                         min.peak.height = 5000, 
#'                         noise.threshold = 500, 
#'                         avoid.processing.mz = c(35:69,73:75,147:149), 
#'                         analysis.time = c(5,15))
#' }
#' @export

setDecPar <- function(min.peak.width, 
                      min.peak.height=2500, 
                      noise.threshold=500, 
                      avoid.processing.mz=c(73:75,147:149), 
                      compression.coef=2, 
                      analysis.time=0){
  softPar <- new("eRahSoftParameters",algorithm="eRah-OSD", min.peak.width = min.peak.width/60, min.peak.height = min.peak.height, noise.threshold = noise.threshold, avoid.processing.mz = avoid.processing.mz, compression.coef = compression.coef, analysis.time=analysis.time)
  softPar
}

#' @name setAlPar
#' @aliases setAlPar
#' @title Set Alignment Parameters
#' @description Setting alignment parameters for eRah.
#' @usage setAlPar(min.spectra.cor, max.time.dist,mz.range = c(70:600))
#' @param min.spectra.cor Minimum spectral correlation value. From 0 (non similar) to 1 (very similar). This value sets how similar two or more compounds have be to be considered for alignment between them.
#' @param max.time.dist Maximum retention time distance. This value (in seconds) sets how far two or more compounds can be to be considered for alignment between them.
#' @param mz.range The range of masses that is considered when comparing spectra.
#' @references [1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{newExp}} \code{\link{setDecPar}} \code{\link{alignComp}}
#' @examples\dontrun{
#' # The following will set eRah for aligning compounds which are
#' # at least 90 (per cent) similar, and which peaks are at a 
#' # maximum distance of 2 seconds. All the masses are considered when
#' # computing the spectral similarity.
#'
#' ex.al.par <- setAlPar(min.spectra.cor=0.90, max.time.dist=2,
#' mz.range=1:600)
#' }
#' @export

setAlPar <- function(min.spectra.cor, max.time.dist, mz.range=c(70:600))
{
  alPar <- new("eRahAlParameters", algorithm="eRah", min.spectra.cor=min.spectra.cor, max.time.dist=max.time.dist/60, mz.range = mz.range, method="eRah")
  alPar
}

#' @name newExp
#' @aliases newExp
#' @title New Experiment
#' @description Sets a new experiment for eRah
#' @param instrumental A data.frame containing the sample instrumental information.
#' @param phenotype (optional) A data.frame containing sample phenotype information.
#' @param info Experiment description
#' @details See eRah vignette for more details. To open the vignette, execute the following code in R:
#' vignette("eRahManual", package="erah")
#' @return \code{newExp} returns an S4 object of the class 'MetaboSet'.
#' @references [1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @examples \dontrun{
#' library(gcspikelite)
#' data(targets)
#' 
#' files <- list.files(system.file('data',package = 'gcspikelite'),full.names = TRUE)
#' files <- files[sapply(files,grepl,pattern = 'CDF')]
#' 
#' instrumental <- createInstrumentalTable(files)
#' phenotype <- createPhenoTable(files,as.character(targets$Group[order(targets$FileName)]))
#' 
#' ex <- newExp(instrumental = instrumental, 
#' phenotype = phenotype, info = "DEMO Experiment")
#' }
#' @seealso \code{\link{createInstrumentalTable}} \code{\link{createPhenoTable}} \code{\link{setDecPar}} \code{\link{setAlPar}}
#' @export

newExp<- function (instrumental, phenotype = NULL, info = character()) 
{ 
  if (is.null(phenotype)) {
    phenotype = as.data.frame(NULL)
    factors.list <- lapply(1:nrow(instrumental), function(x) {
      as.data.frame(NULL)
    })
    names(factors.list) <- as.vector(instrumental$sampleID)
    ident.list <- as.data.frame(matrix(0, ncol = 7, dimnames = list(row = 0, 
                                                                    col = c("AlignID", "tmean", "Name", "MatchFactor", "CAS", 
                                                                            "Formula", "DB.Id"))))
    uni.stats <- as.data.frame(matrix(0, ncol = 3, dimnames = list(row = 0, 
                                                                   col = c("Id", "FoldChangue", "pvalue"))))
    multi.stats <- as.data.frame(matrix(0, ncol = 3, dimnames = list(row = 0, 
                                                                     col = c("Id", "CompoundsInvolved", "pvalue"))))
    al.par <- list()
    id.par <- list()
    soft.par <- list()
    stat.parameters <- new("MSResultsParameters", Alignment = al.par, 
                           Identification = id.par)
    statistics <- new("Statistics", Univariate = uni.stats, 
                      Multivariate = multi.stats)
    MS.Results <- new("Results", Parameters = stat.parameters, 
                      Identification = ident.list, Statistics = statistics)
    MS.Data <- new("Data", FeatureList = list(NULL), FactorList = factors.list, 
                   Parameters = list(NULL))
    MS.MetaData <- new("MetaData", Instrumental = instrumental, 
                       Phenotype = phenotype, DataDirectory = "")
    col.correct <- c("sampleID", "filename", "date", "time")
    for (i in 1:length(col.correct)) if (length(apply(as.matrix(colnames(MS.MetaData@Instrumental)), 
                                                      1, function(x) grep(col.correct[i], x))) == 0) 
      stop("Invalid instrumental file. The file must contain at least the following columns: ", 
           paste(col.correct, collapse = ", "))
    if (!is.null(MS.MetaData@Phenotype)) 
      col.correct <- c("sampleID", "class")
    for (i in 1:length(col.correct)) if (length(apply(as.matrix(colnames(MS.MetaData@Phenotype)), 
                                                      1, function(x) grep(col.correct[i], x))) == 0) 
      sample.container <- new("MetaboSet", Info = info, Data = MS.Data, 
                              MetaData = MS.MetaData, Results = MS.Results)
    sample.container
    
  } else if (!is.null(phenotype)){

    factors.list <- lapply(1:nrow(instrumental), function(x) {
      as.data.frame(NULL)
    })
    names(factors.list) <- as.vector(instrumental$sampleID)
    ident.list <- as.data.frame(matrix(0, ncol = 7, dimnames = list(row = 0, 
                                                                    col = c("AlignID", "tmean", "Name", "MatchFactor", "CAS", 
                                                                            "Formula", "DB.Id"))))
    uni.stats <- as.data.frame(matrix(0, ncol = 3, dimnames = list(row = 0, 
                                                                   col = c("Id", "FoldChangue", "pvalue"))))
    multi.stats <- as.data.frame(matrix(0, ncol = 3, dimnames = list(row = 0, 
                                                                     col = c("Id", "CompoundsInvolved", "pvalue"))))
    al.par <- list()
    id.par <- list()
    soft.par <- list()
    stat.parameters <- new("MSResultsParameters", Alignment = al.par, 
                           Identification = id.par)
    statistics <- new("Statistics", Univariate = uni.stats, 
                      Multivariate = multi.stats)
    MS.Results <- new("Results", Parameters = stat.parameters, 
                      Identification = ident.list, Statistics = statistics)
    MS.Data <- new("Data", FeatureList = list(NULL), FactorList = factors.list, 
                   Parameters = list(NULL))
    MS.MetaData <- new("MetaData", Instrumental = instrumental, 
                       Phenotype = phenotype, DataDirectory = "")
    col.correct <- c("sampleID", "filename", "date", "time")
    for (i in 1:length(col.correct)) if (length(apply(as.matrix(colnames(MS.MetaData@Instrumental)), 
                                                      1, function(x) grep(col.correct[i], x))) == 0) 
      stop("Invalid instrumental file. The file must contain at least the following columns: ", 
           paste(col.correct, collapse = ", "))
    if (!is.null(MS.MetaData@Phenotype)) {
      col.correct <- c("sampleID", "class")
      for (i in 1:length(col.correct)) if (length(apply(as.matrix(colnames(MS.MetaData@Phenotype)), 
                                                        1, function(x) grep(col.correct[i], x))) == 0) 
        stop("Invalid phenotype file. The file must contain at least the following columns: ", 
             paste(col.correct, collapse = ", "))
    }
    sample.container <- new("MetaboSet", Info = info, Data = MS.Data, 
                            MetaData = MS.MetaData, Results = MS.Results)
    sample.container
  }
  
}


#' @rdname deconvolveComp
#' @title Deconvolution of compounds in samples
#' @description Deconvolution of GC-MS data
#' @param Experiment A 'MetaboSet' S4 object containing the experiment data previously created by newExp.
#' @param decParameters The software deconvolution parameters object previously created by setDecPar
#' @param samples.to.process Vector indicating which samples are to be processed.
#' @param down.sample If TRUE, chromatograms are down sampled to define one peak with 10 scan points (according to the minimum peak width). This is to process longer chromatograms with wider peak widths (more than 20 seconds peak width and small scans per second values). See details.
#' @param virtualScansPerSecond A virtual scans per second. If chromatograms are downsampled (for example, for a 1 mean peak width a 1 scans per second sampling frequency was used), eRah could not perform as expected. In these cases, the BEST solution is to re-acquire the samples. However, by selecting a different (virtual) scans per second frequency, eRah can upsample the data and process it more effectively.
#' @details See eRah vignette for more details. To open the vignette, execute the following code in R:
#' vignette("eRahManual", package="erah")
#'
#' eRah uses multivariate methods which run-time performance depend on the amount of data to be analyzed. When peaks are wider and the #' scans per second is also a small value, the number of points (scans) that define a peak might be too many, leading eRah to a poor run#'-time performance. To solve that, use down.sample=TRUE to allow eRah to define a peak with 10 seconds, and analyze the data more #' efficiently.
#' @return The function returns an updated S4 'MetaboSet' class, where the GC-MS samples have been now deconvolved.
#' @references [1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{newExp}} \code{\link{setAlPar}}
#' @examples \dontrun{
#' # Deconvolve data from a created experiment by \code{\link{newExp}}.
#' # ex <- newExp(instrumental="path")
#'
#' # The following will set eRah for analyzing the chromatograms
#' # from minutes 5 to 15, and withouth taking into account the masses
#' # 35:69,73:75,147:149, with a minimum peak width of 0.7 seconds.
#'
#' ex.dec.par <- setDecPar(min.peak.width=0.7, min.peak.height=5000, 
#'                        noise.threshold=500, avoid.processing.mz=c(35:69,73:75,147:149), 
#'                        analysis.time=c(5,15))
#'
#' # An now deconvolve the compounds in the samples:
#' # ex <- deconvolveComp(ex, decParameters=ex.dec.par)
#' }
#' @export

setGeneric('deconvolveComp',function(Experiment, decParameters, samples.to.process=NULL, down.sample=FALSE, virtualScansPerSecond=NULL){
  standardGeneric('deconvolveComp')
})

#' @rdname deconvolveComp
#' @importFrom furrr future_map furrr_options

setMethod('deconvolveComp',signature = 'MetaboSet',
          function(Experiment, decParameters, samples.to.process=NULL, down.sample=FALSE, virtualScansPerSecond=NULL){
            plotting=FALSE
            Number.of.Samples <- nrow(Experiment@MetaData@Instrumental)
            if(is.null(samples.to.process)) samples.to.process <- 1:Number.of.Samples
            stopifnot(samples.to.process>=1, max(samples.to.process)<=Number.of.Samples, length(samples.to.process)<=Number.of.Samples)
            
            soft.par <- list(min.peak.width = decParameters@min.peak.width, min.peak.height = decParameters@min.peak.height, noise.threshold = decParameters@noise.threshold, avoid.processing.mz = decParameters@avoid.processing.mz,  compression.coef = decParameters@compression.coef, analysis.time = decParameters@analysis.time)
            Experiment@Data@Parameters <- soft.par
            
            Experiment@Data@FactorList <- future_map(samples.to.process,~{
              
              k <- which(samples.to.process == .x)
              
              message(paste("\n Deconvolving compounds from",
                  as.character(Experiment@MetaData@Instrumental$filename[k]),
                  "... Processing", 
                  k,"/",
                  length(samples.to.process),"\n"))
              processSample(Experiment, 
                            .x, 
                            plotting, 
                            down.sample, 
                            virtualScansPerSecond)
            },
            .options = furrr_options(seed = TRUE))
            
            names(Experiment@Data@FactorList) <- metaData(Experiment)$sampleID
            Experiment <- scansPerSecond(Experiment)
            message("Compounds deconvolved")
            Experiment	
          }
)

#' @rdname alignComp
#' @title Alignment of compounds
#' @description Alignment of GC-MS deconvolved compounds
#' @usage alignComp(Experiment, alParameters, blocks.size=NULL)
#' @param Experiment A 'MetaboSet' S4 object containing the experiment data previously created by newExp and deconvolved by deconvolveComp.
#' @param alParameters The software alignment parameters object previously created by setAlPar
#' @param blocks.size For experiment of more than 1000 samples, and depending on the computer, alignment can be conducted by block segmentation. See details.
#' @details See eRah vignette for more details. To open the vignette, execute the following code in R:
#' vignette("eRahManual", package="erah")
#'
#' For experiments containing more than 100 (Windows) or 1000 (Mac or Linux) samples (numbers depending on the computer resoures and sample type). In those cases alignment can be conducted by block segmentation. For an experiment of e.g. 1000 samples, the block.size can be set to 100, so the alignment will perform as multiple (ten) 100-samples experiments, to later align them into a single experiment.
#'
#' This parameter is designed to solve the typical problem that appear when aligning under Windows operating system: "Error: cannot allocate vector of size XX Gb". Such a problem will not appear with Mac or Linux, but several hours of computation are expected when aligning a large number of samples. Using block segmentation provides a greatly improved run-time performance.
#' @return The function returns an updated S4 'MetaboSet' class, where the GC-MS samples have been now aligned.
#' @references [1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{newExp}} \code{\link{setDecPar}} \code{\link{deconvolveComp}}
#' @export

setGeneric('alignComp',function(Experiment, alParameters, blocks.size=NULL){
  standardGeneric('alignComp')
})

#' @rdname alignComp

setMethod('alignComp',signature = 'MetaboSet',
          function(Experiment, alParameters, blocks.size=NULL){
            al.par <- list(alignment.algorithm=alParameters@algorithm, min.spectra.cor=alParameters@min.spectra.cor, max.time.dist=alParameters@max.time.dist, mz.range=alParameters@mz.range)
            Experiment@Results@Parameters@Alignment <- al.par
            
            min.spectra.cor <- Experiment@Results@Parameters@Alignment$min.spectra.cor
            max.time.dist <- Experiment@Results@Parameters@Alignment$max.time.dist
            mz.range <- Experiment@Results@Parameters@Alignment$mz.range
            maxMZ <- max(mz.range)
            
            # Experiment@Data@FactorList <- align.factors(Experiment@Data@FactorList, min.spectra.cor, max.time.dist, maxMZ, mz.range)
            # Experiment@Results@Alignment <- create.factorlist.table(Experiment)
            
            if(is.null(blocks.size))
            {
              Experiment@Data@FactorList <- align.factors(Experiment@Data@FactorList, min.spectra.cor, max.time.dist, maxMZ, mz.range)
              Experiment@Results@Alignment <- create.factorlist.table(Experiment)
              #return(Experiment)
            }else{
              
              #blocks.size <- 15
              max.mz <- maxMZ
              Itrt <- length(Experiment@Data@FactorList)/blocks.size
              sequs <- trunc(seq(1, length(Experiment@Data@FactorList), length.out=Itrt))
              sequs[1] <- 0
              
              corresponding.list <- list()
              block.list <- list()
              #i <- 1
              
              for(i in 1:(length(sequs)-1))	
              {
                cat("Aligning block ", i, " of ", length(sequs)-1, "... \n", sep="")
                ghost.object <- Experiment
                ghost.object@Data@FactorList <- Experiment@Data@FactorList[(sequs[i]+1):sequs[(i+1)]]
                factors.list <- ghost.object@Data@FactorList
                ghost.object@Data@FactorList <- align.factors(factors.list, min.spectra.cor, max.time.dist, max.mz, mz.range)
                ghost.factors.list <- create.factorlist.table(ghost.object)
                
                block.list[[i]] <- data.frame(ID=ghost.factors.list$AlignID, RT=ghost.factors.list$tmean, Spectra=ghost.factors.list$Spectra)
                corresponding.list <- c(corresponding.list,lapply(ghost.object@Data@FactorList, function(x) x$AlignID))		
              }
              
              cat("Aligning factors across blocks... \n")
              full.factorlist <- align.factors(block.list, min.spectra.cor, max.time.dist, max.mz, mz.range)
              
              #MaxALID <- max(unlist(lapply(full.factorlist, function(x) x$AlignID)))
              factors.list <- Experiment@Data@FactorList
              if(!(any(unlist(lapply(factors.list,function(x) {is.null(x$AlignID)}))==FALSE)))
              {	
                factors.list <- lapply(factors.list, function(x){
                  outp <- cbind(x,matrix(0,nrow=length(x$ID)))
                  colnames(outp)[ncol(outp)] <- "AlignID"
                  outp
                })
              }else{
                factors.list <- lapply(factors.list, function(x){
                  x$AlignID <- rep(0,length(x$ID))
                  x
                })
              }	
              
              Experiment@Data@FactorList <- factors.list
              
              free.aligned.slots <- list()		
              for(i in 1:length(full.factorlist))
              {
                for(j in (sequs[i]+1):sequs[(i+1)])
                {
                  ID.vct <- sapply(full.factorlist[[i]]$ID, function(x) {x.num <- which(corresponding.list[[j]]==x)
                  if(length(x.num)==0) x.num=0
                  x.num
                  })
                  
                  #full.factorlist[[i]]$AlignID[which(ID.vct!=0)]	
                  #ID.vct[which(ID.vct!=0)]	
                  
                  Experiment@Data@FactorList[[j]]$AlignID[ID.vct[which(ID.vct!=0)]] <- full.factorlist[[i]]$AlignID[which(ID.vct!=0)]	
                  free.aligned.slots[[j]] <- which(full.factorlist[[i]]$AlignID[which(ID.vct!=0)]==0)
                }	
              }
              MaxALID <- max(unlist(lapply(Experiment@Data@FactorList, function(x) x$AlignID)))
              Alid.counter <- MaxALID + 1
              
              for(i in 1:length(free.aligned.slots))
              {
                Experiment@Data@FactorList[[i]]$AlignID[free.aligned.slots[[i]]] <- seq(Alid.counter, Alid.counter + (length(free.aligned.slots[[i]])-1) )
                Alid.counter <- Alid.counter + length(free.aligned.slots[[i]])
              }
              
              cat("Constructing Factor List Table... (This may take a while...)\n")	
              Experiment@Results@Alignment <- create.factorlist.table(Experiment)
              
            }
            
            Experiment
          }
)

#' @rdname identifyComp
#' @title Identification of compounds
#' @description Identification of compounds. Each empirical spectrum is compared against a ms library.
#' @usage identifyComp(Experiment, id.database = mslib,mz.range = NULL, n.putative = 3)
#' @param Experiment A 'MetaboSet' S4 object containing the experiment data previously created by newExp, deconvolved by deconvolveComp and optionally aligned by alignComp.
#' @param id.database The mass-spectra library to be compared with the empirical spectra. By default, the MassBank-[2] - Mass Bank of North America (MoNa) database are employed.
#' @param mz.range The same as in alignComp. If specified already in alignComp, then there is no need to especify it again. If not, it has to be specified.
#' @param n.putative The number of hits (compound candidate names) to be returned for each spectrum found.
#' @return The function returns an updated S4 'MetaboSet' class, where the GC-MS samples have been now aligned.
#' @references 
#' [1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
#'
#' [2] MassBank: A public repository for sharing mass spectral data for life sciences, H. Horai, M. Arita, S. Kanaya, Y. Nihei, T. Ikeda, K. Suwa. Y. Ojima, K. Tanaka, S. Tanaka, K. Aoshima, Y. Oda, Y. Kakazu, M. Kusano, T. Tohge, F. Matsuda, Y. Sawada, M. Yokota Hirai, H. Nakanishi, K. Ikeda, N. Akimoto, T. Maoka, H. Takahashi, T. Ara, N. Sakurai, H. Suzuki, D. Shibata, S. Neumann, T. Iida, K. Tanaka, K. Funatsu, F. Matsuura, T. Soga, R. Taguchi, K. Saito and T. Nishioka, J. Mass Spectrom., 45 (2010) 703-714. 
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{newExp}} \code{\link{alignComp}} \code{\link{setAlPar}} \code{\link{setDecPar}}
#' @export

setGeneric('identifyComp',function(Experiment, id.database=mslib, mz.range=NULL, n.putative=3){
  standardGeneric('identifyComp')
})

#' @rdname identifyComp

setMethod('identifyComp',signature = 'MetaboSet',
          function(Experiment, id.database=mslib, mz.range=NULL, n.putative=3){
            #if(!(any(unlist(lapply(Experiment@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")
            
            if(is.null(Experiment@Results@Parameters@Alignment$mz.range) && is.null(mz.range)) stop("A mz.range has to be specified")
            if(is.null(mz.range)) compare.only.mz <- min(Experiment@Results@Parameters@Alignment$mz.range):max(Experiment@Results@Parameters@Alignment$mz.range)
            if(!is.null(mz.range)) compare.only.mz <- mz.range
            
            
            id.par <- list(database.name = id.database@name, compare.only.mz = compare.only.mz, n.putative = n.putative)
            Experiment@Results@Parameters@Identification <- id.par
            
            avoid.processing.mz <- Experiment@Data@Parameters$avoid.processing.mz
            maxMZ <- max(compare.only.mz)
            Experiment@Results@Identification <- identify.factors(Experiment, maxMZ, compare.only.mz, avoid.processing.mz, id.database@database, n.putative)
            Experiment
          }
)

processSample <- function(Experiment, index, plotting, down.sample, virtual.scans.ps)
{
  if(Experiment@MetaData@DataDirectory=="") {filename <- as.character(Experiment@MetaData@Instrumental$filename[index])
  }else{filename <- paste(Experiment@MetaData@DataDirectory,"/",Experiment@MetaData@Instrumental$filename[index], sep="")}
  
  sampleObject <- NULL
  sampleObject <- load.file(filename)
  
  # file.extension <- strsplit(as.character(Experiment@MetaData@Instrumental$filename[index]), split="\\.")[[1]]
  # file.type <- file.extension[length(file.extension)]	
  # if(file.type=="cdf") sampleObject <- load.ncdf(filename)
  # if(file.type=="mzXML" || file.type=="xml") sampleObject <- load.xml(filename)
  # if(file.type=="MetaboSet")
  # {
  # load(filename)
  # sampleObject <- new("RawDataParameters", data = sampleRD@data, min.mz = sampleRD@min.mz, max.mz = sampleRD@max.mz, start.time = sampleRD@start.time, mz.resolution = 1)
  # }
  
  Experiment@Data@Parameters$scans.per.second <- sampleObject@scans.per.second
  sampleObject@avoid.processing.mz <- Experiment@Data@Parameters$avoid.processing.mz
  sampleObject@min.peak.width <- Experiment@Data@Parameters$min.peak.width*Experiment@Data@Parameters$scans.per.second*60
  sampleObject@min.peak.height <- Experiment@Data@Parameters$min.peak.height
  sampleObject@noise.threshold <- Experiment@Data@Parameters$noise.threshold
  #sampleObject@moving.window.length <- Experiment@Data@Parameters$moving.window.length*Experiment@Data@Parameters$scans.per.second*60
  #sampleObject@moving.window.overlap <- Experiment@Data@Parameters$moving.window.overlap
  sampleObject@compression.coef <- Experiment@Data@Parameters$compression.coef
  #sampleObject@factor.minimum.sd <- Experiment@Data@Parameters$factor.minimum.sd
  
  #sampleObject@filter.matrix <- get.filter.matrix(sampleObject)
  
  
  sampleObject <- avoid.processing(sampleObject)
  factor.list <- try(get.factor.list(sampleObject, analysis.window=Experiment@Data@Parameters$analysis.time, plotting, down.sample, virtual.scans.ps), silent=F)
  if(class(factor.list)=="try-error") {factor.list <- as.data.frame(NULL); warning("Unable to extract factors from ", Experiment@MetaData@Instrumental$filename[index], ". Data may be corrupted.", sep="")}
  
  factor.list		
}

scansPerSecond <- function(Experiment){
  if (Experiment@MetaData@DataDirectory=="") {
    filename <- as.character(Experiment@MetaData@Instrumental$filename[1])
  } else {
    filename <- paste(Experiment@MetaData@DataDirectory,"/",Experiment@MetaData@Instrumental$filename[1], sep="")
  }
  sampleObject <- load.file(filename)
  Experiment@Data@Parameters$scans.per.second <- sampleObject@scans.per.second
  return(Experiment)
}
