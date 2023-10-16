#' @name computeRIerror
#' @aliases computeRIerror
#' @title computeRIerror
#' @description This function uses RI of mslib database and RT of the identified compounds to discrimine proper compound identification.
#' @param Experiment S4 object with experiment Data, Metadata and Results. Results of experiment are used to extract RT and Compound DB Id.
#' @param id.database Name of the preloaded database, in this case the regular db used by erah mslib
#' @param reference.list List with the compounds and their attributes (AlignId...)
#' @param ri.error.type Specify wether absolute or relative RI error is to be computed.
#' @param plot.results Shows the RI/RT graphic (True by default)
#' @details See eRah vignette for more details. To open the vignette, execute the following code in R:
#' vignette("eRahManual", package="erah")
#' @references [1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{showRTRICurve}} 
#' @examples \dontrun{
#' ex <- computeRIerror(
#'   ex, 
#'   mslib, 
#'   reference.list=list(AlignID = c(45,67,92,120)), 
#'   ri.error.type = "relative"
#' )
#' }
#' @export
#' @importFrom stats smooth.spline predict 
#' @importFrom graphics points

computeRIerror <- function(Experiment, id.database=mslib, reference.list, ri.error.type=c('relative','absolute'), plot.results=TRUE){
  
  if(!inherits(reference.list,'list')) stop('The parameter reference.list must be a list')
  if(is.null(reference.list)) stop('A reference list must be provided')
  if(!is.null(reference.list$AlignID) & !(is.null(reference.list$RT) | is.null(reference.list$RI))) stop('reference.list must contain a) user-defined RT and RI  values, or b) the AlignID of the compounds to be used as a reference')
  if(is.null(reference.list$AlignID) & (is.null(reference.list$RT) | is.null(reference.list$RI))) stop('reference.list must contain a) user-defined RT and RI  values, or b) the AlignID of the compounds to be used as a reference')
  if(!is.null(reference.list$RT)) if(length(reference.list$RT)!=length(reference.list$RI)) stop('Both RI and RT vectors in reference list must have the same length! (You provided a different number of RT and RI values, they must have the same length)')
  ri.error.type <- match.arg(ri.error.type, c('relative','absolute'))
  if(is.null(id.database)) stop("A database is needed for spectra comparison. Select a database or set 'compare' parameter to 'False'")
  colnames(Experiment@Results@Identification)
  if(nrow(Experiment@Results@Identification) == 0) stop("Factors must be identified first")
  if(Experiment@Results@Parameters@Identification$database.name != id.database@name) {
    error.msg <- paste("This experiment was not processed with the database selected. Please use ", 
                       Experiment@Results@Parameters@Identification$database.name, 
                       sep = "")
    stop(error.msg)
  }
  
  colRamp <- c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", 
               "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") 
  colSeq <- seq(100,0,length.out=8)
  
  n.putative <- Experiment@Results@Parameters@Identification$n.putative
  idlist <- idList(object=Experiment, id.database=id.database)
  
  genID <- idlist$AlignID
  rtVect.gen <- idlist$tmean	
  
  genInd <- which(Experiment@Results@Identification$AlignID %in% genID)
  riVect.gen <- sapply(as.numeric(as.vector(Experiment@Results@Identification$DB.Id.1[genInd])), function(i) id.database@database[[i]]$RI.VAR5.ALK)
  mf <- as.numeric(as.vector(idlist$MatchFactor.1))
  colVect <- sapply(mf, function(x) colRamp[which.min(abs(x - colSeq))])
  
  if(!is.null(reference.list$AlignID)){
    refInd <- which(Experiment@Results@Identification$AlignID %in% reference.list$AlignID)
    riVect.std <- sapply(as.numeric(as.vector(Experiment@Results@Identification$DB.Id.1[refInd])), function(i) id.database@database[[i]]$RI.VAR5.ALK)
    rtVect.std <- as.numeric(as.vector(Experiment@Results@Identification$tmean[refInd]))
    #plot(riVect.std, rtVect.std, type='n')
    #text(riVect.std, rtVect.std, text=1:length(riVect.std))
  }else{
    riVect.std <- reference.list$RI		
    rtVect.std <- reference.list$RT		
  }
  
  chrom.method <- list(name = 'RI/RT Calibration Curve', method = 'alk.var5', ref.rt = rtVect.std, ref.ri = riVect.std)
  RtIS <- smooth.spline(chrom.method$ref.rt , chrom.method$ref.ri, df=(length(chrom.method$ref.ri)-1))
  
  riMatColnames <- sapply(1:n.putative, function(x) paste('DB.Id.', x, sep=''))
  riMatNewColnames <- sapply(1:n.putative, function(x) paste('RI.error.', x, sep=''))
  refRImatrix <- apply(Experiment@Results@Identification[,riMatColnames], 2, function(x){
    sapply(as.numeric(as.vector(x)), function(i) id.database@database[[i]]$RI.VAR5.ALK)
  })
  
  rtVect <- Experiment@Results@Identification$tmean
  empRI <- predict(RtIS, rtVect)$y
  if(ri.error.type=='absolute') RI.error <- apply(refRImatrix, 2, function(x) abs(x-empRI))	
  if(ri.error.type=='relative') RI.error <- apply(refRImatrix, 2, function(x) round(100*abs(x-empRI)/empRI,2))		
  colnames(RI.error) <- riMatNewColnames
  
  Experiment@Results@Identification <- cbind(Experiment@Results@Identification, RI.error, stringsAsFactors=FALSE)
  
  if(plot.results){
    #colRamp <- RColorBrewer::brewer.pal(n = 11, name = "RdBu")
    #nColors <- 10
    #colRamp <- rainbow(n=nColors, s = 1, v = 1, start = 0, end = max(1, nColors - 1)/nColors, alpha = 0.85)
    #colSeq <- seq(100,0,length.out=8)
    #plot(1:length(colRamp), col=colRamp, pch=20, cex=5)		 
    colRamp <- c("red2", "#B2182B" ,"#D6604D", "#F4A582", "#FDDBC7", "#92C5DE")
    colSeq <- c(100,90,80,70,60,50)
    
    legend_label <- c(colSeq[-length(colSeq)], paste('<', colSeq[length(colSeq)], sep=''))
    
    par(mfrow=c(2,1))
    plot(rtVect.gen, riVect.gen, pch=20, col=colVect, main='Observed RI vs RT w/ Match Factors', xlab='RT (min)', ylab='RI (experimental)')
    legend('topleft', legend=legend_label,col=colRamp, y.intersp=0.8, cex=1, box.lwd=0.5, pt.cex=2, pch=15, title='MF value:', bg='white')	
    plot(rtVect.gen, riVect.gen, pch=20, col='grey', main='Observed RI vs RT w/ calibration curve', xlab='RT (min)', ylab='RI (experimental)')
    lines(RtIS, lwd=3, col='blue3', lty=2)
    points(RtIS, pch=8, col='red2')
  }
  
  Experiment
}

#' @name showRTRICurve
#' @aliases showRTRICurve
#' @title Show RT-RI curve
#' @description This function uses RI of mslib database and RT of the identified compounds to discrimine proper compound identification.
#' @param Experiment S4 object with experiment Data, Metadata and Results. Results of experiment are used to extract RT and Compound DB Id.
#' @param reference.list List with the compounds and their attributes (AlignId...)
#' @param nAnchors The desired equivalent number of degrees of freedom for the smooth.spline function
#' @param ri.thrs Retention Index treshold given by the user to discrimine bewteen identification results
#' @param id.database Name of the preloaded database (mslib by default, the regular db used by erah)
#' @details See eRah vignette for more details. To open the vignette, execute the following code in R:
#' vignette("eRahManual", package="erah")
#' @references [1] Xavier Domingo-Almenara, et al., eRah: A Computational Tool Integrating Spectral Deconvolution and Alignment with Quantification and Identification of Metabolites in GC-MS-Based Metabolomics. Analytical Chemistry (2016). DOI: 10.1021/acs.analchem.6b02927 
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{computeRIerror}} 
#' @examples \dontrun{
#' The following set erah to determine which indetified compounds are in RI treshold
#' RTRICurve <- showRTRICurve(ex, list, nAnchors=4, ri.thrs='1R')
#' }
#' @export

showRTRICurve <- function(Experiment, reference.list, nAnchors=4, ri.thrs='1R', id.database = mslib){
  if (!inherits(reference.list,"list")) 
    stop("The parameter reference.list must be a list")
  if (is.null(reference.list)) 
    stop("A reference list must be provided")
  if (!is.null(reference.list$AlignID) & !(is.null(reference.list$RT) | 
                                           is.null(reference.list$RI))) 
    stop("reference.list must contain a) user-defined RT and RI  values, or b) the AlignID of the compounds to be used as a reference")
  if (is.null(reference.list$AlignID) & (is.null(reference.list$RT) | 
                                         is.null(reference.list$RI))) 
    stop("reference.list must contain a) user-defined RT and RI  values, or b) the AlignID of the compounds to be used as a reference")
  if (!is.null(reference.list$RT)) 
    if (length(reference.list$RT) != length(reference.list$RI)) 
      stop("Both RI and RT vectors in reference list must have the same length! (You provided a different number of RT and RI values, they must have the same length)")
  if (is.null(id.database)) 
    stop("A database is needed for spectra comparison. Select a database or set 'compare' parameter to 'False'")
  if (nrow(Experiment@Results@Identification) == 0) 
    stop("Factors must be identified first")
  if (Experiment@Results@Parameters@Identification$database.name != 
      id.database@name) {
    error.msg <- paste("This experiment was not processed with the database selected. Please use ", 
                       Experiment@Results@Parameters@Identification$database.name, 
                       sep = "")
    stop(error.msg)
  }
  
  idlist <- idList(object = Experiment, id.database = id.database)
  
  if (!is.null(reference.list$AlignID)) {
    refInd <- which(Experiment@Results@Identification$AlignID %in% 
                      reference.list$AlignID)
    riVect.std <- sapply(as.numeric(as.vector(Experiment@Results@Identification$DB.Id.1[refInd])), 
                         function(i) id.database@database[[i]]$RI.VAR5.ALK)
    rtVect.std <- as.numeric(as.vector(Experiment@Results@Identification$tmean[refInd]))
  }
  else {
    riVect.std <- reference.list$RI
    rtVect.std <- reference.list$RT
  }   
  
  relativeThr <- NULL
  if(length(grep('r', tolower(ri.thrs)))==1) {
    relativeThr <- TRUE
    riThrs <- as.numeric(as.vector(gsub('r','', tolower(ri.thrs))))
  }	
  if(length(grep('a', tolower(ri.thrs)))==1) {
    relativeThr <- FALSE
    riThrs <- as.numeric(as.vector(gsub('a','', tolower(ri.thrs))))
  }	
  if(is.null(relativeThr)) stop('Incorrect "ri.thrs" value, please check the documentation for assinging "ri.thrs" values.')
  
  sSp <- smooth.spline(rtVect.std , riVect.std,df=nAnchors)
  riError <- abs(sSp$y - riVect.std)
  if(relativeThr) riError <- 100*(riError/sSp$y)
  
  plot(rtVect.std, riVect.std, type='n', main='RT/RI curve', xlab='RT (min)', ylab='RI')
  points(rtVect.std[which(riError<riThrs)], riVect.std[which(riError<riThrs)], pch=19)
  if(length(which(riError>riThrs))!=0) points(rtVect.std[which(riError>riThrs)], riVect.std[which(riError>riThrs)], pch=19, col='red2')
  lines(sSp)
  
  if(length(which(riError>riThrs))!=0){
    
    if (!is.null(reference.list$AlignID)) {
      cat('The points outside the RI threshold are: ')
      cat(paste0(Experiment@Results@Identification$AlignID[refInd][which(riError>riThrs)], collapse=', '))
    }else{
      cat('The points outside the RI threshold are: \n')
      #cbind(reference.list$RI[which(riError>riThrs)],reference.list$RT[which(riError>riThrs)])
    }
  }else{
    cat('No points outside the selected RI threshold')
  }	
}



