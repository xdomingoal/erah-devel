#' @importFrom stats smooth.spline predict
#' @importFrom graphics points

computeRIerror <- function(Experiment, id.database=mslib, reference.list, ri.error.type=c('relative','absolute'), plot.results=TRUE){
  
  if(class(reference.list)!='list') stop('The parameter reference.list must be a list')
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



