#' @name export2MSP
#' @aliases export2MSP
#' @title Export spectra to MSP
#' @description Export spectra to MSP format for comparison with the NIST library.
#' @usage export2MSP(Experiment, export.id = NULL, 
#' id.database = mslib, store.path = getwd())
#' @param Experiment A 'MetaboSet' S4 object containing the experiment.
#' @param export.id If NULL, all the spectra in the experiment will be exported. Otherwise, only the AlignID in export.id will be exported
#' @param id.database The mass-spectra library used in the experiment.
#' @param store.path The path where the converted files are to be exported.

export2MSP <- function(Experiment, export.id=NULL, id.database = mslib, store.path=getwd())
{
	#Experiment <- ex
	#export.id <- NULL
	#export.id <- c(4,5,10)
	
	#if(dir.exists(paste(c(store.path, "/ExportMSP/"), collapse=""))) stop("Please, delete the following folder before proceeding: ", paste(c(store.path, "/ExportMSP/")))

	
	if(is.null(nrow(Experiment@Results@Identification)) | nrow(Experiment@Results@Identification)==1)
	{
		if(!is.null(export.id)) Experiment@Results@Alignment <- Experiment@Results@Alignment[which(Experiment@Results@Alignment$AlignID %in% export.id),]
			
		SpectList <- sapply(Experiment@Results@Alignment$Spectra, function(x) {
			splitted.spectra.list <- strsplit(as.character(x), split = " ")[[1]]
			splitted.spectra.list <- gsub(",", " ", splitted.spectra.list)
			splitted.spectra.list <- as.character(as.vector(sapply(splitted.spectra.list, function(x) paste(c(x,";"), collapse=""))))
			Npeaks <- length(splitted.spectra.list)
			if(Npeaks>=6)
			{
				sequ <- seq(1, Npeaks, 5)
				if(sequ[length(sequ)]!=Npeaks) sequ <- c(sequ, Npeaks)
				splitted.spectra.list <- paste(unlist(sapply(1:(length(sequ)-1), function(x) paste(c(splitted.spectra.list[sequ[x]:sequ[(x+1)]], "\n")))), collapse=" ")
			}else{
				splitted.spectra.list <- paste(splitted.spectra.list, collapse=" ")
				}
			PeakChar <- paste(gsub(",", " ", splitted.spectra.list), collapse="; ")
			return(list(Npeaks=Npeaks, PeakChar=PeakChar))
		})
		SpectList <- split(SpectList,seq(NROW(SpectList)))

		SpectString <- as.vector(unlist(SpectList[[2]]))
		Npeaks <- as.vector(unlist(SpectList[[1]]))
		
		SpectNames.2 <- as.character(as.vector(Experiment@Results@Alignment$Factor))
		SpectNames.3 <- sapply(as.numeric(as.vector(Experiment@Results@Alignment$tmean)), function(x) paste("Rt:", x))
		
		SpectNames <- apply(cbind(SpectNames.3,SpectNames.2), 1, function(x) paste(x, collapse=" @ "))
	}else{
		
		if(!is.null(export.id)) Experiment@Results@Identification<- Experiment@Results@Identification[which(Experiment@Results@Identification$AlignID %in% export.id),]
		
		SpectList <- sapply(Experiment@Results@Identification$Spectra, function(x) {
			splitted.spectra.list <- strsplit(as.character(x), split = " ")[[1]]
			splitted.spectra.list <- gsub(",", " ", splitted.spectra.list)
			splitted.spectra.list <- as.character(as.vector(sapply(splitted.spectra.list, function(x) paste(c(x,";"), collapse=""))))
			Npeaks <- length(splitted.spectra.list)
			if(Npeaks>=6)
			{
				sequ <- seq(1, Npeaks, 5)
				if(sequ[length(sequ)]!=Npeaks) sequ <- c(sequ, Npeaks)
				splitted.spectra.list <- paste(unlist(sapply(1:(length(sequ)-1), function(x) paste(c(splitted.spectra.list[sequ[x]:sequ[(x+1)]], "\n")))), collapse=" ")
			}else{
				splitted.spectra.list <- paste(splitted.spectra.list, collapse=" ")
				}
			PeakChar <- paste(gsub(",", " ", splitted.spectra.list), collapse="; ")
			return(list(Npeaks=Npeaks, PeakChar=PeakChar))
		})
		SpectList <- split(SpectList,seq(NROW(SpectList)))

		SpectString <- as.vector(unlist(SpectList[[2]]))
		Npeaks <- as.vector(unlist(SpectList[[1]]))
		
        id.found <- as.numeric(as.vector(Experiment@Results@Identification[, "DB.Id.1"]))
		met.name <- unlist(lapply(id.database@database[c(id.found)], function(x) x$Name))
		
		SpectNames.1 <- sapply(as.character(as.vector(Experiment@Results@Identification$AlignID)), function(x) paste("(Align ID: ", x, ")", sep=""))
		SpectNames.2 <- sapply(as.numeric(as.vector(Experiment@Results@Identification$tmean)), function(x) paste("Rt:", x))
		SpectNames.3 <- apply(cbind(SpectNames.2,SpectNames.1), 1, function(x) paste(x, collapse=" "))
		SpectNames <- as.character(as.vector(apply(cbind(SpectNames.3,met.name), 1, function(x) paste(x, collapse=" @ "))))
	}
	
	dir.create(file.path(store.path, "ExportMSP") , showWarnings = FALSE)
	filename <- paste(c(store.path, "/ExportMSP/", "ExportedMSP", ".msp"), collapse="")
		
	fileTag <- character()	
	for(i in 1:length(SpectString))
	{ 
		fileTag[i] <- paste(c("Name: ", SpectNames[i], "\n", "Comments: MSP spectra exported by eRah \n", "Num Peaks: ", Npeaks[i], "\n",SpectString[i]), collapse="")
	}
	fileTagGen <- paste(fileTag, collapse="\n \n")

	writeLines(fileTagGen, filename)	
	cat("Spectra saved at: ", store.path, "/ExportMSP", sep="")
}

#' @name export2CEF
#' @aliases export2CEF
#' @title Export spectra to CEF
#' @description Export spectra to CEF format for comparison with the NIST library through MassHunter interface.
#' @usage export2CEF(Experiment, export.id = NULL, 
#' id.database = mslib, store.path = getwd())
#' @param Experiment A 'MetaboSet' S4 object containing the experiment.
#' @param export.id If NULL, all the spectra in the experiment will be exported. Otherwise, only the AlignID in export.id will be exported
#' @param id.database The mass-spectra library used in the experiment.
#' @param store.path The path where the converted files are to be exported.

export2CEF <- function (Experiment, export.id = NULL, id.database = mslib, store.path = getwd()) 
{
    if (!is.null(export.id)) Experiment@Results@Alignment <- Experiment@Results@Alignment[which(Experiment@Results@Alignment$AlignID %in% export.id), ]
        SpectList <- sapply(Experiment@Results@Alignment$Spectra, function(x) {
                splitted.spectra.list <- strsplit(as.character(x), split = " ")[[1]]
  
                as.vector(sapply(splitted.spectra.list, function(y){
                	y.split <- strsplit(y, ",")[[1]]
                	x.cor <- y.split[1]
                	y.cor <- y.split[2]
                	paste("<p x=", '"', as.numeric(x.cor), '" y="', y.cor, '" />', sep="")
                }))
                })

        SpectNames.2 <- as.character(as.vector(Experiment@Results@Alignment$Factor))
        SpectRT <- as.numeric(as.vector(Experiment@Results@Alignment$tmean))
        SpectNames.3 <- sapply(as.numeric(as.vector(Experiment@Results@Alignment$tmean)), function(x) paste("Rt:", x))
        SpectNames <- apply(cbind(SpectNames.2, SpectNames.3), 1, function(x) paste(x, collapse = " @ "))
    
    dir.create(file.path(store.path, "ExportCEF"), showWarnings = FALSE)
    filename <- paste(c(store.path, "/ExportCEF/", "ExportedCEF", ".cef"), collapse = "")
    fileTag <- character()
    fileTag[1] <- paste("<?xml version=", '"', "1.0", '"' , "encoding=", '"', "utf-8",'"', "?><CEF version=",'"',"1.0.0.0",'"',"><CompoundList>", sep="")
    for (i in 1:length(SpectList)) {
    	fileTag[length(fileTag) + 1] <- ""
    	fileTag[length(fileTag) + 1] <- paste("<Compound mppid=", '"',SpectNames[i], '"'," algo=", '"',"FindByAMDIS", '"'," >", sep="")
    	fileTag[length(fileTag) + 1] <- paste("<Location  rt=", '"',SpectRT[i], '"'," />", sep="")
    	fileTag[length(fileTag) + 1] <- ""

		fileTag[length(fileTag) + 1] <- paste("<Spectrum type=", '"',"AMDIS", '"'," cpdAlgo=", '"',"FindByAMDIS", '"',">", sep="")
		fileTag[length(fileTag) + 1] <- paste("<MSDetails  p=", '"',"+", '"',"  />", sep="")
		fileTag[length(fileTag) + 1] <- "<MSPeaks>" 
    	srtI <- (length(fileTag) + 1)
    		srtE <- (length(fileTag) + 1) + length(SpectList[[i]]) - 1
    	fileTag[srtI:srtE] <- SpectList[[i]]
    		fileTag[length(fileTag) + 1] <- "</MSPeaks>" 
		fileTag[length(fileTag) + 1] <- "</Spectrum>" 
		fileTag[length(fileTag) + 1] <- "</Compound>" 

    }
    fileTag[length(fileTag) + 1] <- ""
    fileTag[length(fileTag) + 1] <- "</CompoundList></CEF>"

    #fileTagGen <- paste(fileTag, collapse = "\n \n")
    writeLines(fileTag, filename)
    cat("Spectra saved at: ", store.path, "/ExportCEF", sep = "")
}








