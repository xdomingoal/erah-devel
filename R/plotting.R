#' @name plotSpectra
#' @aliases plotSpectra
#' @title Plotting Spectra
#' @description Plots the empirical spectra found by eRah, and allows comparing it with the reference spectra.
#' @usage 
#' plotSpectra(Experiment, AlignId, n.putative = 1,
#' compare = T, id.database = mslib, comp.db = NULL, 
#' return.spectra = F, draw.color = "purple", xlim = NULL)
#' @param Experiment A 'MetaboSet' S4 object containing the experiment after being deconolved, aligned and (optionally) identified.
#' @param AlignId the Id identificator for the compound to be shown.
#' @param n.putative The hit number (position) to be returned when comparing the empirical spectrum with the reference. See details
#' @param compare logical. If TRUE, then the reference spectrum from the library is shown for comparison.
#' @param id.database The mass-spectra library to be compared with the empirical spectra. By default, the MassBank-[2] - Mass Bank of North America (MoNa) database are employed.
#' @param comp.db If you want to compare the empirical spectrum with another spectrum from the database, select the comp.db number from the database.
#' @param return.spectra logical. If TRUE, the function returns the empirical spectrum for the selected compound
#' @param draw.color Selects the color for the reference spectrum (see \code{\link{colors}}).
#' @param xlim x axsis (mass - m/z) limits (see \code{\link{plot.default}}).
#' @details When identification is applied (see \code{\link{identifyComp}}), the number of hits to be returned (n.putative) has to be selected. Therefore, here you can compare the empirical spectrum (found by eRah) with each n.putative hit returned (1, 2, ...) by (see \code{\link{identifyComp}}).
#' @return 
#' \code{plotSpectra} returns an vector when return.spectra=TRUE.
#'
#'      \item{x}{vector. Containts the empirical spectrum.}
#' @references 
#' [1] eRah: an R package for spectral deconvolution, alignment, and metabolite identification in GC/MS-based untargeted metabolomics. Xavier Domingo-Almenara, Alexandre Perera, Maria Vinaixa, Sara Samino, Xavier Correig, Jesus Brezmes, Oscar Yanes. (2016) Article in Press.
#'
#' [2] MassBank: A public repository for sharing mass spectral data for life sciences, H. Horai, M. Arita, S. Kanaya, Y. Nihei, T. Ikeda, K. Suwa. Y. Ojima, K. Tanaka, S. Tanaka, K. Aoshima, Y. Oda, Y. Kakazu, M. Kusano, T. Tohge, F. Matsuda, Y. Sawada, M. Yokota Hirai, H. Nakanishi, K. Ikeda, N. Akimoto, T. Maoka, H. Takahashi, T. Ara, N. Sakurai, H. Suzuki, D. Shibata, S. Neumann, T. Iida, K. Tanaka, K. Funatsu, F. Matsuura, T. Soga, R. Taguchi, K. Saito and T. Nishioka, J. Mass Spectrom., 45, 703-714 (2010)
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{plotProfile}} \code{\link{plotAlign}}
#' @export
#' @importFrom graphics text lines

plotSpectra <- function(Experiment, AlignId, n.putative=1, compare=T, id.database=mslib, comp.db=NULL, return.spectra=F, draw.color="purple", xlim=NULL)
{
	if(length(AlignId)!=1) stop("Only one spectrum can be shown at once")
	if(compare==T) if(is.null(id.database)) stop("A database is needed for spectra comparison. Select a database or set 'compare' parameter to 'False'")

	index <- which(as.numeric(as.vector(Experiment@Results@Alignment[,"AlignID"]))==AlignId)
	MSP.spect.emp <- Experiment@Results@Alignment[index,"Spectra"]

	current.column <- paste("DB.Id.",n.putative, sep="")

	if(compare==T) 
	{
		if(is.null(comp.db)) if(nrow(Experiment@Results@Identification)==1) stop("Factors must be identified first")
		if(is.null(comp.db)) MSP.spect.db <- id.database@database[[as.numeric(as.vector(Experiment@Results@Identification[index,current.column]))]]$Spectra
		if(!is.null(comp.db)) MSP.spect.db <- id.database@database[[comp.db]]$Spectra
		
		maxMZ <- max(Experiment@Results@Parameters@Identification$compare.only.mz)
		empiric.name <- id.database@database[[as.numeric(as.vector(Experiment@Results@Identification[index,current.column]))]]$Name
		
		empiric.spectra <- convertMSPspectra(MSP.spect.emp,maxMZ)
	
		delete.mz <- 1:length(empiric.spectra)
		delete.mz <- delete.mz[-Experiment@Results@Parameters@Identification$compare.only.mz]
		if(length(delete.mz)!=0) empiric.spectra[delete.mz] <- 0
		
		empiric.spectra <- empiric.spectra/max(empiric.spectra)*1000
			
	}else{
		maxMZ <- max(Experiment@Results@Parameters@Alignment$mz.range)
		empiric.name <- paste("Factor #",AlignId, sep="")
		
		empiric.spectra <- convertMSPspectra(MSP.spect.emp,maxMZ)
		empiric.spectra <- empiric.spectra/max(empiric.spectra)*1000

	}
		
	mz.len <- 5
	if(length(which(empiric.spectra!=0))<5) mz.len <- length(which(empiric.spectra!=0))
	main_mz.empiric <- sort(empiric.spectra, decreasing=T, index.return=T)$ix[1:mz.len]

	if(compare==F)
	{	
		plot(empiric.spectra, type="h", main=paste(empiric.name, "\n (Empirical Spectra)"), xlab="Mz", ylab="Intensity",xlim=xlim)
		text(main_mz.empiric, empiric.spectra[main_mz.empiric], labels=main_mz.empiric, cex=0.7)	
	}else{
		
		db.spectra <- convertMSPspectra.dot(MSP.spect.db,maxMZ)
		delete.mz1 <- Experiment@Data@Parameters$avoid.processing
		delete.mz2 <- 1:length(empiric.spectra)
		delete.mz2 <- delete.mz2[-Experiment@Results@Parameters@Identification$compare.only.mz]
		delete.mz <- unique(c(delete.mz1,delete.mz2))
		if(length(delete.mz)!=0) db.spectra[delete.mz] <- 0
		db.spectra <- normalize(db.spectra)*(-1000)

		match.factor <- cor.sinus(empiric.spectra,abs(db.spectra))
		match.factor <- round(match.factor*100, digits=1)
		
		if(is.null(comp.db)) main.title <- paste(empiric.name, "\n Match Factor:",match.factor)
		if(!is.null(comp.db)) main.title <- paste("Align ID #", AlignId, " VS ", id.database@database[[comp.db]]$Name, "\n Match Factor:",match.factor, sep="")
		
		plot(empiric.spectra, type="h", main=main.title, ylim=c(-1000,1000), xlab="Mz", ylab="Intensity", xlim=xlim)
		text(main_mz.empiric, empiric.spectra[main_mz.empiric], labels=main_mz.empiric, cex=0.7)	
		lines(db.spectra, col=draw.color, type="h")
		main_mz.db <- sort(abs(db.spectra), decreasing=T, index.return=T)$ix[1:mz.len]
		text(main_mz.db, db.spectra[main_mz.db], labels=main_mz.db, cex=0.7, col=draw.color)	
		if(is.null(comp.db)) legend("topright", legend=c("Empirical","Database"), col=c("black",draw.color), pch=19)
		if(!is.null(comp.db)) legend("topright", legend=c(paste("Align ID #", AlignId, sep=""),id.database@database[[comp.db]]$Name), col=c("black",draw.color), pch=19)
	}
	if(return.spectra==T) return(empiric.spectra)
	
}

#' @name plotProfile
#' @aliases plotProfile
#' @title Plotting chromatographic profile
#' @description Plots the chromatophic profiles of the compounds found by eRah.
#' @usage plotProfile(Experiment,AlignId, per.class = T, xlim = NULL)
#' @param Experiment A 'MetaboSet' S4 object containing the experiment after being deconolved, aligned and (optionally) identified.
#' @param AlignId the Id identificator for the compound to be shown.
#' @param per.class logical. if TRUE the profiles are shown one color per class, if FALSE one color per sample.
#' @param xlim x axsis (retention time) limits (see \code{\link{plot.default}}).
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{plotSpectra}} \code{\link{plotAlign}}
#' @export

plotProfile <- function(Experiment,AlignId, per.class=T, xlim=NULL)
{	
	if(!(any(unlist(lapply(Experiment@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")
	
	if(nrow(Experiment@MetaData@Phenotype)==0) 
		if(per.class==T) {per.class=F; warning("The experiment does not contain phenotypic metadata. Profiles are shown per sample.")}
	
	empty.samples <- which(lapply(Experiment@Data@FactorList,nrow)==0)
	if(length(empty.samples)!=0)
	{
		Experiment@Data@FactorList <- Experiment@Data@FactorList[-empty.samples]
		Experiment@MetaData@Phenotype <- Experiment@MetaData@Phenotype[-empty.samples,]
	}
	
	alignId <- lapply(Experiment@Data@FactorList,function(x){x$AlignID})
	N.groups <- max(unique(unlist(alignId)))	
	N.samples <- length(Experiment@Data@FactorList)
	
	samples.name <- names(Experiment@Data@FactorList)		
	for(i in 1:N.samples) samples.name[i] <- strsplit(as.character(samples.name[i]), split="\\.")[[1]][1]
	
	profile.list <- lapply(Experiment@Data@FactorList,function(x) {
			outp <- as.character(x[which(x$AlignID==AlignId),"Profile"])
			time <- NA
			int <- NA
			if(length(outp)!=0)
			{
				output <- sparse.to.vector(outp)
				time <- output$time
				int <- output$int*as.numeric(as.character(x[which(x$AlignID==AlignId),"Peak Height"]))
			}
			list(time=time,int=int)
			})
	
	profile.len <- max(unlist(lapply(unlist(profile.list, recursive=F),length)))
	
	profile.time = matrix(NA, nrow=profile.len,ncol=N.samples)
	profile.int = matrix(NA, nrow=profile.len,ncol=N.samples)

	for(i in 1:N.samples) {profile.time[1:length(profile.list[[i]]$time),i] <- profile.list[[i]]$time; profile.int[1:length(profile.list[[i]]$int),i] <- profile.list[[i]]$int}
	
	na.samples.i <- which(apply(apply(profile.int,2,is.na),2,all)==T)
	na.samples.t <- which(apply(apply(profile.time,2,is.na),2,all)==T)
	na.samples <- unique(na.samples.i,na.samples.t)
	
	if(length(na.samples)!=0)
	{ 
		samples.name <- samples.name[-na.samples]
		profile.int <- profile.int[,-na.samples]
		profile.time <- profile.time[,-na.samples]
	}
	
	vector.time <- diag(profile.time[apply(profile.int,2,which.max),])
	time.mean <- mean(vector.time)
	
	profile.time <- sweep(profile.time,2,(vector.time-time.mean),"-") 
	
	compound.name <- as.character(Experiment@Results@Identification[which(Experiment@Results@Identification$AlignID==AlignId),"Name"])
	
	if(per.class==F)
	{
		
		matplot(profile.time, profile.int, type="l", lty=1, col=(1:length(samples.name)), main=paste("Profile Comparison \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.name, pch=19, col=(1:length(samples.name)), title="Samples")
		par(font=1)
	}else{
		pn <- Experiment@MetaData@Phenotype
		indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
		class.names <- pn[indx,"class"]
		
		samples.class.type <- levels(pn$class)
			
		matplot(profile.time,profile.int, type="l", lty=1, col=class.names, main=paste("Profile Comparison \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.class.type, pch=19, col=1:length(samples.class.type), title="Classes")
		par(font=1)

	}
	
}

#' @name plotAlign
#' @aliases plotAlign
#' @title Plotting chromatographic profile with and without alignment
#' @description Plots the chromatophic profiles of the compounds found by eRah. Similarly to plotProfile, but with two sub-windows, showing the chromatophic profiles before and after alignment.
#' @usage plotAlign(Experiment,AlignId, per.class = T, xlim = NULL)
#' @param Experiment A 'MetaboSet' S4 object containing the experiment after being deconolved, aligned and (optionally) identified.
#' @param AlignId the Id identificator for the compound to be shown.
#' @param per.class logical. if TRUE the profiles are shown one color per class, if FALSE one color per sample.
#' @param xlim x axsis (retention time) limits (see \code{\link{plot.default}}).
#' @author Xavier Domingo-Almenara. xavier.domingo@urv.cat
#' @seealso \code{\link{plotSpectra}} \code{\link{plotProfile}}
#' @export
#' @importFrom graphics par matplot legend

plotAlign <- function(Experiment, AlignId, per.class=T, xlim=NULL)
{	
	if(!(any(unlist(lapply(Experiment@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")
	
	if(nrow(Experiment@MetaData@Phenotype)==0) 
	{
		per.class=F
		warning("The experiment does not contain phenotypic metadata. Profiles are shown per sample.")	
	}
	
	empty.samples <- which(lapply(Experiment@Data@FactorList,nrow)==0)
	if(length(empty.samples)!=0)
	{
		Experiment@Data@FactorList <- Experiment@Data@FactorList[-empty.samples]
		Experiment@MetaData@Phenotype <- Experiment@MetaData@Phenotype[-empty.samples,]
	}
	
	alignId <- lapply(Experiment@Data@FactorList,function(x){x$AlignID})
	N.groups <- max(unique(unlist(alignId)))	
	N.samples <- length(Experiment@Data@FactorList)
	
	samples.name <- names(Experiment@Data@FactorList)		
	for(i in 1:N.samples) samples.name[i] <- strsplit(as.character(samples.name[i]), split="\\.")[[1]][1]
	
	profile.list <- lapply(Experiment@Data@FactorList,function(x) {
			outp <- as.character(x[which(x$AlignID==AlignId),"Profile"])
			time <- NA
			int <- NA
			if(length(outp)!=0)
			{
				output <- sparse.to.vector(outp)
				time <- output$time
				int <- output$int*as.numeric(as.character(x[which(x$AlignID==AlignId),"Peak Height"]))
			}
			list(time=time,int=int)
			})
	
	profile.len <- max(unlist(lapply(unlist(profile.list, recursive=F),length)))
	
	profile.time = matrix(NA, nrow=profile.len,ncol=N.samples)
	profile.int = matrix(NA, nrow=profile.len,ncol=N.samples)

	for(i in 1:N.samples) {profile.time[1:length(profile.list[[i]]$time),i] <- profile.list[[i]]$time; profile.int[1:length(profile.list[[i]]$int),i] <- profile.list[[i]]$int}
	
	na.samples <- which(apply(apply(profile.time,2,is.na),2,all)==T)
	if(length(na.samples)!=0)
	{ 
		samples.name <- samples.name[-na.samples]
		profile.int <- profile.int[,-na.samples]
		profile.time <- profile.time[,-na.samples]
	}
	
	vector.time <- diag(profile.time[apply(profile.int,2,which.max),])
	time.mean <- mean(vector.time)
	
	profile.un.time <- profile.time
	profile.time <- sweep(profile.time,2,(vector.time-time.mean),"-") 
	
	compound.name <- as.character(Experiment@Results@Identification[which(Experiment@Results@Identification$AlignID==AlignId),"Name"])
	
	if(per.class==F)
	{	
		par(mfrow=c(1,2))
		matplot(profile.un.time, profile.int, type="l", lty=1, col=(1:length(samples.name)), main=paste("Unaligned \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.name, pch=19, col=(1:length(samples.name)), title="Samples")
		par(font=1)
		matplot(profile.time, profile.int, type="l", lty=1, col=(1:length(samples.name)), main=paste("Aligned \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.name, pch=19, col=(1:length(samples.name)), title="Samples")
		par(font=1)
	}else{
		pn <- Experiment@MetaData@Phenotype
		indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
		class.names <- pn[indx,"class"]
		
		samples.class.type <- levels(pn$class)
		
		par(mfrow=c(1,2))
		matplot(profile.un.time, profile.int, type="l", lty=1, col=class.names, main=paste("Unaligned \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.class.type, pch=19, col=1:length(samples.class.type), title="Classes")
		par(font=1)
		matplot(profile.time, profile.int, type="l", lty=1, col=class.names, main=paste("Aligned \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.class.type, pch=19, col=1:length(samples.class.type), title="Classes")
		par(font=1)
	}
	
}

#' @name plotChr
#' @aliases plotChr
#' @title Plotting sample chromatogram
#' @description Plot the sample chromatogram
#' @usage plotChr(Experiment, N.sample = 1, type = c("BIC","TIC","EIC"), xlim = NULL, mz = NULL)
#' @param Experiment A 'MetaboSet' S4 object containing the experiment.
#' @param N.sample Integer. The number of the sample to query.
#' @param type The type of plotting, Base Ion Chromatogram (BIC), Total Ion Chromatogram (TIC), or Extracted Ion Chromatogram (EIC).
#' @param xlim The range in minutes, separated by comas: c(rt.min, rt.max) of the limits of plotting. By default, all the chromatogram is plotted.
#' @param mz Just when EIC is selected. The range separated by comas: c(mz.min, mz.max) or a vector of numbers: c(50,67,80), of the masses to be ploted.
#' @examples \dontrun{
#' # First, an experiment has to be already created by newExp()
#' # then, each sample chromatogram can be plotted by:
#'
#' plotChr(Experiment, 1, "BIC")
#' plotChr(Experiment, 1, "TIC", xlim=c(5,7))  #Plots from
#' # minute 5 to 7.
#' 
#' plotChr(Experiment, 1, "EIC", mz=50:70 xlim=c(5,7))  #Plots from
#' # minute 5 to 7, and only the masses from 50 to 70.
#'
#' plotChr(Experiment, 1, "EIC", xlim=c(7,7.5), mz=c(50,54,70))  #Plots
#' # the EIC from minute 7 to 7.5, and only the masses 50, 54 and 70.
#' }
#' @seealso \code{\link{sampleInfo}}
#' @export
#' @importFrom graphics plot

plotChr <- function(Experiment, N.sample=1, type=c("BIC","TIC","EIC"), xlim=NULL, mz=NULL)
{
	type <- match.arg(type, c("BIC","TIC","EIC"), several.ok = FALSE)
	
	if(Experiment@MetaData@DataDirectory == "") {
        filename <- as.character(Experiment@MetaData@Instrumental$filename[N.sample])
    }else{
        filename <- paste(Experiment@MetaData@DataDirectory, "/", Experiment@MetaData@Instrumental$filename[N.sample], sep = "")
    }
	
	sampleRD <- load.file(filename)
	
	max.rt <- (nrow(sampleRD@data)/(sampleRD@scans.per.second*60)) + sampleRD@start.time/60
	min.rt <- sampleRD@start.time/60
	vect.rt <- seq(min.rt, max.rt, length.out=nrow(sampleRD@data))
	
	if(type=="TIC") plot(vect.rt, rowSums(sampleRD@data), type="l", main="TIC", xlab="RT(min)", ylab="TIC", xlim)
	if(type=="BIC") plot(vect.rt, apply(sampleRD@data,1,max), type="l", main="BIC", xlab="RT(min)", ylab="BIC", xlim)	
	if(type=="EIC") 
	{
		if(is.null(mz)) stop("When plotting the EIC, please, specify the masses or range of masses to be plotted.")
		if(any(mz<sampleRD@min.mz) || any(mz>sampleRD@max.mz)) stop(paste("The adquisition m/z range selected is above or under the selected masses. For this sample, the adquistion range is from ", sampleRD@min.mz, " to ", sampleRD@max.mz, ". Please, change the m/z paramter accordingly, inside the limits.", sep=""))

		mz.rang <- mz - (sampleRD@min.mz - 1)
	    if(!is.null(xlim)) scan.range <- c(which.min(abs(vect.rt - min(xlim))):which.min(abs(vect.rt - max(xlim))))
        if(is.null(xlim)) scan.range <- 1:nrow(sampleRD@data)
		matplot(vect.rt[scan.range], sampleRD@data[scan.range,mz.rang], type="l", lty=1, main="EIC", xlab="RT(min)", ylab="EIC")	
		legend("topright",legend=mz, pch=19, col=1:length(mz), title="m/z")

	}
}

# plotBoxplot <- function(Experiment, AlignId, classes.to.compare=NULL, outline = TRUE, log = "", horizontal = FALSE)
# {
	# if(!(any(unlist(lapply(Experiment@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")
	# if(nrow(Experiment@MetaData@Phenotype)==0) stop("No Phenotype data has been attached to this experiment.")
	
	# align.list <- NULL
	# align.list <- alignList(Experiment)
	# if(is.null(align.list)) return()
	
	# empty.samples <- which(lapply(Experiment@Data@FactorList,nrow)==0)
	# if(length(empty.samples)!=0)
	# {
		# Experiment@Data@FactorList <- Experiment@Data@FactorList[-empty.samples]
		# Experiment@MetaData@Phenotype <- Experiment@MetaData@Phenotype[-empty.samples,]
	# }
	
	# pn <- Experiment@MetaData@Phenotype
	# samples.name <- names(Experiment@Data@FactorList)		
	# indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
	# class.names <- as.vector(pn[indx,"class"])
		
	# samples.class.type <- levels(as.factor(class.names))
	
	# data.list <- align.list[,-c(1:4)]

	# if(!is.null(classes.to.compare))
	# {
		# classes.selected <- unique(unlist(apply(as.matrix(classes.to.compare),1,function(x) which(x==samples.class.type ))))
		# if(length(classes.selected)!=length(classes.to.compare)) stop("Invalid class selected. The class name is invalid or no sample of this class has been processed. Use expClasses() function to obtain more information about the existing classes in this experiment")
		# sel.classes.index <- unlist(apply(as.matrix(classes.to.compare),1,function(x) which(x==class.names)))
		# class.names <- class.names[sel.classes.index]
		# data.list <- data.list[,sel.classes.index]
	# }
	
	# p.value <- round(summary(aov(as.numeric(as.vector(as.matrix(data.list[which(align.list[,"AlignID"]==AlignId),])))~as.factor(class.names)))[[1]][1,5], digits=4)
	# main.title <- paste(align.list[which(align.list[,"AlignID"]==AlignId),"Factor"], "\n p-value = ", round(p.value, digits=4), sep="")
	
	# boxplot(as.numeric(as.vector(as.matrix(data.list[which(align.list[,"AlignID"]==AlignId),])))~as.factor(class.names), outline = outline, log = log, horizontal = horizontal, main=main.title, ylab="Peak Height")

# }

# plotPca <- function(Experiment)
# {
	# if(nrow(Experiment@MetaData@Phenotype)==0) stop("No Phenotype data has been attached to this experiment.")
	
	# align.list <- NULL
	# align.list <- alignList(Experiment)
	# if(is.null(align.list)) return()
	
	# empty.samples <- which(lapply(Experiment@Data@FactorList,nrow)==0)
	# if(length(empty.samples)!=0)
	# {
		# Experiment@Data@FactorList <- Experiment@Data@FactorList[-empty.samples]
		# Experiment@MetaData@Phenotype <- Experiment@MetaData@Phenotype[-empty.samples,]
	# }
	
	# pn <- Experiment@MetaData@Phenotype
	# samples.name <- names(Experiment@Data@FactorList)		
	# indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
	# class.names <- as.vector(pn[indx,"class"])
		
	# samples.class.type <- levels(as.factor(class.names))
	
	# data.list <- align.list[,-c(1:3)]
	
	# PCA.analysis <- prcomp(matrix(as.numeric(as.matrix(data.list)), ncol=ncol(data.list)), scale=F, center=T)
	
	# x.lab <- paste("Scores on PC 1 (", round(summary(PCA.analysis)$importance[2,1]*100, digits=1), "%)", sep="")
	# y.lab <- paste("Scores on PC 2 (", round(summary(PCA.analysis)$importance[2,2]*100, digits=1), "%)", sep="")
	# plot(PCA.analysis$rotation[,1],PCA.analysis$rotation[,2], pch=20, col=as.factor(class.names), xlab=x.lab, ylab=y.lab)
	# legend("topright",legend=samples.class.type, pch=19, col=1:length(samples.class.type), title="Classes")
# }