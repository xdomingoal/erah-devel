#' @name compInfo
#' @aliases compInfo
#' @title Information of a Compound
#' @description Displays basic information of a compound in the MS library.
#' @usage compInfo(comp.id, id.database = mslib)
#' @param comp.id The DB.Id number of the compound.
#' @param id.database The mass-spectra library to be compared with the empirical spectra. By default, the MassBank - Mass Bank of North America (MoNa) database are employed (mslib object).
#' @details Returns details on a given compound such as the synonyms, CAS, KEGG, retention index, among others.
#' @examples 
#' # finding proline
#' findComp("proline")
#'
#' # we see that proline 2TMS has the DB.Id number 42, then:
#' compInfo(42)
#' @seealso \code{\link{findComp}}
#' @export

compInfo <- function(comp.id, id.database=mslib)
{
	comp.info <- id.database@database[[comp.id]]
	comp.info$Synon <- gsub("//",", ", comp.info$Synon)
	cat("Name: ", comp.info$Name, "\n", "Synonyms: ", comp.info$Synon, "\n", "CAS: ", comp.info$CAS, "\n" ,"Formula: ", comp.info$Formula, "\n" ,"MW: ", comp.info$MW, "\n", "KEGG: ", comp.info$KEGG, "\n", "RI (FAME & Var5): ", comp.info$RI.VAR5.FAME, "\n","RI (ALK & Var5): ", comp.info$RI.VAR5.ALK, "\n","RI (FAME & MDN35): ", comp.info$RI.MDN35.FAME, "\n","RI (ALK & MDN35): ", comp.info$MDN35.ALK, "\n", "------------------------------- \n" ,"Comment: ", comp.info$Comment, sep="")
}

#' @name findComp
#' @aliases findComp
#' @title Find a compound
#' @description Finds compounds in the MS library by Name, CAS or chemical formula.
#' @usage findComp(name = NULL, id.database = mslib, CAS = NULL, chem.form = NULL)
#' @param name The name of the compound to be found.
#' @param id.database The mass-spectra library to be compared with the empirical spectra. By default, the MassBank - Mass Bank of North America (MoNa) database are employed (mslib object).
#' @param CAS The CAS number of the compound to be found.
#' @param chem.form The chemical formula of the compound to be found.
#' @return 
#' \code{findComp} returns an S3 object:
#'      \item{DB.Id}{The identification number of the library. Each metbolite in the reference library has a different DB.Id number.}
#'      \item{Compound Name}{Compound Name.}
#'      \item{CAS}{CAS number}
#'      \item{Formula}{Chemical Formula.}
#' @examples 
#' # finding proline
#'
#' findComp("proline")
#'
#' # be careful, exact matches are not supported, 
#' # as well as different names like these cases:
#'
#' findComp("L-proline (2TMS)")
#'
# or
#'
#' findComp("proline 2")
#' @seealso \code{\link{compInfo}}
#' @export

findComp <- function(name=NULL, id.database=mslib, CAS=NULL, chem.form=NULL)
{
	#Only one argument is allowed. If name is introduced, CAS and form, is depreciated, if CAS is introduced, chemical formula is depreciated.
	if(all(is.null(c(name,CAS,chem.form)))) stop("One argument is needed: Name, CAS or formula")
	if(length(which(c(is.null(name),is.null(CAS),is.null(chem.form)))==FALSE)<2) warning("Only one argument will be used! (Name prevails over CAS, and CAS prevails over Chemical Formula)")
	
	db.form <- unlist(lapply(id.database@database, function(x) x$Formula))
	db.cas <- unlist(lapply(id.database@database, function(x) x$CAS))
	db.names <- unlist(lapply(id.database@database, function(x) x$Name))

	if(!is.null(chem.form)) indexes <- grep(chem.form, db.form, ignore.case=T)	
	if(!is.null(CAS)) indexes <- grep(CAS, db.cas, ignore.case=T)	 
	if(!is.null(name)) indexes <- grep(name, db.names, ignore.case=T)

	met.list <- as.data.frame(matrix(c(indexes, db.names[indexes], db.cas[indexes], db.form[indexes]), ncol=4))
	colnames(met.list) <- c("DB.Id","Compound Name","CAS","Formula")
	met.list
}





# seekSimilar <- function(comp.id, id.database=mslib, n=10)
# {
	# spect.list <- lapply(id.database@database, function(x) x$Spectra)
	# splitted.spectra.list <- lapply(spect.list, function(x) strsplit(as.character(x),split=" ")[[1]])
	# splitted.spectra.list <- lapply(splitted.spectra.list, function(c.bin){ 
		# out <- unlist(strsplit(c.bin,split=":"))
		# mz.ind <- seq(from=1, to=(length(out)-1),by=2)
		# int.ind <- seq(from=2, to=length(out),by=2)
		# list(mz=out[mz.ind], int=out[int.ind])
		# })
	# maxMz <- max(unlist(lapply(splitted.spectra.list, function(x) max(as.numeric(as.vector(x$mz))))))

	
	# spect.mat <- unlist(lapply(splitted.spectra.list, function(x) {
		# out <- rep(0,maxMz)
		# out[as.numeric(as.vector(x$mz))] <- as.numeric(as.vector(x$int))
		# out
	# }))
	# spect.mat <- matrix(spect.mat, nrow=maxMz)
	# cor.mat <- cor(spect.mat[,comp.id], spect.mat[,-comp.id])
	
	# corr.vector <- as.vector(cor.mat)
	# sim.index <- order(corr.vector, decreasing=T)[1:n]
	# sim.corr <- round(corr.vector[sim.index]*100, digits=2)
	
	# sim.list <- as.data.frame(NULL)

	# sim.names <- unlist(lapply(id.database@database[sim.index], function(x){x$Name}))
	# sim.cas <- unlist(lapply(id.database@database[sim.index], function(x){x$CAS}))
	# sim.form <- unlist(lapply(id.database@database[sim.index], function(x){x$Formula}))
		
	# sim.list <- as.data.frame(matrix(c(sim.index,sim.names,sim.corr,sim.cas,sim.form), nrow=n))
	# colnames(sim.list) <- c("DB.Id","Name","MatchFactor","CAS","Formula")
	# sim.list		
# }