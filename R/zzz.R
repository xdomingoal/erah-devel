
.onAttach <- function(libname, pkgname) {
	        
	data <- try(XML::xmlParse("http://cran.r-project.org/package=erah", isHTML=T), silent=T)
	currVersion <- as.character(utils::packageVersion("erah"))
	
	msg <- paste("Welcome to eRah. For bugs, problems and issues, please do not hesitate in contacting xdomingo@scripps.edu and describing your problem.", sep="") 
    #packageStartupMessage(msg)   
	
	if(class(data)[1]!="try-error")
	{
		xml_data <- XML::xmlToList(data)
		cranVersion <- xml_data$body$table[[1]][[3]]
		#vers_state <- (cranVersion!=currVersion)
		vers_state <- utils::compareVersion(cranVersion, currVersion)
		if(vers_state==-1) vers_state <- 0
		wMsg <- paste("The current version of eRah (", currVersion, ") is outdated. There is a new version of eRah (", cranVersion, "), available at CRAN. To update your version execute the following: \n install.packages('erah')", sep="")
		if(vers_state) warning(wMsg)
	}else{
		wMsg <- paste("The current available version of eRah in CRAN cannot be checked. Please, make sure that your current version of eRah (", currVersion, ") is the same as in CRAN (http://cran.r-project.org/package=erah). To update your version execute the follwing: \n install.packages('erah')", sep="")
		warning(wMsg)
	}       	        
	 	
}    
