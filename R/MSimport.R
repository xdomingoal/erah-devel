get.compound.info_GMD <- function(k, Spl.List, type)
{	
		#kmain <<-k
		x <- Spl.List[[k]]

		x.split <- unlist(apply(as.matrix(x),1, function(x) {strsplit(x,":")[[1]]}), recursive=T)
	
		Name=" "
		Synon=" "
		RI.VAR5.FAME=0
		RI.VAR5.ALK=0
		RI.MDN35.FAME=0
		RI.MDN35.ALK=0
		CAS=" "
		Formula=" "
		InChi=" "
		KEGG=""
		MW=0
		Comment=" "
		GMD.LINK=" "
		GMD.VERS=" "
		
		Spectra <- ""
		Caps <- apply(as.matrix(x),1, function(x) {strsplit(x,":")[[1]]})
		Sp.Ini <- which(lapply(Caps,function(y)y[1])=="Num Peaks")
		Spectra <- paste(x[(Sp.Ini+1):length(x)],collapse="")
		
		for(j in 1:length(x.split))
		{
			if(x.split[j]==" MST N") Name=x.split[j+1]
			if(x.split[j]==" METB N") {
				if(Synon!=" ") Synon=paste(Synon, "//", x.split[j+1], sep=" ")
				if(Synon==" ") Synon = x.split[j+1]
			}
			if(x.split[j]==" RI")
			{
				if(type=="MDN35.ALK") RI.MDN35.ALK=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
				if(type=="MDN35.FAME") RI.MDN35.FAME=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
				if(type=="VAR5.ALK") RI.VAR5.ALK=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
				if(type=="VAR5.FAME") RI.VAR5.FAME=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
			}
			if(x.split[j]==" MST SEL MASS") SelMZ <- as.vector(sapply(strsplit(x.split[j+1], "\\|"), as.numeric))
			if(x.split[j]=="CAS#") CAS=x.split[j+1]
			if(x.split[j]=="Formula") Formula=x.split[j+1]
			if(x.split[j]==" METB InChI") InChi= strsplit(x.split[j+1], "=")[[1]][2]
			if(x.split[j]==" METB KEGG") KEGG=x.split[j+1]			
			if(x.split[j]=="MW")
			{
				no.dot <- apply(as.matrix(x.split[j+1]), 2, gsub, patt="\\.", replace="")
				MW=as.numeric(apply(as.matrix(no.dot), 2, gsub, patt=",", replace="."))
			}
			if(x.split[j]=="Comment") Comment=x.split[j+1]
			if(x.split[j]==" GMD LINK") GMD.LINK=paste(x.split[j+1],":",x.split[j+2],sep="")
			if(x.split[j]==" GMD VERS") GMD.VERS= paste(x.split[j+1],":",x.split[j+2],sep="")		
		}
		#cat(Name, "\n")
		s.spl <- strsplit(Name," ")[[1]]; Name <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(Synon," ")[[1]]; Synon <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(CAS," ")[[1]]; CAS <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(Formula," ")[[1]]; Formula <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(KEGG," ")[[1]]; KEGG <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(Comment," ")[[1]]; Comment <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(GMD.LINK," ")[[1]]; GMD.LINK <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(GMD.VERS," ")[[1]]; GMD.VERS <- paste(s.spl[-1], collapse=" ")
			
		compound.info <- list(Name=Name, Synon=Synon, RI.VAR5.FAME=RI.VAR5.FAME, RI.VAR5.ALK=RI.VAR5.ALK, RI.MDN35.FAME=RI.MDN35.FAME, RI.MDN35.ALK=RI.MDN35.ALK, SelMZ=SelMZ, CAS=CAS, InChi=InChi, Formula=Formula, MW=MW, KEGG=KEGG, Comment=Comment, GMD.LINK=GMD.LINK,	GMD.VERS=GMD.VERS, Spectra=Spectra)
		
	compound.info
}

get.compound.info_MSP <- function(k, Spl.List, type)
{	
		#kmain <<-k
		x <- Spl.List[[k]]

		x.split <- unlist(apply(as.matrix(x),1, function(x) {strsplit(x,":")[[1]]}), recursive=T)
	
		Name=" "
		Synon=" "
		RI.VAR5.FAME=0
		RI.VAR5.ALK=0
		RI.MDN35.FAME=0
		RI.MDN35.ALK=0
		CAS=" "
		Formula=" "
		InChi=" "
		KEGG=""
		MW=0
		Comment=" "
		GMD.LINK=" "
		GMD.VERS=" "
		SelMZ=NA
		
		Spectra <- ""
		Caps <- apply(as.matrix(x),1, function(x) {strsplit(x,":")[[1]]})
		Sp.Ini <- which(lapply(Caps,function(y)y[1])=="Num Peaks" | lapply(Caps,function(y)y[1])=="Num peaks")
		Spectra.semicolon <- paste(x[(Sp.Ini+1):length(x)],collapse="")
		Spectra.semicolon <- gsub(' {2,}', ' ', Spectra.semicolon)
		
		if(strsplit(Spectra.semicolon, '')[[1]][1]==' ') Spectra.semicolon <- paste(strsplit(Spectra.semicolon, '')[[1]][-1], collapse='')
				
		if(length(grep(";", Spectra.semicolon))!=0){
		 	Spectra.dot <- sapply(strsplit(Spectra.semicolon, ";")[[1]], function(y) paste(strsplit(y, " ")[[1]], collapse=":"))
		}else{
			Spectra.dot <- sapply(x[(Sp.Ini+1):length(x)], function(y) paste(strsplit(y, " ")[[1]], collapse=":"))
		}
		
		Spectra.dot <- sapply(Spectra.dot, function(y) {
			yV <- y
			nVs <- strsplit(y[1], '')[[1]] %in% c(0:9)
			if(nVs[1]==FALSE) yV <- paste(strsplit(y, '')[[1]][-1], collapse='')
			if(nVs[length(nVs)]==FALSE) yV <- paste(strsplit(y, '')[[1]][-length(nVs)], collapse='')
			yV
			})
		
		Spectra <- paste(Spectra.dot , collapse=" ")
		
		for(j in 1:length(x.split))
		{
			if(x.split[j]=="Name") Name=x.split[j+1]
			#if(x.split[j]==" MST SEL MASS") SelMZ <- as.vector(sapply(strsplit(x.split[j+1], "\\|"), as.numeric))
			if(x.split[j]=="CAS#" | x.split[j]=="CASNO") CAS=x.split[j+1]
			if(x.split[j]=="Formula") Formula=x.split[j+1]
			if(x.split[j]=="DB#") InChi= x.split[j+1]
			#if(x.split[j]==" METB KEGG") KEGG=x.split[j+1]		
			
			if(x.split[j]=="RI")
			{
				if(type=="MDN35.ALK") RI.MDN35.ALK=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
				if(type=="MDN35.FAME") RI.MDN35.FAME=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
				if(type=="VAR5.ALK") RI.VAR5.ALK=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
				if(type=="VAR5.FAME") RI.VAR5.FAME=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
			}
				
			if(x.split[j]=="MW")
			{
				no.dot <- apply(as.matrix(x.split[j+1]), 2, gsub, patt="\\.", replace="")
				MW=as.numeric(apply(as.matrix(no.dot), 2, gsub, patt=",", replace="."))
			}
			if(x.split[j]=="Comments" | x.split[j]=="Comment") Comment=x.split[j+1]
			#if(x.split[j]==" GMD LINK") GMD.LINK=paste(x.split[j+1],":",x.split[j+2],sep="")
			#if(x.split[j]==" GMD VERS") GMD.VERS= paste(x.split[j+1],":",x.split[j+2],sep="")		
		}
		#cat(Name, "\n")
		s.spl <- strsplit(Name," ")[[1]]; Name <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(Synon," ")[[1]]; Synon <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(CAS," ")[[1]]; CAS <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(Formula," ")[[1]]; Formula <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(KEGG," ")[[1]]; KEGG <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(Comment," ")[[1]]; Comment <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(GMD.LINK," ")[[1]]; GMD.LINK <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(GMD.VERS," ")[[1]]; GMD.VERS <- paste(s.spl[-1], collapse=" ")
			
		compound.info <- list(Name=Name, Synon=Synon, RI.VAR5.FAME=RI.VAR5.FAME, RI.VAR5.ALK=RI.VAR5.ALK, RI.MDN35.FAME=RI.MDN35.FAME, RI.MDN35.ALK=RI.MDN35.ALK, SelMZ=SelMZ, CAS=CAS, InChi=InChi, Formula=Formula, MW=MW, KEGG=KEGG, Comment=Comment, GMD.LINK=GMD.LINK,	GMD.VERS=GMD.VERS, Spectra=Spectra)
		
	compound.info
}

list.DB <- function(DB.object)
{
		k <- 1
		Cont <- 1
		Spl.List <- list()
		for(i in 1:length(DB.object))
		{
			if(DB.object[i]=="") 
			{	
				Spl.List[Cont] <- list(Compound=DB.object[k:(i-1)])
				k <- i + 1
				Cont <- Cont + 1
			}
		}
		Spl.List
}
	

importGMD <- function(filename, DB.name, DB.version, DB.info, type=c("VAR5.ALK","VAR5.FAME","MDN35.ALK", "MDN35.FAME"))
{
	DB.MSP <- readLines(filename)
		
	Spl.List <- list.DB(DB.MSP)
		
	import.database <- list()
	import.database <- lapply(1:length(Spl.List), function(x) get.compound.info_GMD(x, Spl.List, type))
	
	final.database <- new("eRah_DB",  name=DB.name, version=DB.version, info=DB.info, database=import.database)
	final.database
}	

importMSP <- function(filename, DB.name, DB.version, DB.info)
{
	DB.MSP <- readLines(filename)
		
	Spl.List <- list.DB(DB.MSP)
		
	import.database <- list()
	OutFromList <- unlist(lapply(Spl.List, function(x) all(unlist(x)=="")))
	Spl.List <- Spl.List[!OutFromList]
	import.database <- lapply(1:length(Spl.List), function(x) get.compound.info_MSP(x, Spl.List, type="VAR5.ALK"))
	
	final.database <- new("eRah_DB",  name=DB.name, version=DB.version, info=DB.info, database=import.database)
	final.database
}

