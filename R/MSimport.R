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
    if(x.split[j]==" MST SEL MASS") SelMZ <- as.vector(sapply(strsplit(x.split[j+1], "\\|")[[1]],function(z){
      if (grepl('NA',z)) {
        return(NA)
      } else {
        return(as.numeric(z))
      }
    }))
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
	
#' @name importGMD
#' @aliases importGMD
#' @title Import MSP files from GMD to R
#' @description Import the Golm Metabolome Database.
#' @usage importGMD(filename, DB.name, DB.version, DB.info, 
#' type = c("VAR5.ALK","VAR5.FAME","MDN35.ALK", "MDN35.FAME"))
#' @param filename The filepath containing the GMD database file.
#' @param DB.name The name of the database (each user may chose its own name
#' @param DB.version The version of the database (each user may chose its own version)
#' @param DB.info Some info about the database for further reference
#' @param type The type of RI to be imported from the database
#' @details For more details, please see the eRah manual
#' @export

importGMD <- function(filename, DB.name, DB.version, DB.info, type=c("VAR5.ALK","VAR5.FAME","MDN35.ALK", "MDN35.FAME"))
{
  DB.MSP <- readLines(filename)
  
  Spl.List <- list.DB(DB.MSP)
  
  import.database <- list()
  import.database <- lapply(1:length(Spl.List), function(x) get.compound.info_GMD(x, Spl.List, type))
  
  final.database <- new("eRah_DB",  name=DB.name, version=DB.version, info=DB.info, database=import.database)
  final.database
}	

#' @name importMSP
#' @aliases importMSP
#' @title Import MSP files to R
#' @description Import MS libraries in MSP format to eRah DB format.
#' @usage importMSP(filename, DB.name, DB.version, DB.info)
#' @param filename The filepath containing the MSP library file.
#' @param DB.name The name of the database (each user may chose its own name)
#' @param DB.version The version of the database (each user may chose its own version)
#' @param DB.info Some info about the database for further reference
#' @details 
#' The MSP input file should look like:
#'
#' -----
#'  
#'   Name: Metabolite_name 
#'
#' Formula: H2O 
#'
#' MW: 666
#'
#' ExactMass: 666.266106 
#' 
#' CAS#: 11-22-3 
#'
#' DB#: 1 
#'
#' Comments: Metabolite_name reference standard
#' 
#' Num Peaks: XX
#'
#' 53 1; 54 2; 55 5; 56 2; 57 2; 
#'
#' 58 14; 59 18; 60 1000; 61 2; 67 1; 
#' 
#' Name: Metabolite_name_2
#'
#' Formula: H2O2
#'
#' MW: 999
#'
#' ExactMass: 999.266106
#'
#' CAS#: 22-33-4
#' 
#' DB#: 2
#'
#'Comments: Metabolite_name_"" reference standard
#'
#' Num Peaks: XX
#'
#' 66 10; 67 1000; 155 560; 156 800; 157 2; 
#'
#'158 14; 159 1; 160 100; 161 2; 167 1; 
#'
#' -------	
#'  
#'   OR
#'
#' -----
#'  
#'   Name: Metabolite_name 
#'
#' Formula: H2O 
#'
#' MW: 666
#' 
#' ExactMass: 666.266106 
#'
#' CASNO: 11-22-3 
#'
#' DB#: 1 
#'
#' Comment: Metabolite_name reference standard
#'
#' Num peaks: XX
#'
#' 53 1
#'
#' 54 2
#'
#' 55 5 
#'
#' Name: Metabolite_name_2
#'
#' Formula: H2O2
#'
#' MW: 999
#'
#' ExactMass: 999.266106
#'
#' CASNO: 22-33-4
#'
#' DB#: 2
#'
#' Comment: Metabolite_name_"" reference standard
#' 
#' Num Peaks: XX
#'
#' 66 10
#'
#' 67 1000
#'
#' 155 560
#'
#' -------	
#'  
#'   Or combinations of both.
#'
#'
#' For more details, please see the eRah manual.
#' @export
#' @importFrom methods new

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

