

.onAttach <- function(libname, pkgname) {

  data <- try(XML::xmlParse("http://cran.r-project.org/package=erah", isHTML=T), silent=T)
  currVersion <- as.character(utils::packageVersion("erah"))

  erah.logo <- paste0("\n
      _____        _
     |  __ \\      | |        eRah R package:
  ___| |__) | __ _| |__   	 ----------------
 / _ \\  _   // _` |  _ \\     Untargeted GC-MS metabolomics profiling
|  __/ |  \\ \\ (_| | | | |
 \\___|_|   \\_\\__,_|_| |_|    Version ", currVersion, "

 - Type 'citation('erah')' for citing this R package in publications.
 - Type 'vignette('eRahManual', package='erah')' for a tutorial on eRah's usage.
 - For bugs, problems and issues, please do not hesitate in contacting xavier.domingoa@eurecat.org or opening an issue on the Github repository (https://github.com/xdomingoal/erah-devel/issues).")

  #msg <- paste("Welcome to eRah. For bugs, problems and issues, please do not hesitate in contacting xdomingo@scripps.edu and describing your problem.", sep="")
  packageStartupMessage(erah.logo)
}



