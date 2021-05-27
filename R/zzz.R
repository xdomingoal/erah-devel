

.onAttach <- function(libname, pkgname) {

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

  packageStartupMessage(erah.logo)
}



