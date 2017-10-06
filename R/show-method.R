#' @name show.MetaboSet
#' @aliases show.MetaboSet
#' @title Show MetaboSet object
#' @description Show MetaboSet object

setMethod("show", "MetaboSet", function(object){
  cat("A \"MetaboSet\" object containing", length(object@Data@FactorList), "samples \n \n" ,sep = " ")
  cat("Data processed with", object@Data@Parameters$algorithm, "\n" ,sep = " ")
  cat("Info attached to this experiment: \n", object@Info)
})