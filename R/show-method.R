#' show-MetaboSet
#' @title Show MetaboSet object
#' @description Show MetaboSet object
#' @param object S4 object of class MetaboSet

setMethod("show", "MetaboSet", function(object){
  cat("A \"MetaboSet\" object containing", length(object@Data@FactorList), "samples \n \n" ,sep = " ")
  cat("Data processed with", object@Data@Parameters$algorithm, "\n" ,sep = " ")
  cat("Info attached to this experiment: \n", object@Info)
})

# setMethod("show", "expClasses", function(object) {
# classes.string <- paste(object@classes.type, collapse=", ")
# cat("Experiment containing ", nrow(object@classes.summary), " samples in ", length(object@classes.type), " different type of classes named: ",classes.string, ". \n \n", sep="")
# print(object@classes.summary)
# })