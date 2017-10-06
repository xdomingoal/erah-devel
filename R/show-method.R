
setMethod("show", "MetaboSet", function(object){
  cat("A \"MetaboSet\" object containing", length(object@Data@FactorList), "samples \n \n" ,sep=" ")
  cat("Data processed with", object@Data@Parameters$algorithm, "\n" ,sep=" ")
  cat("Info attached to this experiment: \n", object@Info)
})