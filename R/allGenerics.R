
#' @rdname metaData

setGeneric('metaData',function(object){
  standardGeneric('metaData')
})

#' @rdname phenoData

setGeneric('phenoData',function(object){
  standardGeneric('phenoData')
})

#' @rdname expClasses

setGeneric('expClasses',function(object){
  standardGeneric('expClasses')
})

#' @rdname deconvolveComp

setGeneric('deconvolveComp',function(Experiment, decParameters, samples.to.process=NULL, down.sample=FALSE, virtualScansPerSecond=NULL){
  standardGeneric('deconvolveComp')
})

#' @rdname alignComp

setGeneric('alignComp',function(Experiment, alParameters, blocks.size=NULL){
  standardGeneric('alignComp')
})

#' @rdname recMissComp

setGeneric('recMissComp',function(Experiment, min.samples, free.model=F){
  standardGeneric('recMissComp')
})

#' @rdname identifyComp

setGeneric('identifyComp',function(Experiment, id.database=mslib, mz.range=NULL, n.putative=3){
  standardGeneric('identifyComp')
})