
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

#' @rdname sampleInfo

setGeneric('sampleInfo',function(Experiment, N.sample=1){
  standardGeneric('sampleInfo')
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

#' @rdname plotSpectra

setGeneric('plotSpectra',function(Experiment, AlignId, n.putative=1, compare=T, id.database=mslib, comp.db=NULL, return.spectra=F, draw.color="purple", xlim=NULL){
  standardGeneric('plotSpectra')
})

#' @rdname plotProfile

setGeneric('plotProfile',function(Experiment,AlignId, per.class=T, xlim=NULL){
  standardGeneric('plotProfile')
})

#' @rdname plotAlign

setGeneric('plotAlign',function(Experiment, AlignId, per.class=T, xlim=NULL){
  standardGeneric('plotAlign')
})

#' @rdname plotChr

setGeneric('plotChr',function(Experiment, N.sample=1, type=c("BIC","TIC","EIC"), xlim=NULL, mz=NULL){
  standardGeneric('plotChr')
})