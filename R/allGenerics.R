

#' @rdname plotAlign

setGeneric('plotAlign',function(Experiment, AlignId, per.class=T, xlim=NULL){
  standardGeneric('plotAlign')
})

#' @rdname plotChr

setGeneric('plotChr',function(Experiment, N.sample=1, type=c("BIC","TIC","EIC"), xlim=NULL, mz=NULL){
  standardGeneric('plotChr')
})

#' @rdname idList

setGeneric('idList',function(object, id.database=mslib){
  standardGeneric('idList')
})

#' @rdname alignList

setGeneric('alignList',function(object, by.area=TRUE){
  standardGeneric('alignList')
})

#' @rdname dataList

setGeneric('dataList',function(Experiment, id.database=mslib, by.area=TRUE){
  standardGeneric('dataList')
})