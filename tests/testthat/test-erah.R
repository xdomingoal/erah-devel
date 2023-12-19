
test_that('erah works',{
  skip_if_not_installed('gcspikelite')
  
  library(gcspikelite)
  data('targets')
  
  files <- list.files(system.file('data',package = 'gcspikelite'),full.names = TRUE)
  files <- files[sapply(files,grepl,pattern = 'CDF')][1:2]
  
  instrumental <- createInstrumentalTable(files)
  
  expect_true(identical(class(instrumental),c('tbl_df','tbl','data.frame')))
  expect_true(nrow(instrumental) == 2)
  expect_true(ncol(instrumental) == 4)
  
  phenotype <- createPhenoTable(
    files, 
    cls = as.character(targets$Group[order(targets$FileName)])[1:2]
  )
  
  expect_true(identical(class(phenotype),c('tbl_df','tbl','data.frame')))
  expect_true(nrow(phenotype) == 2)
  expect_true(ncol(phenotype) == 2)
  
  ex <- newExp(instrumental = instrumental,phenotype = phenotype, info = "DEMO Experiment")
  
  expect_true(class(ex) == 'MetaboSet')
  
  ex.dec.par <- setDecPar(min.peak.width = 1,avoid.processing.mz = c(35:69,73:75,147:149))
  expect_true(class(ex.dec.par) == 'eRahSoftParameters')
  
  ex.al.par <- setAlPar(min.spectra.cor = 0.90, max.time.dist = 3, mz.range = 70:600)
  expect_true(class(ex.al.par) == 'eRahAlParameters')
  
  deconvolvedEx <- deconvolveComp(ex,ex.dec.par)
  expect_true(class(deconvolvedEx) == 'MetaboSet')
  
  alignedEx <- alignComp(deconvolvedEx,alParameters = ex.al.par)
  expect_true(class(alignedEx) == 'MetaboSet')
  
  recoveredEx <-   recMissComp(alignedEx,min.samples = 1)
  expect_true(class(recoveredEx) == 'MetaboSet')
  
  identifiedEx <- identifyComp(recoveredEx)
  expect_true(class(identifiedEx) == 'MetaboSet')
  
  expect_no_error(plotProfile(identifiedEx, 84))
  expect_no_error(plotSpectra(identifiedEx, 84))
})
