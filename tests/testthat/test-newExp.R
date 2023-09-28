
context('experiment')

library(gcspikelite)
data("targets")

files <- list.files(system.file('data',package = 'gcspikelite'),full.names = TRUE)
files <- files[sapply(files,grepl,pattern = 'CDF')][1:2]

instrumental <- createInstrumentalTable(files)

phenotype <- createPhenoTable(files, cls = as.character(targets$Group[order(targets$FileName)])[1:2])

ex <- newExp(instrumental = instrumental,phenotype = phenotype, info = "DEMO Experiment")

ex.dec.par <- setDecPar(min.peak.width = 1,avoid.processing.mz = c(35:69,73:75,147:149))
ex.al.par <- setAlPar(min.spectra.cor = 0.90, max.time.dist = 3, mz.range = 70:600)

deconvolvedEx <- deconvolveComp(ex,ex.dec.par)
alignedEx <- alignComp(deconvolvedEx,alParameters = ex.al.par)
recoveredEx <-   recMissComp(alignedEx,min.samples = 1)
identifiedEx <- identifyComp(recoveredEx)

test_that('createInstrumentalTable works',{
  expect_true(identical(class(instrumental),c('tbl_df','tbl','data.frame')))
  expect_true(nrow(instrumental) == 2)
  expect_true(ncol(instrumental) == 4)
})

test_that('createPhenoTable works',{
  expect_true(identical(class(phenotype),c('tbl_df','tbl','data.frame')))
  expect_true(nrow(phenotype) == 2)
  expect_true(ncol(phenotype) == 2)
})

test_that('newExp works',{
  expect_true(class(ex) == 'MetaboSet')
})

test_that('setDecPar works',{
  expect_true(class(ex.dec.par) == 'eRahSoftParameters')
})

test_that('setAlPar works',{
  expect_true(class(ex.al.par) == 'eRahAlParameters')
})

test_that('deconvolveComp works',{
  expect_true(class(deconvolvedEx) == 'MetaboSet')
})

test_that('alignComp works',{
  expect_true(class(alignedEx) == 'MetaboSet')
})

test_that('recMissComp works',{
  expect_true(class(recoveredEx) == 'MetaboSet')
})

test_that('identifyComp works',{
  expect_true(class(identifiedEx) == 'MetaboSet')
})

test_that('plotProfile works', {
  expect_no_error(plotProfile(identifiedEx, 84))
})

test_that('plotSpectra works', {
  expect_no_error(plotSpectra(identifiedEx, 84))
})
