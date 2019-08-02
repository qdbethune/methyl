# -- This version of cpgCollapse from minfi fixes an issue with its use of the getCN function
#    Fixed by: Rafael Irizarry 

new_cpgCollapse <- function (object, what = c("Beta", "M"), maxGap = 500, blockMaxGap = 2.5 * 
            10^5, maxClusterWidth = 1500, dataSummary = colMeans, na.rm = FALSE, 
          returnBlockInfo = TRUE, islandAnno = NULL, verbose = TRUE, 
          ...) 
{
  what <- match.arg(what)
  if (verbose) 
    message("[cpgCollapse] Creating annotation.\n")
  islands <- minfi:::.getIslandAnnotation(object = object, islandAnno = islandAnno)
  relationToIsland <- islands$Relation_to_Island
  islandName <- islands$Islands_Name
  gr <- granges(object)
  anno <- minfi:::cpgCollapseAnnotation(gr = gr, relationToIsland = relationToIsland, 
                                islandName = islandName, maxGap = maxGap, blockMaxGap = blockMaxGap, 
                                maxClusterWidth = maxClusterWidth, verbose = verbose)
  Indexes <- split(seq_along(anno$pns), anno$pns)
  if (verbose) 
    message("[cpgCollapse] Collapsing data")
  meth_signal <- getMethSignal(object, what = what, ...)
  collapsed_meth_signal <- minfi:::.cpgCollapse(x = meth_signal, Indexes = Indexes, 
                                        dataSummary = dataSummary, na.rm = na.rm, verbose = verbose)
  cn <- getCN(object, ...)
  if (!is.null(cn)) {
    collapsed_cn <- .cpgCollapse(x = cn, Indexes = Indexes, 
                                 dataSummary = dataSummary, na.rm = na.rm, verbose = verbose)
  } else collapsed_cn <- NULL
  preproc <- c(collapse = "cpgCollapse", preprocessMethod(object))
  if (what == "M") {
    ret <- GenomicRatioSet(gr = anno$anno, Beta = NULL, M = collapsed_meth_signal, 
                           CN = collapsed_cn, colData = colData(object), annotation = annotation(object), 
                           preprocessMethod = preproc)
  }
  else {
    ret <- GenomicRatioSet(gr = anno$anno, Beta = collapsed_meth_signal, 
                           M = NULL, CN = collapsed_cn, colData = colData(object), 
                           annotation = annotation(object), preprocessMethod = preproc)
  }
  anno <- anno[2:3]
  if (returnBlockInfo) {
    return(list(object = ret, blockInfo = anno))
  }
  ret
}
