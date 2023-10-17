#' Get the estimates from the SMART-covariate model
#'
#' @param model the model output from SMART_covariate().
#' @param stData the stData used in SMART_covariate().
#' @param covarDat the covarDat used in SMART_covariate().
#' @param covar a character string selecting the covariate estimate you would like to get (one covariate at a time).
#'
#' @return a cell type proportion matrix
#' @return a cell type-specific gene expression matrix
#' @import keyATM
#' @import quanteda
#' @export
#'
#' @examples
#' \dontrun{
#' library(SMART)
#' data(MPOA)
#' mod <- SMART_covariate(stData=MPOA$stMat, markerGs=MPOA$markers,
#'                        covarDat=MPOA$covDat,covars='sex',
#'                        noMarkerCts=1,outDir='SMART_results', seed=1)
#' res <- get_est(model=mod$model, stData=MPOA$stMat,
#'                covarDat=MPOA$covDat, covar='sex')
#' }
get_est <- function(model, stData, covarDat, covar){
  stData <- keyATM_read(texts = as.dfm(stData),keep_docnames = T)
  strata_tw <- by_strata_TopicWord(model, stData,
                                   by = as.vector(unlist(covarDat[,covar])))

  colnames(strata_tw$theta)=gsub('[0-9]{1,2}_(.*)','\\1',colnames(strata_tw$theta))
  for(i in 1:length(strata_tw$phi)){
    rownames(strata_tw$phi[[i]]$phi)=gsub('[0-9]{1,2}_(.*)','\\1',rownames(strata_tw$phi[[i]]$phi))
  }

  return(list(ct_proportions = strata_tw$theta, ct_spec_gexp = strata_tw$phi))
}
