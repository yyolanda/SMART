#' Spatial transcriptomics deconvolution with the covariate model of SMART
#'
#' @param stData a spot(rows)-by-gene(columns) spatial transcriptomics matrix.
#' @param markerGs a list with marker genes for each cell type.
#' @param noMarkerCts the number of cell types without marker genes.
#' @param covarDat a data.frame containing the covariate (column) values of each spot (row).
#' @param covars a character vector specifying the covariates to include.
#' @param outDir the path to the output directory.
#' @param seed the starting value.
#' @param iterations the number of iterations (default: 2000).
#' @param priors do not change the priors if you don't know what they are.
#'
#' @return a model object
#' @import keyATM
#' @import quanteda
#' @importFrom stats as.formula
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
SMART_covariate <- function(stData, markerGs, covarDat, covars,
                            noMarkerCts=1, outDir='SMART_results',
                            seed=1, iterations=2000, priors=NULL){
  # read data
  stData <- keyATM_read(texts = as.dfm(stData), keep_docnames = TRUE)

  # options
  my_options <- list(seed          = seed, # automatically generate random seed
                     iterations    = iterations,
                     verbose       = TRUE,
                     llk_per       = 100,
                     use_weights   = TRUE,
                     weights_type  = "inv-freq",
                     prune         = TRUE,
                     thinning      = 10,
                     store_theta   = FALSE,
                     store_pi      = FALSE,
                     parallel_init = FALSE)

  # modelling
  covar_formula = stats::as.formula(paste0("~ ",paste(covars,collapse = ' + ')))
  if(!is.null(priors)){
    out <- keyATM(docs              = stData,
                  no_keyword_topics = noMarkerCts,
                  keywords          = markerGs,
                  model             = "covariates",
                  model_settings    = list(covariates_data    = covarDat,
                                           covariates_formula = covar_formula),
                  options           = my_options,
                  keep              = c("Z", "S"),
                  priors=priors)
  }else{
    out <- keyATM(docs              = stData,
                  no_keyword_topics = noMarkerCts,
                  keywords          = markerGs,
                  model             = "covariates",
                  model_settings    = list(covariates_data    = covarDat,
                                           covariates_formula =  covar_formula),
                  options           = my_options,
                  keep              = c("Z", "S"))
  }

  # save the results
  resDir=paste0(outDir,'/inst_',seed,'/')
  if(!dir.exists(resDir)){
    dir.create(resDir,recursive = T)
  }

  saveRDS(out,paste0(resDir,'covar_model.rds'))
  return(list(model=out))
}
