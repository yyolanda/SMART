#' Spatial transcriptomics deconvolution with the base model of SMART
#'
#' @param stData a spot(rows)-by-gene(columns) matrix.
#' @param markerGs a list with marker genes for each cell type.
#' @param noMarkerCts the number of cell types without marker genes.
#' @param outDir the path to the output directory.
#' @param seed the starting value.
#' @param iterations the number of iterations (default: 2000).
#' @param priors do not change the priors if you don't know what they are.
#'
#' @return a model object
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
#' res <- SMART_base(stData=MPOA$stMat, markerGs=MPOA$markers, noMarkerCts=1,
#'                   outDir='SMART_results', seed=1, iterations=2000, prior=NULL)
#' }
SMART_base <- function(stData, markerGs, noMarkerCts=1,
                       outDir='SMART_results', seed=1,
                       iterations=2000, priors=NULL){
  # read data
  stData <- keyATM_read(texts = as.dfm(stData))

  # options
  my_options <- list(seed          = seed,
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
  if(!is.null(priors)){
    out <- keyATM(docs              = stData,
                  no_keyword_topics = noMarkerCts,
                  keywords          = markerGs,
                  model             = "base",
                  options           = my_options,
                  priors=priors)
  }else{
    out <- keyATM(docs              = stData,
                  no_keyword_topics = noMarkerCts,
                  keywords          = markerGs,
                  model             = "base",
                  options           = my_options)
  }

  # save the results
  resDir=paste0(outDir,'/inst_',seed,'/')
  if(!dir.exists(resDir)){
    dir.create(resDir,recursive = T)
  }

  saveRDS(out,paste0(resDir,'base_model.rds'))
  return(list(model=out, ct_proportions = out$theta, ct_spec_gexp = out$phi))
}
