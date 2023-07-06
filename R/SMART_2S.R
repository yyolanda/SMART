#' Spatial transcriptomics deconvolution with the base model of SMART
#'
#' @param stData a spot(rows)-by-gene(columns) spatial transcriptomics matrix.
#' @param markerGs_S1 a list with marker genes for major cell types.
#' @param markerGs_S2 a list with marker genes for subtypes of the cell type of interest.
#' @param noMarkerCts_S1 the number of cell types without marker genes in stage one.
#' @param noMarkerCts_S2 the number of cell types without marker genes in stage two.
#' @param CT_OI cell type of interest
#' @param CT_similar a cell type that is transcriptomically similar to the cell type of interest (can be a character string indicating a cell type name, 'auto' for automatic selection or NULL)
#' @param outDir the path to the output directory.
#' @param seed the starting value.
#' @param iterations the number of iterations (default: 2000).
#' @param priors_S1 do not change the priors if you don't know what they are.
#' @param priors_S2 do not change the priors if you don't know what they are.
#'
#' @return a model object
#' @return a model
#' @return a cell type proportion matrix
#' @import keyATM
#' @import quanteda
#' @importFrom vegan vegdist
#' @importFrom usedist dist_get
#' @export
#'
#' @examples
#' \dontrun{
#' library(SMART)
#' data(NSCLC)
#' res <- SMART_2S(NSCLC$stData, NSCLC$markers_S1, NSCLC$markers_S2,
#'                 noMarkerCts_S1=1, noMarkerCts_S2=1,
#'                 CT_OI='DC',CT_similar='macrophage',
#'                 outDir='SMART_results')
#' }
SMART_2S <- function(stData, markerGs_S1, markerGs_S2,
                     noMarkerCts_S1=1, noMarkerCts_S2=1,
                     CT_OI,CT_similar=NULL,
                     outDir='SMART_results', seed=1, iterations=2000,
                     priors_S1=NULL, priors_S2=NULL){

  stopifnot(CT_OI %in% names(markerGs_S1))
  stopifnot(is.null(CT_similar) | CT_similar %in% c(names(markerGs_S1),'auto'))

  # stage one modeling
  SMART_res <- SMART_base(stData=stData, markerGs=markerGs_S1, noMarkerCts=noMarkerCts_S1,
                          outDir=outDir, seed=seed, iterations=iterations, priors=priors_S1)

  # get the most similar cell type
  if(!is.null(CT_similar)){
    if(CT_similar=='auto'){
      ct_dist=vegan::vegdist(SMART_res$model$phi, method="euclidean", upper=T, diag=T)
      ctoi_dist=usedist::dist_get(ct_dist,CT_OI,rownames(SMART_res$model$phi))
      names(ctoi_dist)=rownames(SMART_res$model$phi)
      ctoi_dist=ctoi_dist[!grepl('Other',names(ctoi_dist))]
      cts=names(sort(ctoi_dist)[1:2])
    }else if(CT_similar %in% names(markerGs_S1)){
      cts=c(CT_OI,CT_similar)
    }else{
      print('Please use a valid cell type name for CT_similar.')
    }
  }else{
    cts=CT_OI
  }

  # get counts for the cell type of interest
  doccounts=rowSums(stData) # total counts of each spot
  temp=doccounts*SMART_res$model$theta # counts per cell type per spot
  stopifnot(all(abs(rowSums(temp)-doccounts)<0.1))

  newX=list()
  for(i in cts){
    newX[[i]]=temp[,i] %*% t(SMART_res$model$phi[i,]) # counts by gene for a particular cell type
  }
  newX=Reduce('+', newX)
  scalingF=max(1,floor(1e7/sum(newX)))
  newX=newX * scalingF # scale the counts
  newX=apply(newX,c(1,2),round)
  rownames(newX)=rownames(SMART_res$model$theta)

  # stage two modeling
  SMART_res_sub <- SMART_base(stData=newX, markerGs=markerGs_S2,
                              noMarkerCts=noMarkerCts_S2,
                              outDir=paste0(outDir,'/',CT_OI,'_subtype/'),
                              seed=seed, iterations=iterations, priors=priors_S2)

  # calculate the final cell type proportions
  if(length(cts)!=1){
    theta_new=SMART_res_sub$model$theta * rowSums(SMART_res$model$theta[,cts])
  }else{
    theta_new=SMART_res_sub$model$theta * SMART_res$model$theta[,cts]
  }
  stopifnot(all(rownames(SMART_res$ct_proportions)==rownames(theta_new)))
  props = cbind(SMART_res$ct_proportions[,-match(cts,colnames(SMART_res$ct_proportions))],theta_new)

  return(list(S1=SMART_res, S2=SMART_res_sub,
              final_ct_proportions = props))
}
