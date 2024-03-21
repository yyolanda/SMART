#' a scatterpie plot showing the predicted cell type proportions of all cell types
#'
#' @param stLoc a data.frame with spatial coordinates
#' @param CT_prop the cell type proportions
#' @param cols the color scheme (a character vector of multiple colors)
#'
#' @return a ggplot object
#' @import scatterpie
#' @import ggthemes
#' @import ggplot2
#' @import dplyr
#' @import randomcoloR
#' @export
#'
#' @examples
#' \dontrun{
#' SMART_scatterpie(visium_mouse_brain$loc,res$ct_proportions)
#' }
#'
SMART_scatterpie <- function(stLoc,CT_prop,cols=NULL){
  nColor=ncol(CT_prop)
  val1=sort(unique(stLoc[,'x']))
  val2=dplyr::lag(val1,1)
  rds = mean(val1-val2, na.rm = T)

  CT_prop=as.data.frame(CT_prop) %>% mutate(id=rownames(.))
  stLoc=as.data.frame(stLoc) %>% mutate(id=rownames(.))
  CT_prop=left_join(CT_prop,stLoc)
  if(is.null(cols)){
    cols <- distinctColorPalette(nColor)
  }
  p = ggplot() +
    geom_scatterpie(aes(x=x, y=y, group=id, r=rds*0.95), data=CT_prop,
                    cols=colnames(CT_prop)[1:(ncol(CT_prop)-3)], color=NA) +
    coord_equal()+
    theme_few()+
    xlab('')+
    ylab('') +
    scale_fill_manual(values=cols)+
    guides(fill=guide_legend(nrow=floor((ncol(CT_prop)-3)/5)+1,byrow=TRUE))+
    theme(legend.position='bottom',legend.key=element_rect(fill="#deebf7"),
          axis.text = element_blank(),axis.ticks=element_blank(),
          legend.title=element_blank(),legend.text=element_text(size=9),
          legend.background = element_rect(fill="white",
                                           linetype="solid",
                                           colour ="black"))+
    scale_y_reverse()
  return(p)
}
