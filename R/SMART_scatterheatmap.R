#' a scatterheatmap showing the predicted cell type proportion of a selected cell type
#'
#' @param stLoc a data.frame with spatial coordinates
#' @param CT_prop the cell type proportions
#' @param celltype the cell type to plot
#' @param cols the color scheme for the heatmap (a character vector of two colors)
#'
#' @return a ggplot object
#' @import scatterpie
#' @import ggthemes
#' @import ggplot2
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' SMART_scatterheatmap(visium_mouse_brain$loc,res$ct_proportions, celltype = 'Oligos')
#' }
SMART_scatterheatmap <- function(stLoc,CT_prop,celltype,cols=c('white','red')){

  CT_prop=as.data.frame(CT_prop) %>% mutate(id=rownames(.))
  stLoc=as.data.frame(stLoc) %>% mutate(id=rownames(.))
  CT_prop=left_join(CT_prop,stLoc)
  p = ggplot(CT_prop,aes(x=x, y=y, fill=get(celltype))) +
      geom_point(shape=21,color='black',size=3) +
      scale_fill_gradient(low=cols[1],high=cols[2])+
      coord_equal()+
      theme_few()+
      xlab('')+
      ylab('') +
      theme(legend.position=c(0.9,0.1),legend.key=element_rect(fill="#deebf7"),
            axis.text = element_blank(),axis.ticks=element_blank(),
            legend.title=element_blank(),legend.text=element_text(size=13),
            legend.background = element_rect(fill="grey",
                                             linetype="solid",
                                             colour ="black"))+
      scale_y_reverse()

  return(p)
}
