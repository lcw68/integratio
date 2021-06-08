
#' Histogram of log ratio in disease and control group with density curves
#'
#' @param i ith genera or species specific taxonomic units which are used for testing
#' @param simu Operational taxonomic unit table
#' @param disease A list of log ratio in disease group
#' @param conrol A list of log ratio in control group
#' @param name Disease Name, showing in the legend
#' @param position Position adjustment, either as a string, or the result of a call to a position adjustment function.
#' @param binwidth The width of the bins. Can be specified as a numeric value or as a function that calculates width from unscaled x. Here, "unscaled x" refers to the original x values in the data, before application of any scale transformation.
#'
#'
#' @return A ggplot showing the distribution of log ratio
#'
#' @examples
#' ecc.ratio =  sapply(1:nrow(otu),function(x){log_ratio(otu[x,which(meta$cariesfree==0),],version=1,epsilon=1)})
#' free.ratio =  sapply(1:nrow(otu),function(x){log_ratio(otu[x,which(meta$cariesfree==1),],version=1,epsilon=1)})
#' hist.comp(3,otu,ecc.ratio,free.ratio)
#'
#' @importFrom ggplot2 ggplot geom_histogram scale_fill_manual geom_density ggtitle
#'
#' @export
#'
hist.comp <- function(i,simu, disease, control, name = "ECC", position = "identity",binwidth = 0.1)
{
    dtsig = data.frame(value = c(disease[[i]],control[[i]]),type = c(rep(name,length(disease[[i]])),rep("control",length(control[[i]]))))
    ggplot2::ggplot(dtsig,aes(x = value,group = type, fill = type))+
      ggplot2::geom_histogram(position = position,binwidth = binwidth,aes(y = ..count..),alpha = 0.5)+ggplot2::geom_density(aes(y=..count..*0.1,col=type),alpha= 0,size = 1.2)+
      ggplot2::scale_fill_manual(values = c("#F8766D","#00BFC4"))+labs(x="")+
      ggplot2::ggtitle(paste0(rownames(simu)[i]))
}


#' Density curves of log ratio for the whole otu
#'
#' @param simu Operational taxonomic unit table
#' @param id numeric vector, the row number of the genera/species we are interested
#'
#'
#' @return A ggplot showing the distribution of log ratio
#'
#' @examples
#' plot.total(otu,1:109)
#'
#' @importFrom ggplot2 ggplot geom_density scale_fill_manual
#'
#' @export
#'
plot.total <- function(simu,id)
{
  rb <- NULL
  for(i in id)
  {
    rb <- rbind(rb,ratio_factor(simu[i,,]))
  }
  q = ggplot2::ggplot(rb,aes(x = ratio,group = ratio_class,fill=ratio_class))+
    ggplot2::geom_density(aes(y=..count..),alpha = 0.3)+
    ggplot2::scale_fill_manual(values = c("DNA=0 & RNA=0" = "#F8766D","DNA>0 & RNA=0" = "#7CAE00","DNA>0 & RNA>0" = "#00BFC4","DNA=0 & RNA>0" = "#C77CFF"))
  q
}

#' Scatterplot of log DNA versus log RNA
#'
#' @param simu Operational taxonomic unit table
#' @param id numeric vector, the row number of the genera/species we are interested
#'
#' @return A ggplot showing the scatterplot of log DNA versus log RNA
#'
#' @examples
#' scatter.total(otu,1:109)
#'
#' @importFrom ggplot2 ggplot geom_point labs
#'
#' @export
#'

scatter.total <- function(simu,id)
{
  rb <- NULL
  for(i in id)
  {
    v = data.frame(logd = log(1+simu[i,,1]),logr = log(1+simu[i,,2]))
    rb <- rbind(rb,v)
  }

  q=ggplot2::ggplot(rb,aes(x=logd,y=logr))+ggplot2::geom_point()+ggplot2::labs(x="log(1+DNA)",y="log(1+RNA)")
  q
}



