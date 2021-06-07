#' OTU preprocessing function
#'
#'
#' @param data A list of microbiome data containing the information of otu and meta data
#'
#'
#' @return A list of microbiome data after preprocessing
#'
#' @examples
#'
#' data(otu.data)
#' otu <- preprocess_otu(otu.data)
#'
#'
#'
#' @export
#'
preprocess_otu <- function(data){
  otu <- data[['otu']]
  meta <- data[['meta']]
  delete_id <- which(meta$exclusion!='included')
  otu <- otu[,-c(delete_id),]
  meta <- meta[-c(delete_id),]
  otu_dna <- otu[,,1]

  bi_mat <- array(as.numeric(otu_dna!=0), dim(otu_dna))
  delete_f <- which(rowSums(bi_mat)<30)
  otu <- otu[-c(delete_f),,]

  data_list <- list(otu=otu, meta=meta)
  return (data_list)
}



#' Calculate ratio and type
#'
#'
#' @param species_data p*2 matrix, p represents the samples, DNA and RNA columns
#' @param epsilon the constant we add in both denominator and numerator
#'
#' @return A data frame contaning ratio and which type it belongs to for every sample
#'
#' @examples
#'
#' ratio_factor(otu[1,,])
#'
#'
#'
#'
#' @export
#'
ratio_factor <- function(species_data,epsilon = 1e-6){

  ratio <- log10((species_data[,2]+epsilon)/(species_data[,1]+epsilon))
  ratio_class <- rep(0, length(ratio))
  ratio_class[which(species_data[,1]==0 & species_data[,2]==0)]='DNA=0 & RNA=0'
  ratio_class[which(species_data[,1]!=0 & species_data[,2]==0)]='DNA>0 & RNA=0'
  ratio_class[which(species_data[,1]!=0 & species_data[,2]!=0)]='DNA>0 & RNA>0'
  ratio_class[which(species_data[,1]==0 & species_data[,2]!=0)]='DNA=0 & RNA>0'

  #d1_r1 <- as.numeric(ratio_class=='Q3')
  #d1_r0 <- as.numeric(ratio_class=='Q2')
  df <- data.frame(ratio=ratio, ratio_class)
  return (df)
}

#' Two version ratio calculation
#'
#'
#' @param species_data p*2 matrix, p represents the samples, DNA and RNA columns
#' @param version which ratio version should be calculated
#' @param epsilon the constant we add in two versions ratio
#'
#' @return A data frame contaning ratio and which type it belongs to for every sample
#'
#' @examples
#'
#' log_ratio(otu[1,,],version = 1, epsilon = 1)
#'
#'
#'
#'
#' @export
#'
log_ratio <- function(species_data,version = 1, epsilon = 1)
{
  ###remove zero DNA
  species_data <- species_data[which(species_data[,1]!=0), ]
  if(version == 1)
  {
    if(is.null(dim(species_data))) {
      ratio = log10((species_data[2]+epsilon)/(species_data[1]+epsilon))
    }
    else{
      ratio <- log10((species_data[,2]+epsilon)/(species_data[,1]+epsilon))
    }
  }
  else if(version == 2)
  {
    if(is.null(dim(species_data))) {
      ratio = log10(species_data[2]/species_data[1]+epsilon)
    }
    else{
      ratio <- log10(species_data[,2]/species_data[,1]+epsilon)
    }
  }

  return(ratio)
}


#' Construct data-set including meta information and ratio
#'
#'
#' @param i ith genera or species specific taxonomic units we used
#' @param simu Operational taxonomic unit table
#' @param all A list of log ratio for all the species/genes
#' @param df original meta data (including age, caries, race, reads and t3c)
#'
#' @return A data frame contaning ratio and which type it belongs to for every sample
#'
#' @examples
#'
#' ratio.otu = sapply(1:nrow(otu),function(x){log_ratio(otu[x,,],version=1,epsilon=1)})
#' newdf(1,otu,ratio.otu)
#'
#'
#'
#'
#' @export
#'
newdf <- function(i,simu,all,df=df)
{
  df1 <- df[which(simu[i,,1]!=0), ]
  df1 <- data.frame(df1,ratio = all[[i]])
  return(df1)
}




