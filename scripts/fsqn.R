# script from jenniferfranks/FSQN: Feature Specific Quantile Normalization

library(preprocessCore)

quantileNormalizeByFeature <- function(matrix_to_normalize,
                                       target_distribution_matrix){
  
  if (ncol(matrix_to_normalize) != ncol(target_distribution_matrix)){
    cat("ERROR: Data matrices are not compatible - column lengths differ!")
  }
  else{
    
    data.qn <- matrix(0, nrow = nrow(matrix_to_normalize),
                      ncol = ncol(matrix_to_normalize))
    
    for (i in 1:ncol(matrix_to_normalize)){
      feature.to.normalize <- matrix_to_normalize[,i]
      target.feature.dist <- target_distribution_matrix[,i]
      result <- normalize.quantiles.use.target(
        x = as.matrix(feature.to.normalize),
        target = target.feature.dist,
        copy = TRUE)
      data.qn[,i] <- result
    }
    rownames(data.qn) = rownames(matrix_to_normalize)
    colnames(data.qn) = colnames(matrix_to_normalize)
    return(data.qn)
  }
}
