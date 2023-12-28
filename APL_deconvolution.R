
library(dplyr)
library(e1071)

# deconvolution analysis of APL bulk RNA-seq data
# exp.mat: TPM matrix of APL bulk RNA-seq data
# sig.mat: Signature matrix of APL branches obtained from single cell RNA-seq
# ncores: cores used for analysis
deconvolution <- function(exp.mat, sig.mat, ncores) {

	message("Version 1.0.1")

  merge.mat <- list(y = as.matrix(exp.mat[rownames(sig.mat), ]), x = as.matrix(sig.mat) )
  nu <- c(0.25, 0.5, 0.75, 1)

  res <- parallel::mclapply(seq_len(ncol(merge.mat$y)), function(z) {
    lapply(nu, function(i) {
      suppressMessages(
        suppressWarnings(
          svm(x = merge.mat$x,
              y = as.vector(merge.mat$y[,z]), 
              nu = i, 
              type = 'nu-regression',
              scale = FALSE, 
              kernel = 'linear')
        ))
    })
  }, mc.cores = ncores)
  
  rmse <-  lapply(res, function(z) { 
    lapply(z, function(i){
      sqrt(sum(i$residuals ^ 2)/length(i$residuals))
    })
  })
  rmse <- unlist(lapply(rmse, function(z) which.min(unlist(z))))
  # extract models with min rmse
  mod <- lapply(seq_along(res), function(z) {res[[z]][[rmse[z]]] })
  
  # return coefficients
  coefs <- lapply(mod, function(z) {as.numeric(t(z$coefs) %*% z$SV)})
  coefs <- as.data.frame(t(data.frame(coefs)))
  rownames(coefs) <- colnames(merge.mat$y)
  colnames(coefs) <- colnames(merge.mat$x)
  
  return(coefs)
}

