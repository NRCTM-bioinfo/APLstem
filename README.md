
# APL stemness

This is the source code for APL deconvolution analysis in APL stemness project.

This script needs to be run in the R environment and is applicable to various platforms.

## 1 Installation

### GitHub

To download code from GitHub: `git clone `

Librarying the source code:

``` {r eval = FALSE}

source("APL_deconvolution.R")

```

## 2 Running demos

``` {r eval = FALSE}

## Read expression data
exp.tpm <- readRDS("demos/test_sample.deconvolution.TPM.rds")

## Read APL cell signature 
apl.signature <- read.table("demos/CIBERSORTx_sig_1w_branch_mix_out.txt", header = T, row.names = 1)

exp.tpm %>% dim()

exp.tpm %>% head()

decon.bulk <- deconvolution(exp.tpm, apl.signature, 2)

decon.bulk


# Read model
linear.model <- readRDS("demos/PredModel.rds")

branch.name <- c("Stem_like","CMP_like","S100_GMP_like","GMP_like","Cycling_GMP_like","MDP_like","Bcell","Tcell","Ery_HBB")
pred.score <- NULL
for (i in 1:length(branch.name)) {
  exp = as.matrix(decon.bulk[, i])
  
  lm.model <- linear.model[[i]]
  pred.obv <-  as.numeric(lm.model$coefficients[1]) + 
    as.numeric(lm.model$coefficients[2]) * as.numeric(decon.bulk[, i])
  pred.score <- cbind(pred.score, pred.obv)
}
rownames(pred.score) <- rownames(decon.bulk)
colnames(pred.score) <- branch.name

# Coefficeint matrix
pred.score


pred.percent <- pred.score
for (i in 1:nrow(pred.percent)) {
  if (sum(pred.percent[i, ] < 0) > 0) {
    pred.percent[i, ] <- pred.percent[i, ] - min(pred.percent[i, ]) 
    pred.percent[i, ] <- pred.percent[i, ] / sum(pred.percent[i, ]) * 100
  } else {
    pred.percent[i, ] <- pred.percent[i, ] / sum(pred.percent[i, ]) * 100
  }
}

# Proportion estimation
pred.percent



```

## 3 Running time

For samples less than 50, all code will be executed within 10 minutes. And for 500 samples, the estimated run time is approximately one hour. 

## 4 Current version

Version 1.0.0 


## 5 Reported bugs and solutions

If there is any error in installing the scripts, please contact us via e-mail dyt12423@rjh.com.cn


## 6 Reference

[1] Zeng, A.G.X. et al. A cellular hierarchy framework for understanding heterogeneity and predicting drug response in acute myeloid leukemia. Nat Med (2022).

[2] Tan, Y. et al. A PML/RARα direct target atlas redefines transcriptional deregulation in acute promyelocytic leukemia. Blood 137, 1503-1516 (2021).










