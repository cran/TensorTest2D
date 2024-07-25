[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/TensorTest2D?color=green)](https://cran.r-project.org/package=TensorTest2D)
[![](http://cranlogs.r-pkg.org/badges/grand-total/badger?color=green)](https://cran.r-project.org/package=badger)
[![Downloads](https://cranlogs.r-pkg.org/badges/TensorTest2D)](https://CRAN.R-project.org/package=TensorTest2D)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/TensorTest2D?color=orange)](https://CRAN.R-project.org/package=TensorTest2D)

# TensorTest2D
An implementation of fitting generalized linear models on second-order tensor type data. The functions within this package mainly focus on parameter estimation, including parameter coefficients and standard deviation.

# Installation
```
git clone https://github.com/yuting1214/TensorTest2D
R CMD INSTALL TensorTest2D
```
or in R console window type the following
```
install.packages("TensorTest2D")
```

# Quick start
```
library(TensorTest2D)

# Simulate data
n <- 500 # number of observations
n_P <- 3; n_G <- 64 # dimension of 3-D tensor variables.
n_d <- 1 # number of numerical variable, if n_d == 1,  numerical variable equals to intercept.
beta_True <- rep(1, n_d)
B_True <- c(1,1,1)%*%t(rnorm(n_G)) + c(0, .5, .5)%*%t(rnorm(n_G))
B_True <- B_True / 10
W <- matrix(rnorm(n*n_d), n, n_d); W[,1] <- 1
X <- array(rnorm(n*n_P*n_G), dim=c(n_P, n_G, n))

## Regression Data
y_R<- as.vector(W%*%beta_True + X%hp%B_True + rnorm(n))
DATA_R <- list(y = y_R, X = X, W = W)

# Execution (Regression)
result_R <- tensorReg2D(y = DATA_R$y, X = DATA_R$X, W=NULL, n_R = 1, family = "gaussian",
                        opt = 1, max_ite = 100, tol = 10^(-7) )
# Visualization
image(B_True);image(result_R$B_EST)
head(predict(result_R, DATA_R$X))
```

# Relevant Packages
* [tensor](https://CRAN.R-project.org/package=tensor): The tensor product of two arrays is notionally an outer product of the arrays collapsed in specific extents by summing along the appropriate diagonals.
* [rTensor](https://CRAN.R-project.org/package=rTensor): Tools for Tensor Analysis and Decomposition
* [tensorregress](https://CRAN.R-project.org/package=tensorregress): Implement the alternating algorithm for supervised tensor decomposition with interactive side information.

# Publications

* Ping-Yang Chen/Hsing-Ming Chang/Yu-Ting Chen/Jung-Ying Tzeng/Sheng-Mao Chang* (2022) ,TensorTest2D: Fitting Generalized Linear Models with Matrix Covariates,The R Journal,14,152-163,SSCI
