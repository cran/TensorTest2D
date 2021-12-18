#' Fitting Second-order Tensor Generalized Regression
#'
#' \kbd{tensorReg2D} is used to fit second-order tensor generalized regression model. It mainly
#' focus on parameter estimation, including parameter coefficients and standard
#' deviation. The function is built upon \strong{Alternating Least Square
#' Algorithm}, so we provide two criterion to determine optimal result (see more
#' details below in Arguments). Also, we offer model complexity
#' measurement,including AIC and BIC.
#'
#' @details \kbd{tensorReg2D} focuses on second-order tensor generalized regression problems.
#' To be more specific, it provides statistical inference for input variables.
#' Moreover, the function isn't restricted to second-order tensor input \kbd{X};
#' it could combine with other meaningful numerical variables \kbd{W}.
#'
#' Since \kbd{tensorReg2D} is based on \strong{Alternating Least Square
#' Algorithm}, we need to pre-define following arguments to meet favorable
#' optimization result.
#'
#' \kbd{n_R}: In the case of regression with the order 2, P-by-G-by-n tensor, we
#' can break a unknown parameter matrix \strong{B}(P-by-G) into multiplication
#' of two matrix \strong{B_1}(P-by-R) and \strong{t(B_2)} (R-by-G), which means
#' that we can estimate the original matrix \strong{B} by iteratively updating
#' \strong{B_1} and \strong{B_2}. In this scenario, \kbd{n_R} equals to the rank
#' of these two approximate matrix \strong{B_1} and \strong{B_2}. Conceivably,
#' \kbd{1 <= n_R <= min(P,G)}, and by properly pre-appointing \kbd{n_R}, we can
#' estimate a unknown parameter matrix. By default, \kbd{n_R = 1}.
#'
#' \kbd{opt}: In optimization algorithm, we have to determine stopping
#' criterion. In \kbd{tensorReg2D}, we offer two criteria. \kbd{If opt = 1}, the
#' criterion is that we stop our execution when the maximum difference between
#' the elements among an estimated parameter matrix \strong{B} with an estimated
#' parameter vector \strong{b} and preceding ones is less than predefined
#' tolerance (\kbd{tol}) . \kbd{If opt = 2}, the criterion is that we stop our
#' execution when the maximum difference between the elements among an estimated
#' approximate parameter matrix \strong{B_1} , \strong{B_2} with an estimated
#' parameter vector \strong{b} and preceding ones is less than predefined
#' tolerance (\kbd{tol}).
#'
#' \kbd{family}: In \kbd{tensorReg2D}, we provide three options for specific generalized regression
#' problem. First, \kbd{family = "gaussian"} using \kbd{identity} link function corresponds to linear regression
#' model, where dependent variable is real number. Next, \kbd{family = "binomial"} based on \kbd{logit} link function
#' corresponds to logistic regression, where dependent variable is restricted to zero or one binary variable. Finally,
#' \kbd{family = "poisson"} built upon \kbd{log} link function corresponds to poisson regression, where dependent
#' variable is non-negative integer.
#'
#' \kbd{max_ite}: In optimization algorithm, we have to beforehand determine maximum iteration beforehand.
#' By default, \kbd{max_ite = 100}.
#'
#' \kbd{tol}: In optimization algorithm, we have to beforehand determine maximum tolerance to cooperate with
#' stopping criterion(\kbd{opt}).
#'
#' @importFrom stats coefficients glm pnorm pt rnorm symnum deviance
#'
#' @param y A numerical vector. Dependent variable.
#' @param X A numerical 3-D array Independent variable(3-D tensor).
#' @param W A numerical matrix. Independent variable.
#' @param n_R A numerical constant. A predefined value determines the rank of
#'   the approximate matrix
#' @param family Family of \kbd{generalized linear model}. Provide three options for model.(see more details in
#'   \strong{Details})
#' @param opt Optimization options. Provide two options for optimization
#'   stopping criterion. \strong{opt = 1 or 2}. (see more details in
#'   \strong{Details})
#' @param max_ite Maximum iteration. The value of maximum iterations for the
#'   algorithm.
#' @param tol Tolerance. The value of tolerance with respect to optimization.
#'
#' @return tensorReg2D returns an object of \kbd{"tsglm"}.
#'
#'  The function, \code{\link{summary.tsglm}} a customized method from generic function
#'  \code{\link[base]{summary}}, can be used to obtain and print a summary and analysis
#'  of variance table of the results.
#'
#'  An object of class \kbd{tsglm} is a list containing at least the following components:
#'
#'   \kbd{ite}: The number of executed times when stopping the function.
#'
#'   \kbd{b_EST}: The estimated coefficients for numerical variables.
#'
#'   \kbd{b_SD}: The estimated standard deviation for numerical variables.
#'
#'   \kbd{b_PV}: The p-value for numerical variables.
#'
#'   \kbd{B_EST}: The estimated coefficients for 3-D tensor variables.
#'
#'   \kbd{B_SD}: The estimated standard deviation for 3-D tensor variables.
#'
#'   \kbd{B_PV}: The p-value for 3-D tensor variables.
#'
#'   \kbd{Residuals}: The differences between true values and prediction values. Provide for
#'   \kbd{family = "gaussian"}.
#'
#'   \kbd{Dev_res}: Deviance residuals for glm. Provide for model except \kbd{family = "gaussian"}.
#'
#'   \kbd{Dev}: The value of Null deviances and Residual deviance. Provide for model except
#'   \kbd{family = "gaussian"}.
#'
#'   \kbd{IC}: The value of AIC and BIC.
#'
#'   \kbd{DoF}: Degree of freedom.
#'
#'   \kbd{call}: The formula of fitted model.
#'
#'   \kbd{family}: The family for model.
#'
#' @examples
#' # Simulation data
#' n <- 500 # number of observations
#' n_P <- 3; n_G <- 64 # dimension of 3-D tensor variables.
#' n_d <- 1 # number of numerical variable, if n_d == 1,  numerical variable equals to intercept.
#' beta_True <- rep(1, n_d)
#' B_True <- c(1,1,1)%*%t(rnorm(n_G)) + c(0, .5, .5)%*%t(rnorm(n_G))
#' B_True <- B_True / 10
#' W <- matrix(rnorm(n*n_d), n, n_d); W[,1] <- 1
#' X <- array(rnorm(n*n_P*n_G), dim=c(n_P, n_G, n))
#' ## Regression
#' y_R<- as.vector(W%*%beta_True + X%hp%B_True + rnorm(n))
#' DATA_R <- list(y = y_R, X = X, W = W)
#' ## Binomial
#' p_B <- exp(W%*%beta_True + X%hp%B_True); p_B <- p_B/(1+p_B)
#' y_B <- rbinom(n, 1, p_B)
#' DATA_B <- list(y = y_B, W = W, X = X)
#' ## Poisson
#' p_P <- exp(W%*%beta_True + X%hp%B_True)
#' y_P <- rpois(n, p_P)
#' y_P[which(y_P > 170)] <- 170 # If y_P > 170, factorial(y_P) == inf.
#' DATA_P <- list(y = y_P, W = W, X = X)
#'
#' # Execution
#' ## Regression
#' result_R <- tensorReg2D(y = DATA_R$y, X = DATA_R$X, W=NULL, n_R = 1, family = "gaussian",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' ## Visualization
#' image(B_True);image(result_R$B_EST)
#' head(predict(result_R, DATA_R$X))
#'
#' ## Binomial
#' result_B <- tensorReg2D(y = DATA_B$y, X = DATA_B$X, W=NULL, n_R = 1, family = "binomial",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' ## Visualization
#' image(B_True);image(result_B$B_EST)
#' head(predict(result_B, DATA_B$X))
#'
#' ## Poisson
#' result_P <- tensorReg2D(y = DATA_P$y, X = DATA_P$X, W=NULL, n_R = 1, family = "poisson",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' ## Visualization
#' image(B_True);image(result_P$B_EST)
#' head(predict(result_P, DATA_P$X))
#'
#' @references
#'   Mengyun Wu, Jian Huang, and Shuangge Ma (2017). Identifying gene-gene
#'   interactions using penalized tensor regression.
#' @references Sheng-Mao Chang, Meng Yang, Wenbin Lu, Yu-Jyun Huang, Yueyang Huang, Hung Hung,
#' Jeffrey C Miecznikowski, Tzu-Pin Lu, Jung-Ying Tzeng,
#' Gene-set integrative analysis of multi-omics data using tensor-based association test,
#' Bioinformatics, 2021;, btab125,
#' (\href{https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab125/6154849}{Link}))
#'
#' @author Sheng-Mao Chang
#'
#' @export
tensorReg2D <- function (y, X, W = NULL, n_R, family, opt = 1, max_ite = 100,
                         tol = 10^(-7))
{
  Tidy_data <- Check_tidy_input(y, X, W, n_R, family)
  result <- ALS(Tidy_data$DATA, n_R, family, opt = opt, max_ite = max_ite,
                tol = tol)
  result$call <- Tidy_data$fm
  result$family <- family
  class(result) <- "tsglm"
  return(result)
}

#' Computation of the vector of Hadamard product values of the matrices
#' \code{X[,,i]} and \code{B}.
#'
#' @param X A numerical 3D array. Each slice is a matrix of size the same as \code{B}.
#' @param B A numerical matrix.
#'
#' @return numerical vector.
#'
#' @author Sheng-Mao Chang
#'
#' @export
`%hp%` <- function(X, B) {
  return(
    sapply(1:dim(X)[3], function(i) sum(X[,,i] * B))
  )
}

#' Computation of two matrices: the column of the output matrix is the
#' Kronecker product of two columns of each input matrix.
#'
#' @param A A numerical matrix.
#' @param B A numerical matrix.
#'
#' @return numerical matrix.
`%b%` <- function(A, B) {
  n_A <- ncol(A)
  n_B <- ncol(B)
  O <- matrix(0, nrow(A)*nrow(B), n_A*n_B)
  # %x%: Kronecker product
  for (i in 1:n_B) {
    for (j in 1:n_A) {
      O[,(i-1)*n_A + j] <- A[,j] %x% B[,i]
    }
  }
  return(O)
}

#' Computation of the matrix with rows being the linearized matrix products of
#' \code{X[,,i]} and \code{B}.
#'
#' @param X A numerical 3D array. Each slice is a matrix of column size
#' the same as the row size of \code{B}.
#' @param B A numerical matrix.
#'
#' @return numerical 3D array.
`%w%` <- function(X, B) {
  return(
    sapply(1:dim(X)[3], function(i) X[,,i] %*% B)
  )
}

#' Computation of the matrix with rows being the linearized matrix products of
#' transposed \code{X[,,i]} and \code{B}.
#'
#' @param X A numerical 3D array. Each slice is a matrix of row size
#' the same as the row size of \code{B}.
#' @param B A numerical matrix.
#'
#' @return numerical vector.
`%wt%` <- function(X, B) {
  return(
    sapply(1:dim(X)[3], function(i) t(X[,,i]) %*% B)
  )
}

#' getGLMCoef: Computing the regression coefficients of generalized linear model.
#'
#' @param y A numerical vector. Dependent variable.
#' @param X A numerical matrix. Independent variable.
#' @param family Family of \kbd{generalized linear model}. Provide three options for model.(see more details in
#'   \strong{Details})
#' @param offset_vec A numerical vector. Prior knowledge to be included in the predictors.
#'
#' @return list.
getGLMCoef <- function(X, y, family, offset_vec) {
  if (family == "gaussian") {
    return(
      solve(t(X) %*% X, t(X) %*% (y - offset_vec))
    )
  } else if (family == "binomial") {
    return(
      as.vector(
        coefficients(glm(y ~ -1 + X, offset = offset_vec, family = "binomial"))
      )
    )
  } else if (family == "poisson") {
    return(
      as.vector(
        coefficients(glm(y ~ -1 + X, offset = offset_vec, family = "poisson"))
      )
    )
  }
}

#' The function confirming the input variables fit the
#' requirements of the tensorReg2D function.
#'
#' @param y A numerical vector. Dependent variable.
#' @param X A numerical 3-D array. Independent variable(3-D tensor).
#' @param W A numerical matrix. Independent variable.
#' @param n_R A numerical constant. A predefined value determines the rank of
#'   the approximate matrix
#' @param family Family of \kbd{generalized linear model}. Provide three options for model.(see more details in
#'   \strong{Details})
#'
#' @return A list.
#'   \kbd{DATA}: list of input data.\code{DATA$y} is the dependent variable.
#'     \code{DATA$X} is the 3-D tensor independent variables.
#'     \code{DATA$W} is other independent variables.
#'
#'   \kbd{fm}: The model expression shown in \code{summary.tsglm}.
Check_tidy_input <- function(y, X, W, n_R, family) {
  Default_family <- c("gaussian", "binomial", "poisson")
  if (!(family %in% Default_family)) {
    stop("Mismatch with default family!")
  }
  else {
    if (!is.vector(y)) {
      stop("The response variable y is not a vector!\n")
    }
    n <- length(y)
    if (family == "binomial") {
      if (sum(!(y %in% c(0, 1))) != 0) {
        stop("The response variable y is not a binary variable!\n")
      }
    }
    if (family == "poisson") {
      if (!all(y >= 0)) {
        stop("The response variable y is not positive!\n")
      }
      if (!all((y - round(y)) == 0)) {
        stop("The response variable y is not a integer variable!\n")
      }
    }
    if (length(dim(X)) != 3) {
      stop("The matrix independent variable X is not a 3-D array!\n")
    }
    n_vec <- dim(X)
    if (n_vec[3] != n) {
      stop("y and X are not conformable!\n")
    }
    if (is.null(W)) {
      fm <- "y ~ I + X"
      W <- matrix(1, n, 1)
    }
    else {
      if (all(W == 1)) {
        fm <- "y ~ I + X"
      }
      else {
        fm <- "y ~ W + X"
      }
    }
    if (length(dim(W)) != 2) {
      stop("W has wrong dimension!\n")
    }
    if (nrow(W) != n) {
      stop("y and W are not conformable!\n")
    }
    if (n_R > min(n_vec[1:2])) {
      stop("n_R must be less or equal than min(dim(X_matrix))! \n")
    }
  }
  return(list(DATA = list(y = y, X = X, W = W), fm = fm))
}

#' The function computing the covariance matrices of the tensor regression parameters.
#'
#' @param DATA A list. The input data. \code{DATA$y} is the dependent variable.
#'   \code{DATA$X} is the 3-D tensor independent variables.
#'   \code{DATA$W} is other independent variables.
#' @param n_R A numerical constant. A predefined value determines the rank of
#'   the approximate matrix
#' @param B1 A numerical matrix. Parameter matrix B1 of the tensor regression.
#' @param B2 A numerical matrix. Parameter matrix B2 of the tensor regression.
#' @param beta A numerical vector. Parameter vector of the covariates, W.
#' @param family Family of \kbd{generalized linear model}. Provide three options for model.(see more details in
#'   \strong{Details})
#'
#' @return A list.
#'   \kbd{V_B}: A numerical matrix. Covariance matrix of the vectorized tensors' parameters.
#'
#'   \kbd{V_b}: A numerical matrix. Covariance matrix of the covariates' parameters.
VAR_ALS <- function(DATA, n_R, B1, B2, beta, family) {
  y <- DATA$y
  X <- DATA$X
  W <- DATA$W
  n_vec <- dim(X)
  n <- n_vec[3]
  n_P <- n_vec[1]
  n_G <- n_vec[2]
  n_d <- ncol(W)
  B <- B1 %*% t(B2)
  C <- matrix(B1[-(1:n_R), ], n_P - n_R, n_R)
  H <- matrix(B1[1:n_R, 1:n_R], n_R, n_R)
  df <- (n_P - n_R + n_G) * n_R + n_d
  w_seq <- X %hp% B + W %*% beta
  res <- y - w_seq
  if (family == "gaussian")
    var_vec <- rep((n - df)/sum(res^2), n)
  if (family == "binomial") {
    var_vec <- exp(w_seq)
    var_vec <- var_vec/(1 + var_vec)
    var_vec <- var_vec * (1 - var_vec)
  }
  if (family == "poisson")
    var_vec <- exp(w_seq)
  I <- matrix(0, df, df)
  for (i in 1:n) {
    tKi2 <- t(diag(n_R) %x% matrix(X[(n_R + 1):n_P, ,
                                     i], n_P - n_R, n_G))
    tKi1 <- t(diag(n_R) %x% matrix(X[1:n_R, , i], n_R,
                                   n_G))
    V <- c(W[i, ], t(tKi2) %*% matrix(c(B2), n_G * n_R,
                                      1), tKi2 %*% matrix(c(C), (n_P - n_R) * n_R,
                                                          1) + tKi1 %*% matrix(c(H), n_R * n_R, 1))
    I <- I + V %*% t(V) * var_vec[i]
  }
  Iinv <- try(solve(I), TRUE)
  V_beta <- as.matrix(Iinv[1:n_d, 1:n_d])
  Iinv <- Iinv[-c(1:n_d), -c(1:n_d)]
  if (!is.null(attr(Iinv, "class"))) {
    df <- prod(dim(B))
    V_B <- matrix(NA, df, df)
  }
  else {
    A1 <- cbind(matrix(0, n_G * n_R, (n_P - n_R) * n_R),
                diag(n_G * n_R))
    A2 <- cbind(diag(n_P - n_R) %b% B2, C %x% diag(n_G))
    A3 <- C %x% diag(n_G)
    A4 <- H %x% diag(n_G)
    Vff <- A1 %*% Iinv %*% t(A1)
    V11 <- A4 %*% Vff %*% t(A4)
    V12 <- A4 %*% Vff %*% t(A3)
    V22 <- A2 %*% Iinv %*% t(A2)
    V_B <- rbind(cbind(V11, V12), cbind(t(V12), V22))
  }
  return(list(V_B = V_B, V_b = V_beta))
}


#' The function computing information criterion values and deviances.
#'
#' @param y A numerical vector. Dependent variable.
#' @param w_seq A numerical vector. Fitted values of the tensor regression model.
#' @param df integer. Degree of freedom.
#' @param family Family of \kbd{generalized linear model}. Provide three options for model.(see more details in
#'   \strong{Details})
#'
#' @return A list.
#'
#'   \kbd{IC}: The values of Akaike information criterion (AIC) and
#'   Bayesian information criterion (BIC).
#'
#'   \kbd{DoF}: The values of the residual degrees of freedom for the null model and
#'   the residual degrees of freedom.
#'
#'   \kbd{Dev_res}: The deviance values of residuals.
#'
#'   \kbd{Dev}: The deviance values for the null model and the working model.
Calculate_IC_Dev <- function(y, w_seq, df, family) {
  n <- length(w_seq)
  if (family == "gaussian") {
    res <- y - w_seq
    s2_MLE <- mean(res * res)
    AIC <- n * log(s2_MLE) + 2 * df
    BIC <- n * log(s2_MLE) + log(n) * df
    I_C <- c(AIC, BIC)
    names(I_C) <- c("AIC", "BIC")
    Dof <- c(n - 1, n - df)
    names(Dof) <- c("Null", "Model")
    return(list(IC = I_C, DoF = Dof, Residuals = res))
  }
  if (family == "binomial") {
    logLik_vec <- w_seq * y - log(1 + exp(w_seq))
    logLik <- sum(logLik_vec)
    AIC <- -2 * logLik + 2 * df
    BIC <- -2 * logLik + log(n) * df
    I_C <- c(AIC, BIC)
    names(I_C) <- c("AIC", "BIC")
    Null_deviance <- deviance(glm(y ~ 1, family = "binomial"))
    Residual_deviance <- -2 * logLik
    sign <- y
    sign[sign == 0] <- -1
    deviance_residual <- sign * sqrt(logLik_vec * -2)
    Dev <- c(Null_deviance, Residual_deviance)
    names(Dev) <- c("Null_deviance", "Residual_deviance")
    Dof <- c(n - 1, n - df)
    names(Dof) <- c("Null", "Model")
    return(list(IC = I_C, DoF = Dof, Dev_res = deviance_residual,
                Dev = Dev))
  }
  if (family == "poisson") {
    logLik_vec <- y * w_seq - exp(w_seq) - log(factorial(y))
    logLik <- sum(logLik_vec)
    logLik_saturated_vec <- ifelse(y == 0, 0, y * log(y)) -
      y - log(factorial(y))
    logLik_saturated <- sum(logLik_saturated_vec)
    AIC <- -2 * logLik + 2 * df
    BIC <- -2 * logLik + log(n) * df
    I_C <- c(AIC, BIC)
    names(I_C) <- c("AIC", "BIC")
    Null_deviance <- deviance(glm(y ~ 1, family = "poisson"))
    Residual_deviance <- -2 * (logLik - logLik_saturated)
    y_hat <- exp(w_seq)
    sign <- sign(y - y_hat)
    deviance_residual <- sign * sqrt(-2 * (logLik_vec -
                                             logLik_saturated_vec))
    Dev <- c(Null_deviance, Residual_deviance)
    names(Dev) <- c("Null_deviance", "Residual_deviance")
    Dof <- c(n - 1, n - df)
    names(Dof) <- c("Null", "Model")
    return(list(IC = I_C, DoF = Dof, Dev_res = deviance_residual,
                Dev = Dev))
  }
}

#' The function performing the Alternating Least Square (ALS) Algorithm.
#'
#' @param DATA A list. The input data. \code{DATA$y} is the dependent variable.
#'   \code{DATA$X} is the 3-D tensor independent variables.
#'   \code{DATA$W} is other independent variables.
#' @param n_R A numerical constant. A predefined value determines the rank of
#'   the approximate matrix
#' @param family Family of \kbd{generalized linear model}. Provide three options for model.(see more details in
#'   \strong{Details})
#' @param opt Optimization options. Provide two options for optimization
#'   stopping criterion. \strong{opt = 1 or 2}. (see more details in
#'   \strong{Details})
#' @param max_ite Maximum iteration. The value of maximum iterations for the
#'   algorithm.
#' @param tol Tolerance. The value of tolerance with respect to optimization.
#'
#' @return A list. See \code{tensorReg2D}.
ALS <- function(DATA, n_R, family, opt = opt, max_ite = max_ite,
                tol = tol) {

  # Preprocess data
  y <- DATA$y
  X <- DATA$X
  W <- DATA$W
  n_vec <- dim(DATA$X)
  n <- n_vec[3]
  n_P <- n_vec[1]
  n_G <- n_vec[2]
  n_d <- ncol(W)
  Z <- matrix(0, n, n_P * n_G)
  for (i in 1:n) Z[i, ] <- c(X[, , i])
  Q <- cbind(W, Z)
  B_1_temp <- matrix(rnorm(n_P * n_R), n_P, n_R)
  test <- 1
  # Two kinds of iterations
  ## n_R == n_P
  if (n_R == n_P) {
    temp <- summary(glm(y ~ Q, family = family))$coefficients
    sel <- 1:n_d
    Std_Bhat <- temp[, 2]
    Std_b <- as.matrix(Std_Bhat[sel])
    rownames(Std_b) <- NULL
    Std_B <- matrix(Std_Bhat[-sel], n_P, n_G)
    V_b <- Std_b^2
    V_B <- Std_B^2
    beta <- as.matrix(temp[1:n_d, 1])
    rownames(beta) <- NULL
    B <- matrix(temp[-c(1:n_d), 1], n_P, n_G)
    w_seq <- X %hp% B + W %*% beta
    res <- y - w_seq
    IC_Dev <- Calculate_IC_Dev(y, w_seq, n_P * n_G + n_d,
                               family)
    ite_index <- 1
  }
  # n_R != n_P
  else {
    offset_vec <- rep(0, n)
    beta <- getGLMCoef(W, y, family, offset_vec)
    B_1 <- matrix(B_1_temp[, 1:n_R], n_P, n_R)
    H <- matrix(B_1[1:n_R, 1:n_R], n_R, n_R)
    C <- matrix(B_1[-c(1:n_R), 1:n_R], n_P - n_R, n_R)
    B_2 <- matrix(0, n_G, n_R)
    B <- B_1 %*% t(B_2)
    for (ite_index in 1:max_ite) {
      offset_vec <- W %*% beta
      B_2_new <- matrix(getGLMCoef(t(X %wt% B_1), y, family, offset_vec),
                        n_G, n_R)
      offset_vec <- W %*% beta
      B_1_new <- matrix(getGLMCoef(t(X %w% B_2_new), y, family, offset_vec),
                        n_P, n_R)
      G <- matrix(B_1_new[n_R, n_R], n_R, n_R)
      B_new <- B_1_new %*% t(B_2_new)
      offset_vec <- X %hp% B_new
      beta_new <- getGLMCoef(W, y, family, offset_vec)
      if (opt == 1)
        test <- max(max(abs(B_new - B)), max(abs(beta_new - beta)))
      if (opt == 2)
        test <- max(c(abs(B_1 - B_1_new)), c(abs(B_2 - B_2_new)), abs(beta_new - beta))
      B_1 <- B_1_new
      B_2 <- B_2_new
      beta <- beta_new
      B <- B_new
      if (test < tol)
        break
    }
    beta <- as.matrix(beta)
    w_seq <- X %hp% B + W %*% beta
    df <- n_d + (n_G + n_P - n_R) * n_R
    V <- VAR_ALS(DATA, n_R, B_1, B_2, beta, family)
    if (is.na(V$V_B[1, 1])) {
      df <- prod(dim(B))
      Std_B <- matrix(NA, df, df)
      Std_b <- matrix(NA, n_d, n_d)
    }
    else {
      V_B <- V$V_B
      V_b <- V$Vb
      Std_B <- sqrt(t(matrix(diag(V$V_B), n_G, n_P)))
      Std_b <- sqrt(as.matrix(diag(V$V_b)))
    }
    IC_Dev <- Calculate_IC_Dev(y, w_seq, df, family)
  }
  # Calculate P value
  if (is.na(Std_B[1])) {
    B_PV <- Std_B
    b_PV <- Std_b
  }
  else {
    if(family == "gaussian"){
      B_PV <- pt(-abs(B/Std_B), IC_Dev$DoF[2]) * 2
      b_PV <- pt(-abs(beta/Std_b), IC_Dev$DoF[2]) * 2
    }
    else {
      B_PV <- pnorm(-abs(B/Std_B)) * 2
      b_PV <- pnorm(-abs(beta/Std_b)) * 2
    }
  }
  # Organize results
  X_names <- dimnames(X)[1:2]
  dimnames(B) <- X_names
  dimnames(Std_B) <- X_names
  dimnames(B_PV) <- X_names
  W_names <- colnames(W)
  rownames(beta) <- W_names
  rownames(Std_b) <- W_names
  rownames(b_PV) <- W_names
  if (family == "gaussian") {
    result <- list(ite = ite_index, b_EST = beta, b_SD = Std_b,
                   b_PV = b_PV, B_EST = B, B_SD = Std_B,
                   B_PV = B_PV, Residuals = IC_Dev$Residuals,
                   IC = IC_Dev$IC, DoF = IC_Dev$DoF)
  }
  else {
    result <- list(ite = ite_index, b_EST = beta, b_SD = Std_b,
                   b_PV = b_PV, B_EST = B, B_SD = Std_B,
                   B_PV = B_PV, Dev_res = IC_Dev$Dev_res,
                   Dev = IC_Dev$Dev, IC = IC_Dev$IC, DoF = IC_Dev$DoF)
  }
  return(result)
}
