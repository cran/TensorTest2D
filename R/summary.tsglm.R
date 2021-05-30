#' Summarizing Second-order Tensor Generalized Regression Fits
#'
#' \kbd{summary} method for self-defined class \kbd{"tsglm"}.
#'
#' @details \kbd{summary.tsglm} is combined with \code{\link[base]{print}} to provide
#' formatting the coefficients, standard errors, etc. and additionally gives 'significance stars'
#'
#' @importFrom stats coefficients glm pnorm pt rnorm symnum
#'
#' @param object an object of class \kbd{"tsglm"}.
#' @param ... further arguments passed to or from other methods.

#' @return \kbd{summary.tsglm} returns a printout similar to \code{\link[stats]{summary.glm}}.
#'
#' The printout contains the following components:
#'
#' \kbd{Call}: The formula for fitted model.
#'
#' \kbd{Deviance Residuals}: The summary statistics of deviance residuals. Provide for model
#' except \kbd{family = "gaussian"}.
#'
#' \kbd{Residuals}: The summary statistics of residuals. Provide for \kbd{family = "gaussian"}.
#'
#' \kbd{Coefficients}: The coefficient table includes estimation, standard deviation,
#' test statistics, and p-value of parameters and significance stars.
#'
#' \kbd{Deviance}: The deviance of a fitted model. Provide for model except \kbd{family = "gaussian"}.
#'
#' \kbd{AIC}: Akaike information criterion.
#'
#' @seealso \code{\link[base]{summary}, \link{predict.tsglm}}
#'
#' @examples
#' # Predefined function: sum of hadamard product in each array
#' `%i%` <- function(X, B) sapply(1:dim(X)[3], function(i) sum(X[,,i]*B))
#'
#' # Simulation data
#' n <- 1000 # number of observations
#' n_P <- 3; n_G <- 64 # dimension of 3-D tensor variables.
#' n_d <- 1 # number of numerical variable, if n_d == 1,  numerical variable equals to intercept.
#' beta_True <- rep(1, n_d)
#' B_True <- c(1,1,1)%*%t(rnorm(n_G)) + c(0, .5, .5)%*%t(rnorm(n_G))
#' B_True <- B_True / 10
#' W <- matrix(rnorm(n*n_d), n, n_d); W[,1] <- 1
#' X <- array(rnorm(n*n_P*n_G), dim=c(n_P, n_G, n))
#' ## Regression
#' y_R<- as.vector(W%*%beta_True + X%i%B_True + rnorm(n))
#' DATA_R <- list(y = y_R, X = X, W = W)
#' ## Binomial
#' p_B <- exp(W%*%beta_True + X%i%B_True); p_B <- p_B/(1+p_B)
#' y_B <- rbinom(n, 1, p_B)
#' DATA_B <- list(y = y_B, W = W, X = X)
#' ## Poisson
#' p_P <- exp(W%*%beta_True + X%i%B_True)
#' y_P <- rpois(n, p_P)
#' y_P[which(y_P > 170)] <- 170 # If y_P > 170, factorial(y_P) == inf.
#' DATA_P <- list(y = y_P, W = W, X = X)
#'
#' # Execution
#' ## Regression
#' result_R <- TRtest(y = DATA_R$y, X = DATA_R$X, W=NULL, n_R = 1, family = "gaussian",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' summary(result_R)
#' head(predict(result_R, DATA_R$X))
#'
#' ## Binomial
#' result_B <- TRtest(y = DATA_B$y, X = DATA_B$X, W=NULL, n_R = 1, family = "binomial",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' summary(result_B)
#' head(predict(result_B, DATA_B$X))
#'
#' ## Poisson
#' result_P <- TRtest(y = DATA_P$y, X = DATA_P$X, W=NULL, n_R = 1, family = "poisson",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' summary(result_P)
#' head(predict(result_P, DATA_P$X))
#'
#'
#' @author Mark Chen
#'
#' @export
summary.tsglm <- function(object, ...){
  if(object$family == "gaussian"){
    # Formula
    cat("Call:","\n")
    cat("formula = ",object$call)
    cat("\n\n")
    # Residuals
    cat("Residuals:","\n")
    print(summary(as.vector(object$Residuals)))
    cat("\n")
    # Coefficients
    b <- dim(object$b_EST)[1]; p <- dim(object$B_EST)[1]; g <- dim(object$B_EST)[2]
    if( b == 1){
      names1 <- "(Intercept)"
    }else{
      names1 <- c("(Intercept)", paste0("X", 1:(b-1)))
    }
    names2 <- paste0("X", 1:p, ".", rep(1:g, each = p))
    names <- c(names1, names2)
    mystarformat <- function(x) symnum(x, corr = FALSE, na = FALSE,
                                       cutpoints = c(0, 0.01, 0.05, 0.1, 1),
                                       symbols = c("***", "**", "*", " "))
    df <- object$DoF
    est <- c(object$b_EST, as.vector(object$B_EST))
    se <- c(object$b_SD, as.vector(object$B_SD))
    t_value <- est / se
    pval <- c(object$b_PV, as.vector(object$B_PV))
    sign <- mystarformat(pval)
    coef_matrix <- round(cbind(Estimate = est, `Std. Error` = se, `t value` = t_value), 5)
    coef_df <- as.data.frame(cbind(coef_matrix, `Pr(>|t|)` = format.pval(pval), ` ` = sign))
    rownames(coef_df) <- names
    cat("Coefficients:","\n")
    print(format(coef_df, trim = T))
    cat("---\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1","\n")
  }else{
    # Formula
    cat("Call:","\n")
    cat("formula = ",object$call)
    cat("\n\n")
    # Residuals
    cat("Deviance Residuals:","\n")
    print(summary(as.vector(object$Dev_res)))
    cat("\n")
    # Coefficients
    b <- dim(object$b_EST)[1]; p <- dim(object$B_EST)[1]; g <- dim(object$B_EST)[2]
    if( b == 1){
      names1 <- "(Intercept)"
    }else{
      names1 <- c("(Intercept)", paste0("X", 1:(b-1)))
    }
    names2 <- paste0("X", 1:p, ".", rep(1:g, each = p))
    names <- c(names1, names2)
    mystarformat <- function(x) symnum(x, corr = FALSE, na = FALSE,
                                       cutpoints = c(0, 0.01, 0.05, 0.1, 1),
                                       symbols = c("***", "**", "*", " "))
    df <- object$DoF
    est <- c(object$b_EST, as.vector(object$B_EST))
    se <- c(object$b_SD, as.vector(object$B_SD))
    t_value <- est / se
    pval <- c(object$b_PV, as.vector(object$B_PV))
    sign <- mystarformat(pval)
    coef_matrix <- round(cbind(Estimate = est, `Std. Error` = se, `t value` = t_value), 5)
    coef_df <- as.data.frame(cbind(coef_matrix, `Pr(>|t|)` = format.pval(pval), ` ` = sign))
    rownames(coef_df) <- names
    cat("Coefficients:","\n")
    print(format(coef_df, trim = T))
    cat("---\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1","\n")
    cat("\n")
    cat("    Null deviance:", object$Dev[1], " on", object$DoF[1]," degrees of freedom", "\n")
    cat("Residual deviance:", object$Dev[2], " on", object$DoF[2]," degrees of freedom", "\n")
    cat("AIC:", object$IC[1], "\n\n")
    cat("Number of calculating iterations:", object$ite)
  }
}

