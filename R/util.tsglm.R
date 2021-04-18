#' Predict by Second-order Tensor Generalized Regression 
#'
#' \kbd{predict} method for self-defined class \kbd{"tsglm"}.
#'
#' @importFrom stats coefficients glm pnorm pt rnorm symnum
#'
#' @param object Fitted \kbd{"tsglm"} object.
#' @param newx a 3-dimensional array for \code{x} for which the prediction are interested. 
#' @param type the type of prediction required. The default is \code{type = "link"} that returns 
#' prediction values on the scale of the linear predictors (eta).
#' Alternatively, set \code{type = "response"} for returning predictions on the scale of the response variable. 
#' @param ... further arguments passed to or from other methods.
#'
#' @return There are two types of the output of \kbd{predict.tsglm} function.
#' By setting \code{type = "link"}, it returns the values of the linear predictors; 
#' and by setting \code{type = "response"}, it returns the the expected values of response variable.
#' For example, for a binomial model, the predictions are log-odds (probabilities on logit scale) 
#' if \code{type = "link"}, and \code{type = "response"} gives the predicted probabilities of Y=1.
#'
#' @seealso \code{\link{TRtest.omics}, \link{summary.tsglm}}
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
#' result_R <- TRtest.omics(y = DATA_R$y, X = DATA_R$X, W=NULL, n_R = 1, family = "gaussian",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' ## Visualization
#' image(B_True);image(result_R$B_EST)
#' head(predict(result_R, DATA_R$X))
#'
#' ## Binomial
#' result_B <- TRtest.omics(y = DATA_B$y, X = DATA_B$X, W=NULL, n_R = 1, family = "binomial",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' ## Visualization
#' image(B_True);image(result_B$B_EST)
#' head(predict(result_B, DATA_B$X))
#'
#' ## Poisson
#' result_P <- TRtest.omics(y = DATA_P$y, X = DATA_P$X, W=NULL, n_R = 1, family = "poisson",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' ## Visualization
#' image(B_True);image(result_P$B_EST)
#' head(predict(result_P, DATA_P$X))
#'
#' @author Ping-Yang Chen
#'
#' @export
predict.tsglm <- function(object, newx, type = c("link", "response"), ...){
  
  type <- match.arg(type)
  
  `%i%` <- function(X, B) sapply(1:dim(X)[3], function(i) sum(X[,,i]*B))
  eta <- rep(1, dim(newx)[3]) %*% object$b_EST + newx %i% object$B_EST
  
  if(object$family == "gaussian"){
    return(eta)
  }else if(object$family == "binomial"){
    if(type == "link") {
      return(eta)
    }else if(type == "response") {
      mu <- exp(eta)/(1 + exp(eta))
      return(mu)
    }
  }else if(object$family == "poisson"){
    if(type == "link") {
      return(eta)
    }else if(type == "response") {
      mu <- exp(eta)
      return(mu)
    }
  }
}

#' Plot Effective Image Pixels for A \kbd{"tsglm"} Object
#'
#' \kbd{plot} method for self-defined class \kbd{"tsglm"}.
#'
#' @importFrom stats p.adjust p.adjust.methods
#' @importFrom grDevices rgb
#'
#' @param x an object of class \kbd{"tsglm"}.
#' @param method p-value correction method. See \code{\link[stats]{p.adjust}}.
#' @param background an image data that used as the background of the effectiveness markers. 
#' If \code{background = NULL}, the background color shows the effect size of the each pixel.
#' @param ... further arguments passed to the \code{\link[graphics]{image}} function.
#'
#' @seealso \code{\link{TRtest.omics}, \link{drawpixelmarks}}
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
#'
#' # Binomial Responses
#' p_B <- exp(W%*%beta_True + X%i%B_True); p_B <- p_B/(1+p_B)
#' y_B <- rbinom(n, 1, p_B)
#' DATA_B <- list(y = y_B, W = W, X = X)
#'
#' # Binomial Model
#' result_B <- TRtest.omics(y = DATA_B$y, X = DATA_B$X, W=NULL, n_R = 1, family = "binomial",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' 
#' # Plot the effect size of the effective pixels
#' plot(result_B, "fdr")
#'
#' # Plot the effective pixels with data image as the background
#' x0 <- DATA_B$X[,,which(DATA_B$y == 0)]
#' m0 <- matrix(0, dim(DATA_B$X)[1], dim(DATA_B$X)[2])
#' for (i in 1:dim(x0)[3]) m0 <- m0 + x0[,,i]/dim(x0)[3]
#' plot(result_B, "fdr", m0, col = gray(seq(0, 1, 0.05)))
#'
#' @author Ping-Yang Chen
#'
#' @export
plot.tsglm <- function(x, method = p.adjust.methods, background = NULL, ...){
  
  adjp <- matrix(p.adjust(as.vector(x$B_PV), method = method),
                 nrow(x$B_PV), ncol(x$B_PV))
  marks <- x$B_EST*(adjp < 0.05)
  
  if (is.null(background)) {
    cL <- 20
    colormap <- rgb(c(rep(0.00, 1*cL), seq(0.00, 1.00, length = 4*cL), rep(1, 3), #R
                      rep(1.00, 4*cL), seq(1.00, 0.40, length = 1*cL)),
                    c(rep(0.00, 1*cL), seq(0.00, 1.00, length = 4*cL), rep(1, 3),	#G
                      seq(1.00, 0.00, length = 4*cL), rep(0.00, 1*cL)),
                    c(seq(0.40, 1.00, length = 1*cL), rep(1.00, 4*cL), rep(1, 3), #B
                      seq(1.00, 0.00, length = 4*cL), rep(0.00, 1*cL)),
                    maxColorValue = 1)
    
    cM <- ifelse(max(abs(marks)) == 0, 1, max(abs(marks)))
    cM_digit <- floor(log10(cM))
    if (cM_digit < 0) {
      imgval <- round(marks, -cM_digit + 3)
      cM <- round(cM, -cM_digit + 3)
    } else {
      imgval <- marks
    }
    drawpixelmarks(imgval, marks, grids = TRUE, col = colormap, zlim = c(-1,1)*cM)
    
  } else {
    drawpixelmarks(background, marks, grids = FALSE, ...)
  }
  
}

#' Marking Specific Pixels on the Given Image Plot
#'
#' @import graphics
#'
#' @param img a matrix of image data.
#' @param marks a matrix of the same size as \code{img}. 
#' On the image plot, "red" rectangles are marked on the pixels in which the cells in \code{marks} are positive, 
#' and, "blue" rectangles are marked on the pixels in which the cells in \code{marks} are negative. 
#' For zero-valued cells of \code{marks}, there is no mark on the corresponding pixel for the image plot.
#' @param grids boolean. If \code{grids = TRUE}, grid lines are added for the image plot.
#' @param ... further arguments passed to the \code{\link[graphics]{image}} function.
#'
#' @seealso \code{\link{plot.tsglm}}
#' @examples
#' # 
#'
#' @author Ping-Yang Chen
#'
#' @export
drawpixelmarks <- function(img, marks, grids = FALSE, ...){
  
  stopifnot(all(dim(img) == dim(marks)))
  
  image(1:nrow(img), 1:ncol(img), img, 
        xlab="", ylab="", axes=FALSE, ...)
  if (grids) {
    abline(h = (0:ncol(img)) + .5, col = '#66666666')
    abline(v = (0:nrow(img)) + .5, col = '#66666666')
  }
  for (i in 1:nrow(img)) {
    for (j in 1:ncol(img)) {
      bcol <- ifelse(marks[i,j] == 0, NA, ifelse(marks[i,j] < 0, "blue", "red"))
      rect(i-.5, j-.5, i+.5, j+.5, col = NA, border = bcol, lwd = 1)
    }
  }
  
}


