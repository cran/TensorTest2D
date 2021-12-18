#' Predict by Second-order Tensor Generalized Regression
#'
#' \kbd{predict} method for self-defined class \kbd{"tsglm"}.
#'
#' @importFrom stats coefficients glm pnorm pt rnorm symnum
#'
#' @param object Fitted \kbd{"tsglm"} object.
#' @param newx a 3-dimensional array for \code{x} for which the prediction are interested.
#' @param neww a numerical matrix for \code{W} for which the prediction are interested.
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
#' @seealso \code{\link{tensorReg2D}, \link{summary.tsglm}}
#'
#' @examples
#' # Predefined function: sum of hadamard product in each array
#' `%hp%` <- function(X, B) sapply(1:dim(X)[3], function(i) sum(X[,,i]*B))
#'
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
#' ## Prediction
#' head(predict(result_R, DATA_R$X))
#'
#' ## Binomial
#' result_B <- tensorReg2D(y = DATA_B$y, X = DATA_B$X, W=NULL, n_R = 1, family = "binomial",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' ## Prediction
#' head(predict(result_B, DATA_B$X))
#'
#' ## Poisson
#' result_P <- tensorReg2D(y = DATA_P$y, X = DATA_P$X, W=NULL, n_R = 1, family = "poisson",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' ## Prediction
#' head(predict(result_P, DATA_P$X))
#'
#' @author Ping-Yang Chen
#'
#' @export
predict.tsglm <- function(object, newx, neww = NULL, type = c("link", "response"), ...){

  type <- match.arg(type)
  Default_type <- c("link", "response")
  if (!(type %in% Default_type)) {
    stop("Mismatch with type!")
  }
  if (length(dim(newx)) != 3) {
    stop("The matrix independent variable X is not a 3-D array!\n")
  }
  if (is.null(neww)) {
    neww <- matrix(1, dim(newx)[3], 1)
  }
  if (nrow(neww) != dim(newx)[3]) {
    stop("neww and newx are not conformable!\n")
  }
  #`%hp%`: see tensorReg2D.R
  eta <- neww %*% object$b_EST + newx %hp% object$B_EST

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
#' Its default value is \code{"none"}.
#' @param alpha double. The value of significance level.
#' Its default value is \code{NULL}. If specify real-valued \code{alpha}, for example 0.05,
#' the image plot marks the pixels with p-values smaller than \code{alpha}.
#' @param type string. The type of values shown on the image pixels when \code{background = NULL}.
#' Set \code{type = 'coef'} for showing the values of estimated coefficients of the \(B\) matrix.
#' Set \code{type = 'tval'} for showing the t-statistics of the coefficients of the \(B\) matrix.
#' If \code{background} is not \code{NULL}, the plot will neglect the choice for
#' \code{type} and show the background image as per user's interest.
#' @param background an image data that used as the background of the effectiveness markers.
#' If \code{background = NULL}, the background color shows the effect size of
#' the each pixel according to the setting in \code{type}.
#' @param showlabels boolean. For \code{showlabels = TRUE}, if the row and column names of the image exist, the row and column names are
#' shown on the sides of the image plot; otherwise, the row and column indices are shown.
#' @param plot.legend boolean. Set \code{plot.legend = TRUE} if the colorbar legend is needed. The default is \code{TRUE}.
#' @param ... further arguments passed to the \code{\link[graphics]{image}} function.
#'
#' @seealso \code{\link{tensorReg2D}, \link{draw.coef}}
#'
#' @examples
#' # Simulation data
#' n <- 500 # number of observations
#' n_P <- 3; n_G <- 16 # dimension of 3-D tensor variables.
#' n_d <- 1 # number of numerical variable, if n_d == 1,  numerical variable equals to intercept.
#' beta_True <- rep(1, n_d)
#' B_True <- c(1, 1, 1)%*%t(rnorm(n_G)) + c(0, .5, .5)%*%t(rnorm(n_G))
#' B_True <- B_True / 10
#' W <- matrix(rnorm(n*n_d), n, n_d); W[,1] <- 1
#' X <- array(rnorm(n*n_P*n_G), dim=c(n_P, n_G, n))
#'
#' # Binomial Responses
#' p_B <- exp(W%*%beta_True + X%hp%B_True); p_B <- p_B/(1+p_B)
#' y_B <- rbinom(n, 1, p_B)
#' DATA_B <- list(y = y_B, W = W, X = X)
#'
#' # Binomial Model
#' result_B <- tensorReg2D(y = DATA_B$y, X = DATA_B$X, W=NULL, n_R = 1,
#' family = "binomial", opt = 1, max_ite = 100, tol = 10^(-7) )
#'
#' # Plot the effect size of the pixels
#' plot(result_B, method = "fdr", alpha = 0.05, type = "coef")
#' # Plot the t-statistics of the coefficients of the pixels
#' plot(result_B, method = "fdr", alpha = 0.05, type = "tval")
#'
#' # Plot the effective pixels with data image as the background
#' x0 <- DATA_B$X[,,which(DATA_B$y == 0)]
#' x1 <- DATA_B$X[,,which(DATA_B$y == 1)]
#' m0 <- m1 <- matrix(0, dim(DATA_B$X)[1], dim(DATA_B$X)[2])
#' for (i in 1:dim(x0)[3]) m0 <- m0 + x0[,,i]/dim(x0)[3]
#' for (i in 1:dim(x1)[3]) m1 <- m1 + x1[,,i]/dim(x1)[3]
#' par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
#' plot(result_B, method = "fdr", alpha = 0.05,
#' background = m0, col = gray(seq(0, 1, 0.05)))
#' title("Category 0")
#' plot(result_B, method = "fdr", alpha = 0.05,
#' background = m1, col = gray(seq(0, 1, 0.05)))
#' title("Category 1")
#'
#' @author Ping-Yang Chen
#'
#' @export
plot.tsglm <- function(x, method = p.adjust.methods, alpha = NULL,
                       type = c("coef", "tval"), background = NULL,
                       showlabels = TRUE, plot.legend = TRUE, ...){
  type <- match.arg(type)
  if (length(method > 1)) { method <- "none" }
  stopifnot(is.null(alpha) | (alpha > 0) | (alpha < 1))
  adjp <- matrix(p.adjust(as.vector(x$B_PV), method = method),
                 nrow(x$B_PV), ncol(x$B_PV))
  if (is.null(alpha)) {
    marks <- matrix(0, nrow(x$B_EST), ncol(x$B_EST))
  } else {
    marks <- x$B_EST*(adjp < alpha)
  }
  if (type == "coef") {
    imgval <- x$B_EST
  } else {
    imgval <- x$B_EST/x$B_SD
  }

  if (is.null(background)) {
    cL <- 20
    cLim <- .3
    #color goes to black at extreme values (clim lower than 1.0, the lower the darker)
    # colormap <- rgb(c(rep(0.00, 1*cL), seq(0.00, 1.00, length = 4*cL), rep(1, 3), #R
    #                   rep(1.00, 4*cL), seq(1.00, cLim, length = 1*cL)),
    #                 c(rep(0.00, 1*cL), seq(0.00, 1.00, length = 4*cL), rep(1, 3), #G
    #                   seq(1.00, 0.00, length = 4*cL), rep(0.00, 1*cL)),
    #                 c(seq(cLim, 1.00, length = 1*cL), rep(1.00, 4*cL), rep(1, 3), #B
    #                   seq(1.00, 0.00, length = 4*cL), rep(0.00, 1*cL)),
    #                 maxColorValue = 1)
    #
    #color goes to brighter red and blue at extreme values (clim lower but not equal to 1.0, the lower the darker)
    colormap <- rgb(c(seq(cLim, 1.00, length = 5*cL), rep(1, 3), rep(1.00, 5*cL)),
                    c(seq(cLim, 1.00, length = 5*cL), rep(1, 3), seq(1.00, cLim, length = 5*cL)),
                    c(rep(1.00, 5*cL),                rep(1, 3), seq(1.00, cLim, length = 5*cL)),
                    maxColorValue = 1)
    #
    cM <- ifelse(max(abs(imgval)) == 0, 1, max(abs(imgval)))
    cM_digit <- floor(log10(cM))
    if (cM_digit < 0) {
      imgval <- round(imgval, -cM_digit + 3)
      cM <- round(cM, -cM_digit + 3)
    }
    draw.coef(imgval, marks, markstyle = 'black', showlabels = showlabels, grids = TRUE, plot.legend = plot.legend, col = colormap, zlim = c(-1,1)*cM)

  } else {

    draw.coef(background, marks, markstyle = 'bi-dir', showlabels = showlabels, grids = FALSE, plot.legend = plot.legend, ...)
  }

}

#' Marking Specific Pixels on the Given Image Plot
#'
#' @import graphics
#'
#' @param img a matrix of image data.
#' @param marks a matrix of the same size as \code{img}.
#' On the image plot, the pixels are marked if the corresponding cells in \code{marks} are non-zero.
#' The user can specify the style of the marks through \code{markstyle}.
#' @param markstyle string. The style of pixels' marks. If \code{markstyle = 'black'}, the rectangles
#' are marked by black edges for non-zero cells in \code{marks}.
#' If \code{markstyle = 'bi-dir'}, "red" rectangles are marked on the pixels in which the cells in \code{marks} are positive,
#' and, "blue" rectangles are marked on the pixels in which the cells in \code{marks} are negative.
#' @param showlabels boolean. For \code{showlabels = TRUE}, if \code{dimnames(img)} exists, the row and column names are
#' shown on the sides of the image plot; otherwise, the row and column indices are shown.
#' @param plot.legend boolean. Set \code{plot.legend = TRUE} if the colorbar legend is needed. The default is \code{TRUE}.
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
draw.coef <- function(img, marks, markstyle = c("black", "bi-dir"),
                      showlabels = TRUE, plot.legend = TRUE, grids = FALSE, ...){
  stopifnot(all(dim(img) == dim(marks)))
  markstyle <- match.arg(markstyle)
  B_names <- dimnames(img)
  if(is.null(B_names)){
    B_names <- list(paste0(1:nrow(img)), paste0(1:ncol(img)))
  }

  if (plot.legend){
    layout(matrix(c(1, 2), 1, 2), widths = c(ncol(img)*.9, ncol(img)*.1))
    par(mar = c(5, 4, 4, .1)+.1, las = 1)
  }
  image(1:nrow(img), 1:ncol(img), img,
        xlab="", ylab="", xlim=c(.49,nrow(img)+.51), ylim=c(.49,ncol(img)+.51), axes=FALSE, ...)
  if (grids) {
    abline(h = (0:ncol(img)) + .5, col = '#66666666')
    abline(v = (0:nrow(img)) + .5, col = '#66666666')
  }
  for (i in 1:nrow(img)) {
    for (j in 1:ncol(img)) {
      if (markstyle == "black") {
        if (marks[i,j] != 0) {
          rect(i-.5, j-.5, i+.5, j+.5, col = NA, border = "black", lwd = 2)
          arrows(x0 = c(i-.5, i-.5), y0 = c(j-.5, j+.5),
                 x1 = c(i+.5, i+.5), y1 = c(j+.5, j-.5), col = "black", code = 0, lwd = 2)
        }
      } else {
        bcol <- ifelse(marks[i,j] == 0, NA, ifelse(marks[i,j] < 0, "blue", "red"))
        rect(i-.5, j-.5, i+.5, j+.5, col = NA, border = bcol, lwd = 2)
      }
    }
  }
  if (showlabels) {
    axis(1, at = 1:nrow(img), labels = B_names[[1]], tick = 0, line = -.7)
    axis(2, at = 1:ncol(img), labels = B_names[[2]], las = 2, tick = 0, line = -.7)
  }
  if (plot.legend){
    par(mar = c(5, 0, 4, 1)+.1)
    clen <- 50
    image(0, 1:clen, matrix(seq(min(img), max(img), length = clen), 1, clen),
          xlab="", ylab="", axes=FALSE, ...)
    mtext(signif(max(img), 2), 3); mtext(signif(min(img), 2), 1)
  }

}
