#' Read the pre-processed MNIST dataset
#'
#' A pre-processed MNIST dataset with each image of 10*10 = 100 pixels.
#'
#' @details The original MNIST handwritten digit is the image of 28*28 = 784 pixels.
#' The pre-processing procedure is as follows.
#' First, the original 28 by 28 image is separated into 14 by 14 clusters, 
#' each one is a 2 by 2 block. For each cluster, the maximal value in its cells 
#' is taken, and therefore, the original MNIST image is reduced into 14 by 14 image. 
#' Because the surrounding cells of the reduced image are usually zero value, 
#' we only take the center 10 by 10 sub-image by cutting the edge cells. 
#' Thus, the pre-processed MNIST dataset has images of 10*10 = 100 pixels.
#'
#' @docType data
#'
#' @usage data(mnist_mp2c2)
#'
#' @format A list of two sublists, \code{mnist_mp2c2$train} and \code{mnist_mp2c2$test}.
#' In each sublist, the data is stored as a list of length two, image and label. 
#' The image is a 3-dimensional array of size (10, 10, n), where n represents the data size.
#' For i = 1, ..., n the i-th slice of image is an integer matrix with elements in [0, 255]
#' representing the image of 10*10 = 100 pixels in grey scale.
#' The label is a vector of length n. The i-th value is the digit of the i-th slice of image.
#'
#' @keywords datasets
#'
#' @references LeCun, Y., Bottou, L., Bengio, Y., & Haffner, P. (1998). Gradient-based learning applied to document recognition. Proceedings of the IEEE, 86(11), 2278-2324.
#' (\href{https://ieeexplore.ieee.org/abstract/document/726791}{URL})
#' 
#' @references Rafael A. Irizarry and Amy Gill (2019). dslabs: Data Science Labs. R package version 0.7.3.
#' (\href{https://CRAN.R-project.org/package=dslabs}{URL})
#' 
#' @references LeCun, Y. http://yann.lecun.com/exdb/mnist/
#'
#' @source \href{https://CRAN.R-project.org/package=dslabs}{https://CRAN.R-project.org/package=dslabs}
#'
#' @examples
#' data(mnist_mp2c2)
#' dim(mnist_mp2c2$train$image)
#' # 10    10 60000
#' image(mnist_mp2c2$train$image[,,1])
#' mnist_mp2c2$train$label[1]
"mnist_mp2c2"
