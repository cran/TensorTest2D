context("test - TensorTest2D")

set.seed(5)
# Simulation data
n <- 500 # number of observations
n_P <- 3; n_G <- 64 # dimension of 3-D tensor variables.
n_d <- 2 # number of numerical variable
beta_True <- rep(1, n_d)
B_True <- c(1, 1, 1) %*% t(rnorm(n_G)) + c(0, .5, .5) %*% t(rnorm(n_G))
B_True <- B_True / 10
W <- matrix(rnorm(n*n_d), n, n_d); W[,1] <- 1
X <- array(rnorm(n*n_P*n_G), dim = c(n_P, n_G, n))
## Regression
y_R <- as.vector(W %*% beta_True + X %hp% B_True + rnorm(n))
DATA_R <- list(y = y_R, X = X, W = W)
## Binomial
p_B <- exp(W %*% beta_True + X %hp% B_True); p_B <- p_B/(1 + p_B)
y_B <- rbinom(n, 1, p_B)
DATA_B <- list(y = y_B, W = W, X = X)
## Poisson
p_P <- exp(W %*% beta_True + X %hp% B_True)
y_P <- rpois(n, p_P)
y_P[which(y_P > 170)] <- 170
DATA_P <- list(y = y_P, W = W, X = X)

# Execution
## Regression
result_R <- tensorReg2D(y = DATA_R$y, X = DATA_R$X, W = DATA_R$W,
                        n_R = 1, family = "gaussian")
## Binomial
result_B <- tensorReg2D(y = DATA_B$y, X = DATA_B$X, W = DATA_B$W,
                        n_R = 1, family = "binomial")
## Poisson
result_P <- tensorReg2D(y = DATA_P$y, X = DATA_P$X, W = DATA_P$W,
                        n_R = 1, family = "poisson")
## Regression
result_RI <- tensorReg2D(y = DATA_R$y, X = DATA_R$X, W = NULL,
                         n_R = 1, family = "gaussian")
## Binomial
result_BI <- tensorReg2D(y = DATA_B$y, X = DATA_B$X, W = NULL,
                         n_R = 1, family = "binomial")
## Poisson
result_PI <- tensorReg2D(y = DATA_P$y, X = DATA_P$X, W = NULL,
                         n_R = 1, family = "poisson")

data(omics)
data(mnist_mp2c2)

## Begin Tests
test_that("tensorReg2D", {

  expect_error(tensorReg2D(y = DATA_R$y, X = DATA_R$X, W = DATA_R$W, n_R = 1, family = "somethingelse"))
  expect_error(tensorReg2D(y = DATA_R$y, X = DATA_B$X, W = DATA_B$W, n_R = 1, family = "binomial"))
  expect_error(tensorReg2D(y = DATA_R$y, X = DATA_P$X, W = DATA_P$W, n_R = 1, family = "poisson"))
  expect_error(tensorReg2D(y = DATA_R$y, X = DATA_R$X, W = rep(1, n), n_R = 1, family = "gaussian"))
  expect_error(tensorReg2D(y = DATA_R$y, X = DATA_R$X, W = DATA_R$W, n_R = max(c(n_P, n_G)) + 1, family = "gaussian"))

  all.equal(dim(result_R$b_EST), c(2L, 1)); all.equal(dim(result_R$B_EST), c(n_P, n_G))
  all.equal(dim(result_R$b_SD),  c(2L, 1)); all.equal(dim(result_R$B_SD),  c(n_P, n_G))
  all.equal(dim(result_R$b_PV),  c(2L, 1)); all.equal(dim(result_R$B_PV),  c(n_P, n_G))
  all.equal(length(result_R$Residuals), n)
  all.equal(length(result_R$IC), 2L);       all.equal(length(result_R$DoF), 2L)
  all.equal(result_R$call,   "y ~ W + X");  all.equal(result_R$family, "gaussian")

  all.equal(dim(result_B$b_EST), c(2L, 1)); all.equal(dim(result_B$B_EST), c(n_P, n_G))
  all.equal(dim(result_B$b_SD),  c(2L, 1)); all.equal(dim(result_B$B_SD),  c(n_P, n_G))
  all.equal(dim(result_B$b_PV),  c(2L, 1)); all.equal(dim(result_B$B_PV),  c(n_P, n_G))
  all.equal(length(result_B$Dev_res), n);   all.equal(length(result_B$Dev), 2L)
  all.equal(length(result_B$IC), 2L);       all.equal(length(result_B$DoF), 2L)
  all.equal(result_B$call, "y ~ W + X");    all.equal(result_B$family, "binomial")

  all.equal(dim(result_P$b_EST), c(2L, 1)); all.equal(dim(result_P$B_EST), c(n_P, n_G))
  all.equal(dim(result_P$b_SD),  c(2L, 1)); all.equal(dim(result_P$B_SD),  c(n_P, n_G))
  all.equal(dim(result_P$b_PV),  c(2L, 1)); all.equal(dim(result_P$B_PV),  c(n_P, n_G))
  all.equal(length(result_P$Dev_res), n);   all.equal(length(result_P$Dev), 2L)
  all.equal(length(result_P$IC), 2L);       all.equal(length(result_P$DoF), 2L)
  all.equal(result_P$call, "y ~ W + X");    all.equal(result_P$family, "poisson")

  expect_error(predict(result_R, DATA_R$X, DATA_R$W, type = "somethingelse"))
  expect_error(predict(result_R, DATA_R$X, type = "link"))
  expect_error(predict(result_B, DATA_R$X, type = "link"))
  expect_error(predict(result_P, DATA_R$X, type = "link"))
  expect_error(predict(result_R, DATA_R$X, type = "response"))
  expect_error(predict(result_B, DATA_R$X, type = "response"))
  expect_error(predict(result_P, DATA_R$X, type = "response"))

  all.equal(length(predict(result_R, DATA_R$X, DATA_R$W, type = "link")), n)
  all.equal(length(predict(result_B, DATA_R$X, DATA_R$W, type = "link")), n)
  all.equal(length(predict(result_P, DATA_R$X, DATA_R$W, type = "link")), n)
  all.equal(length(predict(result_R, DATA_R$X, DATA_R$W, type = "response")), n)
  all.equal(length(predict(result_B, DATA_R$X, DATA_R$W, type = "response")), n)
  all.equal(length(predict(result_P, DATA_R$X, DATA_R$W, type = "response")), n)

  # Intercept only
  all.equal(dim(result_RI$b_EST), c(2L, 1)); all.equal(dim(result_RI$B_EST), c(n_P, n_G))
  all.equal(dim(result_RI$b_SD),  c(2L, 1)); all.equal(dim(result_RI$B_SD),  c(n_P, n_G))
  all.equal(dim(result_RI$b_PV),  c(2L, 1)); all.equal(dim(result_RI$B_PV),  c(n_P, n_G))
  all.equal(length(result_RI$Residuals), n)
  all.equal(length(result_RI$IC), 2L);       all.equal(length(result_RI$DoF), 2L)
  all.equal(result_RI$call,   "y ~ W + X");  all.equal(result_RI$family, "gaussian")

  all.equal(dim(result_BI$b_EST), c(2L, 1)); all.equal(dim(result_BI$B_EST), c(n_P, n_G))
  all.equal(dim(result_BI$b_SD),  c(2L, 1)); all.equal(dim(result_BI$B_SD),  c(n_P, n_G))
  all.equal(dim(result_BI$b_PV),  c(2L, 1)); all.equal(dim(result_BI$B_PV),  c(n_P, n_G))
  all.equal(length(result_BI$Dev_res), n);   all.equal(length(result_BI$Dev), 2L)
  all.equal(length(result_BI$IC), 2L);       all.equal(length(result_BI$DoF), 2L)
  all.equal(result_BI$call, "y ~ W + X");    all.equal(result_BI$family, "binomial")

  all.equal(dim(result_PI$b_EST), c(2L, 1L)); all.equal(dim(result_PI$B_EST), c(n_P, n_G))
  all.equal(dim(result_PI$b_SD),  c(2L, 1L)); all.equal(dim(result_PI$B_SD),  c(n_P, n_G))
  all.equal(dim(result_PI$b_PV),  c(2L, 1L)); all.equal(dim(result_PI$B_PV),  c(n_P, n_G))
  all.equal(length(result_PI$Dev_res), n);    all.equal(length(result_PI$Dev), 2L)
  all.equal(length(result_PI$IC), 2L);        all.equal(length(result_PI$DoF), 2L)
  all.equal(result_PI$call, "y ~ W + X");     all.equal(result_PI$family, "poisson")

  all.equal(length(predict(result_RI, DATA_R$X, type = "link")), n)
  all.equal(length(predict(result_BI, DATA_R$X, type = "link")), n)
  all.equal(length(predict(result_PI, DATA_R$X, type = "link")), n)
  all.equal(length(predict(result_RI, DATA_R$X, type = "response")), n)
  all.equal(length(predict(result_BI, DATA_R$X, type = "response")), n)
  all.equal(length(predict(result_PI, DATA_R$X, type = "response")), n)

  all.equal(names(omics), c("Y", "omics"))
  all.equal(names(mnist_mp2c2$train), c("image", "label"))
  all.equal(names(mnist_mp2c2$test), c("image", "label"))
  all.equal(dim(omics$omics), c(3L, 10L, 68L))
  all.equal(dim(mnist_mp2c2$train$image), c(10L, 10L, 60000L))
  all.equal(dim(mnist_mp2c2$test$image), c(10L, 10L, 10000L))
})


