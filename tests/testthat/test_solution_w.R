## c:/Dropbox/Rpackages/SynthHelpers/tests/testthat/test_solution_w.R

##    Chandler Lutz
##    Questions/comments: cl.eco@cbs.dk
##    $Revisions:      1.0.0     $Date:  2016-11-21

##To test solution_w and make sure it outputs similar data
##as Synth

library(kernlab); library(Rcpp); library(RcppArmadillo); library(microbenchmark)

set.seed(1234)

##sourceCpp("../../src/SynthCpp_functions.cpp")

##load some data from a synth example to check that everything is okay
## solution.v <- readRDS("../../data-raw/solution_v_forc_example.rds")[[1]]
## X.scaled <- readRDS("../../data-raw/X_forc_example.rds")
##Z <- readRDS("../../data-raw/Z_forc_example.rds")

##Run test using synth
data(solution_v_forc_example)
data(X_scaled_forc_example)
data(Z_forc_example)

## treated <- 1
## X0.scaled <- X_scaled_forc_example[, -treated, drop = FALSE]
## X1.scaled <- X_scaled_forc_example[, treated, drop = FALSE]
## Z0 <- Z_forc_example[, -treated, drop = FALSE]
## Z1 <- Z_forc_example[, treated, drop = FALSE]

############# Test solution_w_cpp #############

tol <- 1e-2

##The typical R synth code
f.reg.synth <- function(solution.v, X0.scaled, X1.scaled) {
    ##Synth implementation
    nvarsV <- length(solution.v)

    ##Now that we have V, Solve the linear quadratic programming problem
    ##and get the weights
    V <- diag(x = as.numeric(solution.v), nrow = nvarsV, ncol = nvarsV)

    ##t(X0.scaled) %*% V %*% (X0.scaled)
    ##H <- t(self$X.scaled[,temp.control.ids]) %*% V %*% (self$X.scaled[,temp.control.ids])
    ##Use crossprod to speed things up
    H <- crossprod(X0.scaled,V) %*% X0.scaled
    ##X1.scaled
    a <- X1.scaled

    ##-1 * c(t(a) %*% V %*% (X0.scaled))
    ##which equals
    ##-1 * c(t(X1.scaled) %*% V %*% (X0.scaled))
    c <- -1 * c( crossprod(X1.scaled, V) %*% X0.scaled)

    A <- t(rep(1, length(c)))
    b <- 1;
    l <- rep(0, length(c));
    u <- rep(1, length(c));
    r <- 0;

    res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u,
                r = r, bound = 10, margin = 5e-04,
                maxiter = 1000, sigf = 5)
    ## res2 <- ipopCpp(c = c, H = H, A = A, b = b, l = l, u = u,
    ##             r = r, bound = 10, margin = 5e-04,
    ##             maxiter = 1000, sigf = 5)
    ## print("ipop: ")
    ## print(res@primal)
    ## print("ipopCpp: ")
    ## print(res2$primal)

    ##res2 <- LowRankQP(Vmat=H,dvec=c,Amat=A,bvec=1,uvec=rep(1,length(c)),method="LU")
    ##print(res2$alpha)
    solution.w <- as.matrix(primal(res))  ##row vector
    ##print(dual(res))
    ##print(how(res))
    return(solution.w)
}


## -- Use equal weights for the starting value -- ##

solution.v <- rep(1 / 12, 12)

for (i in 1:ncol(X_scaled_forc_example)){

    treated <- i
    X0.scaled <- X_scaled_forc_example[, -treated, drop = FALSE]
    X1.scaled <- X_scaled_forc_example[, treated, drop = FALSE]

    test_that(paste0("solution_w_cpp() produces the same results as the original Synth Code using equal weights for ", i), {
        expect_equal(f.reg.synth(solution.v, X0.scaled, X1.scaled),
                     solution_w_cpp(as.matrix(solution.v), X0.scaled, X1.scaled),
                     tolerance = tol)
    })

}

## -- Use random starting values -- ##

solution.v <- runif(length(solution_v_forc_example))
solution.v <- solution.v / sum(solution.v)

for (i in 1:ncol(X_scaled_forc_example)) {

    treated <- i
    X0.scaled <- X_scaled_forc_example[, -treated, drop = FALSE]
    X1.scaled <- X_scaled_forc_example[, treated, drop = FALSE]

    test_that(paste0("solution_w_cpp() produces the same results as the original Synth Code using random weights for ", i), {
        expect_equal(f.reg.synth(solution.v, X0.scaled, X1.scaled),
                     solution_w_cpp(as.matrix(solution.v), X0.scaled, X1.scaled),
                     tolerance = tol)
    })

}

## -- Using the actual solution from the Synth package -- ##

solution.v <- solution_v_forc_example

for (i in 1:ncol(X_scaled_forc_example)) {

    treated <- i
    X0.scaled <- X_scaled_forc_example[, -treated, drop = FALSE]
    X1.scaled <- X_scaled_forc_example[, treated, drop = FALSE]

    test_that(paste0("solution_w_cpp() produces the same results as the original Synth Code using the real solution.v for ", i), {
        expect_equal(f.reg.synth(solution.v, X0.scaled, X1.scaled),
                     solution_w_cpp(as.matrix(solution.v), X0.scaled, X1.scaled),
                     tolerance = tol)
    })

}




## microbenchmark(
##     f.reg.synth(solution.v, X0.scaled, X1.scaled),
##     solution_w_cpp(as.matrix(solution.v), X0.scaled, X1.scaled)
## )

############# Test fn_v_cpp #############

tol <- 1e-03

##the original R code
fn_V <- function(
    variables.v = stop("variables.v missing"), X0.scaled = stop("X0.scaled missing"),
    X1.scaled = stop("X1.scaled missing"), Z0 = stop("Z0 missing"),
    Z1 = stop("Z1 missing")) {

    ##Use the ipop() function to minimize (X1 - X0*W)'V(X1-X0*W)
    V <- diag(x = as.numeric(abs(variables.v) / sum(abs(variables.v))))
    ##H <- t(X0.scaled) %*% V %*% (X0.scaled)
    H <- crossprod(X0.scaled,V) %*% X0.scaled  ##Use crossprod to speed things up
    a <- X1.scaled
    c <- -1 * c( (t(a) %*% V) %*% X0.scaled)
    A <- t(rep(1, length(c)))
    b <- 1
    l <- rep(0, length(c))
    u <- rep(1, length(c))
    r <- 0

    res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u,
                r = r, bound = 10, margin = 5e-04,
                maxiter = 1000, sigf = 5)
    solution.w <- as.matrix(primal(res))

    ##loss.v <- as.numeric(t(Z1 - (Z0 %*% solution.w)) %*% (Z1 - (Z0 %*% solution.w)))
    ##Use crossprod() to speed things up
    loss.v <- as.numeric(crossprod(Z1 - (Z0 %*% solution.w)))

    loss.v <- loss.v/nrow(Z0)
    return(loss.v)
}

## -- Use equal weights for the starting value -- ##

solution.v <- rep(1 / 12, 12)

for (i in 1:ncol(X_scaled_forc_example)) {

    treated <- i
    X0.scaled <- X_scaled_forc_example[, -treated, drop = FALSE]
    X1.scaled <- X_scaled_forc_example[, treated, drop = FALSE]
    Z0 <- Z_forc_example[, -treated, drop = FALSE]
    Z1 <- Z_forc_example[, treated, drop = FALSE]

    result <- tryCatch(
        fn_V(solution.v, X0.scaled, X1.scaled, Z0, Z1)
    )
    if (inherits(result, "error")) {
        test_that(paste0("fn_v_cpp() produces the same results as the original Synth Code using equal weights for ", i), {
            expect_error(fn_v_cpp(as.matrix(solution.v), X0.scaled, X1.scaled, Z0, Z1))
        })
    }else {
        test_that(paste0("fn_v_cpp() produces the same results as the original Synth Code using equal weights for ", i), {
            expect_equal(fn_V(solution.v, X0.scaled, X1.scaled, Z0, Z1),
                         fn_v_cpp(as.matrix(solution.v), X0.scaled, X1.scaled, Z0, Z1),
                         tolerance = tol)
        })
    } #end of inherits(result, "error") if
}



## -- Use random weights for the starting value -- ##

## solution.v <- runif(length(solution_v_forc_example))
## solution.v <- solution.v / sum(solution.v)
## solution.v.random <- solution.v

## for (i in 1:ncol(X_scaled_forc_example)) {

##     treated <- i
##     X0.scaled <- X_scaled_forc_example[, -treated, drop = FALSE]
##     X1.scaled <- X_scaled_forc_example[, treated, drop = FALSE]
##     Z0 <- Z_forc_example[, -treated, drop = FALSE]
##     Z1 <- Z_forc_example[, treated, drop = FALSE]

##     result <- tryCatch(
##         fn_V(solution.v, X0.scaled, X1.scaled, Z0, Z1)
##     )
##     if (inherits(result, "error")) {
##         test_that(paste0("fn_v_cpp() produces the same results as the original Synth Code using random weights for ", i), {
##             expect_error(fn_v_cpp(as.matrix(solution.v), X0.scaled, X1.scaled, Z0, Z1))
##         })
##     }else {
##         test_that(paste0("fn_v_cpp() produces the same results as the original Synth Code using random weights for ", i), {
##             expect_equal(fn_V(solution.v, X0.scaled, X1.scaled, Z0, Z1),
##                          fn_v_cpp(as.matrix(solution.v), X0.scaled, X1.scaled, Z0, Z1),
##                          tolerance = tol)
##         })
##     } #end of inherits(result, "error") if
## }


## -- Use real solution for the starting value -- ##

solution.v <- solution_v_forc_example

for (i in 1:ncol(X_scaled_forc_example)) {

    treated <- i
    X0.scaled <- X_scaled_forc_example[, -treated, drop = FALSE]
    X1.scaled <- X_scaled_forc_example[, treated, drop = FALSE]
    Z0 <- Z_forc_example[, -treated, drop = FALSE]
    Z1 <- Z_forc_example[, treated, drop = FALSE]

    result <- tryCatch(
        fn_V(solution.v, X0.scaled, X1.scaled, Z0, Z1)
    )
    if (inherits(result, "error")) {
        test_that(paste0("fn_v_cpp() produces the same results as the original Synth Code using true weights for ", i), {
            expect_error(fn_v_cpp(as.matrix(solution.v), X0.scaled, X1.scaled, Z0, Z1))
        })
    }else {
        test_that(paste0("fn_v_cpp() produces the same results as the original Synth Code using true weights for ", i), {
            expect_equal(fn_V(solution.v, X0.scaled, X1.scaled, Z0, Z1),
                         fn_v_cpp(as.matrix(solution.v), X0.scaled, X1.scaled, Z0, Z1),
                         tolerance = tol)
        })
    } #end of inherits(result, "error") if
}



test_that("fn_v_cpp() produces the same results as the original Synth code", {
    expect_equal(fn_V(solution.v, X0.scaled, X1.scaled, Z0, Z1),
                 fn_v_cpp(as.matrix(solution.v), X0.scaled, X1.scaled, Z0, Z1),
                 tolerance = 1e-5
                 )
})


## treated <- 4
## solution.v <- rep(1 / 12, 12)
## X0.scaled <- X_scaled_forc_example[, -treated, drop = FALSE]
## X1.scaled <- X_scaled_forc_example[, treated, drop = FALSE]
## Z0 <- Z_forc_example[, -treated, drop = FALSE]
## Z1 <- Z_forc_example[, treated, drop = FALSE]
## microbenchmark(
##     fn_V(solution.v, X0.scaled, X1.scaled, Z0, Z1),
##     fn_v_cpp(as.matrix(solution.v), X0.scaled, X1.scaled, Z0, Z1)
## )


