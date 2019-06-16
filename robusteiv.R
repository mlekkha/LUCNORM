# Tukey's loss function. Note that the function listed in Eq. 5 in Zamar
# (1989) is incorrect! This is from Microsoft's site
# (http://research.microsoft.com/en-us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html).
tukey <- function(t,c){

   x <- t/c
   xsq <- x*x
   xfourth <- xsq*xsq

   tukey1 <- c^2*xsq*(3 - 3*xsq + xfourth)/6

   return(min(c^2/6, tukey1))

}   

# Derivative of Tukey's loss function. Note that the function listed in Eq.
# 5 in Zamar (1989) is incorrect! This is from Microsoft's site
# (http://research.microsoft.com/en-us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html).
tukey_deriv <- function(t, c) { 

    
   x <- t/c
   xsq <- x*x
   xfourth <- xsq*xsq

   tukey_d_1 <- c*x*(1 - 2*xsq + xfourth)

   if (abs(x) > 1)
        return(0)
   else 
        return(tukey_d_1)

}

# Function to get orthogonal distance of points from line specified by the
# unit vector specified in polar coordinates
getOrthDist <- function(theta, y) {

    # calculate orthogonal distances of the points y to a, the unit vector
    a <- c(cos(theta), sin(theta))
    e <- y %*% a

    return(e)

}

# Given S, a vector of the dot product of each data point with the unit
# vector, e, c, and beta, this function computes function O(s) = sum_i_n
# (chi(e_i/S) - n*beta, which is Eq. 4 in Zamar (1989).
obj_s <- function(S, theta, y, c, beta) {

    n <- dim(y)[1]

    # calculate orthogonal distances of the points y to a, the unit vector
    e <- getOrthDist(theta, y)

    # get first term in above equation
    chi_terms <- sapply(e/S, tukey, c = c)
    chi_sum <- sum(chi_terms)

    return(chi_sum/n - beta)

}    

# Given the angle theta of the unit vector, the points y, c, and beta, this
# function solves Eq. 4 in Zamar (1989), sum_i_n (chi(e_i/S) = n*beta, for
# S.
getS <- function(theta, y, c, beta, maxiter = 1000) {

    # calculate orthogonal distances of the points y to a, the unit vector
    e <- getOrthDist(theta, y)

    # First set the lower and upper limits for S. LHS of Eq. 4 is greater
    # than zero when S << e_i and chi(e_i/S) ~ 1. We choose S_lower =
    # 0.01*min(abs(e_i)). RHS of Eq. 4 is less than zero when S >> e_i and
    # chi(e_i/S) ~ 0. We choose S_upper = 100*max(abs(e_i)).

    S_lower <- 0.01*min(abs(e[abs(e) > 0]))
    S_upper <- 100*max(abs(e))

    # call uniroot to solve for S
    sol <- uniroot(obj_s, interval = c(S_lower, S_upper), theta = theta, 
                                        y = y, c = c,
                                        beta = beta, maxiter = maxiter)

    if (sol$iter == maxiter)
       stop("uniroot failed to converge in getS")

    return(sol$root)

}    

# The derivative of S(theta) w.r.t. theta. According to the displayed
# equation after Eq. 4 in Zamar (1989).
getSDeriv <- function(theta, y, c, beta) {

    # calculate orthogonal distances of the points y to a, the unit vector
    e <- getOrthDist(theta, y)

    # Get the S
    S <- getS(theta = theta, y = y, c = c, beta = beta)
    
    # get the tukey derivative at e and S.
    td <- sapply(e/S, tukey_deriv, c = c)

    # get e derivative (z_ij in Zamar)
    a_deriv <- c(-1*sin(theta), cos(theta))
    e_deriv <- y %*% a_deriv

    # calculate numerator
    numerator_terms <- td * e_deriv
    numerator <- sum(numerator_terms)*S

    # denominator
    denominator_terms <- td * e
    denominator <- sum(denominator_terms)

    return(numerator/denominator)

}

# function that returns a list of objective and gradient function w.r.t.
# theta
obj_deriv_S <- function(theta, y, c, beta){

    return(list("objective" = getS(theta = theta, y = y, c = c, 
                                                    beta = beta),
                "gradient" = getSDeriv(theta = theta, y = y, c = c, 
                                                    beta = beta)))

}

# Function to calculate the derivative of the robust EIV objective function
# (Eq. 3 of Zamar) w.r.t. theta.
getObjThetaDeriv <- function(theta, Sn, y, c) {

    n <- dim(y)[1]

    # calculate orthogonal distances of the points y to a, the unit vector
    e <- getOrthDist(theta, y)

    # get the tukey derivative at e and S.
    td <- sapply(e/Sn, tukey_deriv, c = c)

    # get e derivative (z_ij in Zamar)
    a_deriv <- c(-1*sin(theta), cos(theta))
    e_deriv <- y %*% a_deriv

    # calculate numerator
    numerator_terms <- td * e_deriv
    numerator <- sum(numerator_terms)/(Sn*n)

    return(numerator)
    

}

# First step of Zamar's algorithm: find the unit vector that minimizes the
# scale parameter S. This is done in two stages. First, a grid search
# provides a suitable starting point. 
findOptimalThetaS <- function(y, c = 1.56, beta, ngrid = 21,
                        opts = list( "algorithm" = "NLOPT_LD_SLSQP",
                                      "xtol_rel"  = 1.0e-7,
                                      "maxeval"   = 1000) ) {

    library(nloptr)

    # Grid search
    thetas <- seq(0, pi, by = pi/(ngrid-1))
    Ss <- sapply(thetas, getS, y, c = c, beta = beta)
    initial_theta_pos <- match(min(Ss), Ss)
    initial_theta <- thetas[initial_theta_pos]

    # Use nloptr to find theta that minimizes S.
    lb <- 0
    ub <- pi

    res <- nloptr(x0 = initial_theta,
                  eval_f = obj_deriv_S,
                  lb = lb,
                  ub = ub,
                  opts = opts,
                  y = y,
                  c = c,
                  beta = beta)

    if (res$status <= 1) 
        stop("findOptimalThetaS: optimization for theta that minimizes S
                failed. Nloptr status ", res$status)

    return(list(theta = res$solution, Sn = res$objective))

}

# The objective function for robust EIV. Eq. 3 from Zamar (1989). Note that
# e = y*a = x*a - b, since y = x - median(x). This is the same as the
# objective function for determing S, with beta set to 0. The gradient
# (derivative) is different.
objRobustEIV <- function(theta, Sn, y, c) {


        return(list("objective" = obj_s(S = Sn, theta = theta, 
                                        y = y, c = c, beta = 0.0),
                "gradient" = getObjThetaDeriv(theta, Sn = Sn, y = y, 
                                                            c = c)))


}

# Find the optimal theta by minimizing the objective function for theta
# (Eq. 3 of Zamar)
findOptimalEIVTheta <- function(Sn, y, c = 4.7, initial_theta, 
                        opts = list( "algorithm" = "NLOPT_LD_SLSQP",
                                      "xtol_rel"  = 1.0e-7,
                                      "maxeval"   = 1000) ) {

    library(nloptr)

    # Use nloptr to find theta that minimizes the EIV objective function
    lb <- 0
    ub <- pi

    res <- nloptr(x0 = initial_theta,
                  eval_f = objRobustEIV,
                  lb = lb,
                  ub = ub,
                  opts = opts,
                  Sn = Sn,
                  y = y,
                  c = c)

    if (res$status <= 1) 
        stop("findOptimalEIVTheta: optimization for theta that minimizes 
                the objective function failed. Nloptr status ", res$status)

    return(list(theta = res$solution, objective = res$objective))

}

# robusteiv_pe: perform robust EIV regression according to Zamar (1989) to
# give a point estimate of the slope and intercept
robusteiv_pe <- function(data, inds = NULL, calcIntercept = FALSE,
                        robeiv_opts = list("c1" = 1.56,
                                           "c2" = 4.7,
                                           "beta" = 0.05,
                                           "ngrid" = 21),
                        optim_opts = list( "algorithm" = "NLOPT_LD_SLSQP",
                                      "xtol_rel"  = 1.0e-7,
                                      "maxeval"   = 1000) ) {


    n <- dim(data)[1]

    if (is.null(inds)) inds <- seq(1,n)

    x1 <- data[inds,1]
    x2 <- data[inds,2]

    intercept = 0

    # transform x1 and x2 to center them at their medians if calculating
    # intercepts. Not centering them forces intercept to be zero.
    if (calcIntercept)
    {

        y1 <- x1 - median(x1)
        y2 <- x2 - median(x2)

    } else {

        y1 <- x1 
        y2 <- x2

    }

    # Make an nx2 matrix with x1 in first column and x2 in the second
    # column
    y <- matrix(c(y1,y2), nrow = n, ncol = 2)

    # Perform the minimization of S w.r.t. theta to determine the scale
    # factor to use and an initial guess for theta for minimizing the loss
    # function
    res <- findOptimalThetaS(y, c = robeiv_opts$c1, 
                                beta = robeiv_opts$beta,
                                ngrid = robeiv_opts$ngrid,
                                opts = optim_opts)

    # Perform the minimization of the objective function w.r.t. theta
    res1 <- findOptimalEIVTheta(Sn = res$Sn, y = y,
                                            c = robeiv_opts$c2,
                                            initial_theta = res$theta,
                                            opts = optim_opts)

    if (res1$theta == 0)
        res1$theta <- pi


    slope <- -cos(res1$theta)/sin(res1$theta)

    if (calcIntercept) {
        x <- matrix(c(x1,x2), nrow = n, ncol = 2)

        # get orthogonal distances of the uncenteres points x to a, 
        # the unit vector
        e <- getOrthDist(res1$theta, x)
        intercept = median(e)/sin(res1$theta)

    }    

    return(c(intercept, slope, res1$theta))



}

# gencontam: generate EIV gaussian contaminated deviates
gencontam <- function(beta, alpha = 0, contamination = 0.05, 
                    truestd = 0.5, sigma = 5, tau = 0.5, 
                    npoints = 20) {

        X <- rnorm(npoints, 0, 1)
        Y <- alpha*sqrt(1+beta^2) + beta*X

        trueerror <- rnorm(npoints, 0, truestd)
        contam <- rnorm(npoints, 0, sigma)
        trueorcontam <- rbinom(npoints, size = 1, prob = contamination)
        x <- X + (1-trueorcontam)*trueerror + trueorcontam*contam

        trueerror <- rnorm(npoints, 0, truestd)
        contam <- rnorm(npoints, 0, tau)
        trueorcontam <- rbinom(npoints, size = 1, prob = contamination)
        y <- Y + (1-trueorcontam)*trueerror + trueorcontam*contam

        return(matrix(c(x,y), nrow = npoints, ncol = 2))

}

# Test the robust EIV regression w.r.t. OLS and EIV a la Zamar (1989)
testperf <- function(beta_range = c(-5,5), alpha = 0, contamination = 0.05, 
                    truestd = 0.5, sigma = 5, tau = 0.5, numreps = 100, 
                    npoints = 20) {


    source("~/lib/luc_analysis/eiv.R")
    perf <- matrix(nrow = numreps, ncol = 3)

    for (i in 1:numreps) {

        beta <- runif(1, min = beta_range[1], max = beta_range[2])

        d <- gencontam(beta = beta, alpha = alpha, 
                       contamination = contamination, 
                       truestd = truestd, sigma = sigma, 
                       tau = tau, npoints = npoints) 

        if (alpha == 0) {

            fit <- lm(d[,2] ~ d[,1])
            t <- fit$coefficients[[2]]
            perf[i,1] <- 1 - abs(1+t*beta)/sqrt((1+t^2)*(1+beta^2))

            fit <- eiv(d[,1],d[,2])
            t <- fit$b
            if (!is.nan(t))    
                perf[i,2] <- 1 - abs(1+t*beta)/sqrt((1+t^2)*(1+beta^2))

            fit <- robusteiv_pe(d)
            t <- fit[2]
            perf[i,3] <- 1 - abs(1+t*beta)/sqrt((1+t^2)*(1+beta^2))

        } else {
    
            fit <- lm(d[,2] ~ d[,1])
            t <- fit$coefficients[[2]]
            perf[i,1] <- 1 - abs(1+t*beta)/sqrt((1+t^2)*(1+beta^2))

            fit <- eiv(d[,1],d[,2], calcIntercept = TRUE)
            t <- fit$b
            if (!is.nan(t))    
                perf[i,2] <- 1 - abs(1+t*beta)/sqrt((1+t^2)*(1+beta^2))

            fit <- robusteiv_pe(d, calcIntercept = TRUE)
            t <- fit[2]
            perf[i,3] <- 1 - abs(1+t*beta)/sqrt((1+t^2)*(1+beta^2))

        }

    }

    return(list("OLS" = median(perf[,1]), "EIV" = median(perf[,2]), "ROBEIV" = median(perf[,3])))



}

# robusteiv: perform robust EIV regression according to Zamar (1989).
robusteiv <- function(x1, x2, alpha = 0.05, 
                        calcIntercept = FALSE, normalize = FALSE,
                        robeiv_opts = list("c1" = 1.56,
                                           "c2" = 4.7,
                                           "beta" = 0.05,
                                           "ngrid" = 21),
                        optim_opts = list( "algorithm" = "NLOPT_LD_SLSQP",
                                      "xtol_rel"  = 1.0e-7,
                                      "maxeval"   = 1000),
                        ci_opts = list("sim" = "ordinary",
                                       "parallel" = "multicore",
                                       "numreps" = 999,
                                       "type" = "basic") ) {


    library(parallel)
    library(boot)
    
    # normalize if requested. This is useful for very steep slopes where
    # very small changes in the angle cause big changes in the slope. Even
    # small errors in the determination of theta have a big impact on the
    # slope, especially during the calculation of the confidence intervals.
    # Normalizing x1 and x2 to their respective scales puts the slope near
    # 1, where the slope is not so sensitive to theta.
    if (!normalize)
        data <- matrix(c(x1,x2), nrow = length(x1), ncol = 2)
    else {

        x1scale <- max(abs(x1)) 
        x2scale <- max(abs(x2)) 
        data <- matrix(c(x1/x1scale,x2/x2scale), 
                                    nrow = length(x1), ncol = 2)

    }


    res <- robusteiv_pe(data = data,
                        calcIntercept = calcIntercept,
                        robeiv_opts = robeiv_opts,
                        optim_opts = optim_opts)
    

    # do bootstrap estimation of confidence intervals if options are
    # supplied

    if (!is.null(ci_opts))
    {
        ncpus = detectCores()
        boot.out <- boot(data, robusteiv_pe, R = ci_opts$numreps,
                               parallel = ci_opts$parallel,
                               ncpus = ncpus,
                               sim = ci_opts$sim)

        cires <- boot.ci(boot.out = boot.out, conf = 1 - alpha,
                                              type = ci_opts$type, 
                                              index = c(3))

        t_lower <- cires$basic[4] 
        t_upper <- cires$basic[5]

        if (t_upper < pi) 
            s_upper <- -cos(t_upper)/sin(t_upper)
        else
            s_upper <- NA

        if (t_lower < pi) 
            s_lower <- -cos(t_lower)/sin(t_lower)
        else
            s_lower <- NA


        boot_ci <- c(s_lower, s_upper)


        if ((t_upper < pi) && (t_upper >= pi/2))
            s_upper <- -cos(t_upper)/sin(t_upper)
        else
            s_upper <- NA

        if ((t_lower < pi) && (t_lower >= pi/2))
            s_lower <- -cos(t_lower)/sin(t_lower)
        else
            s_lower <- NA

        boot_positive_ci <- c(s_lower, s_upper)

    } else
    {
        boot_ci <- c(0, 0)
        boot_positive_ci <- c(0, 0)
    }

    slope <- res[2]
    intercept <- res[1]

    # correct factors if data were normalized
    if (normalize)
    {

        slope <- slope*x2scale/x1scale
        boot_ci <- boot_ci*x2scale/x1scale
        boot_positive_ci <- boot_positive_ci*x2scale/x1scale

        intercept <- intercept*x2scale

    }

    return(list("a" = intercept, "b" = slope, 
                "boot_ci" = boot_ci,
                "boot_positive_ci" = boot_positive_ci))                        


}
