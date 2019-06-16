# genSyntheticData()
#
# Arguments:
#
#           params: a vector containing the parameters of the contaminated
#                   errors in variables model:
#
#                   numSamples: number of reporter experiments
#                   R: maximum Renilla luminescence
#                   A: activity of the Firefly construct
#                   beta_a: alpha parameter of the Beta distribution
#                   beta_b: beta parameter of the Beta distribution
#                   s11: std dev of gaussian Renilla errors
#                   s12: std dev of the contamination in Renilla
#                   s21: std dev of gaussian Firefly errors
#                   s22: std dev of the contamination in Firefly
#                   e: probability of contamination
#
# Returns:
#
#           numSamples x 2 matrix of synthetic Renilla/Firefly luminescence
#           data, with Renilla in first column and Firefly in the second
#           one.
#
#
# DESCRIPTION
# -----------
# Generate synthetic Firefly/Renilla luminescence data with numSamples
# samples. 
#
# Renilla (r_i) and Luciferase (l_i) luminescence data are generated using a
# contaminated (for outliers) errors-in-variables model. Transfection
# efficiency is modeled as a Beta-distributed random variable.
#
#           t_i = Beta(a,b),    a > 1, b > 1
#
# Renilla luminescence is 
#
#           r_i = R*t_i + CN(s_{11}^2, s_{12}^2, e),
# 
# where R is the maximum Renilla luminescence achieved at 100% transfection
# efficiency and 
#
#           CN(s_1^2, s_2^2, e) = (1-e)*N(0,s_1^2) + e*N(0,s_2^2)
# 
# represents contaminated Normal errors introduced by factors such as
# viability. 
#
# The Firefly luminescence is 
#
#           l_i = L*t_i + CN(s_{21}^2, s_{22}^2, e),
#
# where L is the maximum Firefly luminescence. The activity of the Firefly
# construct is given as
#
#           A = L/R.
#
# The objective of the normalization is to estimate A.
#

genSyntheticData <- function(params)
{

    # Convert parameters
    numSamples <- params[1]
    R <- params[2]
    A <- params[3]
    beta_a <- params[4]
    beta_b <- params[5]
    sd11 <- params[6]
    sd12 <- params[7]
    sd21 <- params[8]
    sd22 <- params[9]
    e <- params[10]

    # generate transfection efficiencies
    t = rbeta(numSamples, beta_a, beta_b)

    # generate contaminated gaussian errors for Ren

    # Gaussian errors
    err11 <- rnorm(numSamples, mean=0, sd=sd11)

    # contaminating Gaussian errors
    err12 <- rnorm(numSamples, mean=0, sd=sd12)

    # matrix of errors, allowing us to choose according to probability e
    errs <- cbind(err11, err12)

    # do numSamples Bernoulli trials (p = e) to generate column numbers to 
    # choose between the errors and the contamination
    echoose <- sample(c(1,2), numSamples, replace=TRUE, 
                                    prob=c(1-e,e))

    # the contaminated normal for Ren
    cni <- errs[cbind(seq(1,numSamples), echoose)]

    # the random Renilla luminescence deviates
    r <- R * t + cni

    # generate contaminated gaussian errors for Firefly Luc

    # Gaussian errors
    err21 <- rnorm(numSamples, mean=0, sd=sd21)

    # contaminating Gaussian errors
    err22 <- rnorm(numSamples, mean=0, sd=sd22)

    # matrix of errors, allowing us to choose according to probability e
    errs <- cbind(err21, err22)

    # do numSamples Bernoulli trials (p = e) to generate column numbers to 
    # choose between the errors and the contamination
    echoose <- sample(c(1,2), numSamples, replace=TRUE, 
                                    prob=c(1-e,e))

    # the contaminated normal for Firefly luc
    cni <- errs[cbind(seq(1,numSamples), echoose)]

    # the random Renilla luminescence deviates
    l <- R * A * t + cni

    # the synthetic data vector (numSamples x 2)
    rl <- cbind(r,l)

    # filter to remove negative values
    rl_filt <- rl[rl[,1] > 0 & rl[,2] > 0,]

    # return numSamples x 2 matrix
    return(rl_filt)

}
