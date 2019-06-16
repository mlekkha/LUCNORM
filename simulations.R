# SIMULATIONS
# v0.1
# 
# Author: Manu
# 05/17/19

# A script to test the performace of different methods of normalizing
# Firefly Luciferase and Renilla Luciferase luminescence data. Generates
# synthetic data repeatedly according to specified parameters (sample size,
# number of experiments, standard deviations etc) and then determines the
# Firefly Luciferase activity according to each method:
#
#       1. Ratio of Firefly/Renilla (common method)
#
#       2. Least-squares regression
#
#       3. Errors-in-variables (EIV) regression 
#
#       4. Robust EIV regression (according to Zamar, 1989)
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

# source the functions to be used
source("simulateReporterExperiment.R")
source("genSyntheticData.R")
source("robusteiv.R")
source("eiv.R")

# output file name
outfname <- "results.csv"

# initialize the random seed
set.seed(seed=NULL)

# number of simulations
numSim <- 300

# number of samples in each simulation
numSamples <- 10

# Set R, since it's actual value doesn't matter. 
R <- 50

# Random activities (chosen uniformly between 1 and 50)
#As <- runif(numSim, min = 1, max = 50)

# Choose the activity level
A <- 10

# sds for the contaminated Gaussian model. The standard deviations are
# scaled by the activity to maintain the same level of relative error in
# Firefly and Renilla synthetic data
sd11 <- 3
sd12 <- 3
sd21 <- A * sd11
sd22 <- A * sd12

# how much contamination
e <- 0.05

# different transfection efficiencies correspond to different Beta
# distributions shape parameter combinations. In the beta_shapes matrix,
# parameters a and b are in the 1st and 2nd column respectively. The mean
# of the Beta, and hence the transfection efficiency, is a/(a+b). The
# values below correspond to mean efficiencies of 0.1, 0.25, 0.5, 0.75, and
# 0.9.
beta_shapes <- cbind(c(2,2,2,6,18),c(18,6,2,2,2))

# resulting transfection efficiencies
teffs <- beta_shapes[,1]/rowSums(beta_shapes)

# column names for results
results_colnames <- c("SimType", "Measure", "TxnEff", "N", "R", "A",
                      "sd11", "sd12", "sd21", "sd22", "Contamination",
                      "Ratio", "OLS", "EIV", "REIV")

# FIRST SIMULATION: vary transfection efficiency

# Initialize variables to store results
crit <- data.frame("Efficiency", "Zamar", teffs, numSamples, R, A, sd11, 
                   sd12, sd21, sd22, e, NA, NA, NA, NA)
colnames(crit) <- results_colnames

err <- data.frame("Efficiency", "Error", teffs, numSamples, R, A, sd11, 
                   sd12, sd21, sd22, e, NA, NA, NA, NA)
colnames(err) <- results_colnames

precision <- data.frame("Efficiency", "Precision", teffs, numSamples, R, A,
                        sd11, sd12, sd21, sd22, e, NA, NA, NA, NA)
colnames(precision) <- results_colnames

for (j in seq(1,5))
{

    beta_a <- beta_shapes[j,1]
    beta_b <- beta_shapes[j,2]

    simParams <- matrix(
                        rep(c(numSamples, 
                               R, 
                               A,
                               beta_a,
                               beta_b,
                               sd11,
                               sd12,
                               sd21,
                               sd22,
                               e), 
                          numSim), 
                  nrow=numSim, byrow=TRUE)


    # carry out simulations and get activities
    Aps <- apply(simParams, 1, simulateReporterExperiment)

    # compute the performance criterion according to Zamar (1989)
    Asr <- t(replicate(4, rep(A, numSim)))

    mmatrix <- abs(1 + Aps * Asr)
    mmatrix <- mmatrix/sqrt(1 + Aps * Aps)
    mmatrix <- mmatrix/sqrt(1 + Asr * Asr)
    mmatrix <- 1 - mmatrix

    # The Zamar performance criterion. 
    crit[j,12:15] <- rowSums(mmatrix)/numSim

    # The mean absolute difference as percentage of theoretical activity
    # A
    accuracy <- abs(Aps - Asr)/Asr
#    err[j,12:15] <- rowMeans(accuracy)
    err[j,12:15] <- apply(accuracy, 1, quantile, probs=0.9)
    precision[j,12:15] <- apply(Aps, 1, mad)/A

}

# combine into a single data frame and write out
results <- rbind(crit, err, precision)

# SECOND SIMULATION: vary relative activity of the Firefly reporter. Fix
# transfection efficiency to 0.25

Alevels <- c(1.0, 3.0, 10.0, 30.0, 100.0)
beta_a <- beta_shapes[2,1]
beta_b <- beta_shapes[2,2]
teff <- beta_a/(beta_a + beta_b)

# Initialize variables to store results
crit <- data.frame("Activity", "Zamar", teff, numSamples, R, Alevels,
                   sd11, sd12, sd11*Alevels, sd12*Alevels, e, NA, NA, NA, NA)

colnames(crit) <- results_colnames

err <- data.frame("Activity", "Error", teff, numSamples, R, Alevels,
                  sd11, sd12, sd11*Alevels, sd12*Alevels, e, NA, NA, NA, NA)

colnames(err) <- results_colnames

precision <- data.frame("Activity", "Precision", teff, numSamples, R,
                        Alevels, sd11, sd12, sd11*Alevels, sd12*Alevels, e,
                        NA, NA, NA, NA)

colnames(precision) <- results_colnames

for (j in seq(1,5))
{

    A <- Alevels[j]

    # the noise in the Firefly luciferase data scales with activity
    sd21 <- A * sd11
    sd22 <- A * sd12

    simParams <- matrix(
                        rep(c(numSamples, 
                               R, 
                               A,
                               beta_a,
                               beta_b,
                               sd11,
                               sd12,
                               sd21,
                               sd22,
                               e), 
                          numSim), 
                  nrow=numSim, byrow=TRUE)

    # carry out simulations and get activities
    Aps <- apply(simParams, 1, simulateReporterExperiment)

    # compute the performance criterion according to Zamar (1989)
    Asr <- t(replicate(4, rep(A, numSim)))

    mmatrix <- abs(1 + Aps * Asr)
    mmatrix <- mmatrix/sqrt(1 + Aps * Aps)
    mmatrix <- mmatrix/sqrt(1 + Asr * Asr)
    mmatrix <- 1 - mmatrix

    # The Zamar performance criterion. 
    crit[j,12:15] <- rowSums(mmatrix)/numSim

    # The mean absolute difference as percentage of theoretical activity
    # A
    accuracy <- abs(Aps - Asr)/Asr
#    err[j,12:15] <- rowMeans(accuracy)
    err[j,12:15] <- apply(accuracy, 1, quantile, probs=0.9)
    precision[j,12:15] <- apply(Aps, 1, mad)/A

}

# combine into a single data frame and write out
results <- rbind(results, crit, err, precision)

# THIRD SIMULATION: vary sample size of the Firefly reporter. Fix
# transfection efficiency to 0.25, and activity to 10.

# Reset activity level to 10. The noise level of Firefly has to be reset as
# well, since it depends on A.
A <- 10.0
sd21 <- A * sd11
sd22 <- A * sd12

beta_a <- 2
beta_b <- 6
teff <- beta_a/(beta_a + beta_b)

# The sample sizes we want to test
#SampleSizes <- c(3, 10, 20, 30, 100)
SampleSizes <- seq(3, 100, by=5)

# Initialize variables to store results
crit <- data.frame("SampleSize", "Zamar", teff, SampleSizes, R, A,
                   sd11, sd12, sd21, sd22, e, NA, NA, NA, NA)

colnames(crit) <- results_colnames

err <- data.frame("SampleSize", "Error", teff, SampleSizes, R, A,
                   sd11, sd12, sd21, sd22, e, NA, NA, NA, NA)

colnames(err) <- results_colnames

precision <- data.frame("SampleSize", "Precision", teff, SampleSizes, R, A,
                   sd11, sd12, sd21, sd22, e, NA, NA, NA, NA)

colnames(precision) <- results_colnames


for (j in seq(1,20))
{

    numSamples <- SampleSizes[j]

    simParams <- matrix(
                        rep(c(numSamples, 
                               R, 
                               A,
                               beta_a,
                               beta_b,
                               sd11,
                               sd12,
                               sd21,
                               sd22,
                               e), 
                          numSim), 
                  nrow=numSim, byrow=TRUE)

    # carry out simulations and get activities
    Aps <- apply(simParams, 1, simulateReporterExperiment)

    # compute the performance criterion according to Zamar (1989)
    Asr <- t(replicate(4, rep(A, numSim)))

    mmatrix <- abs(1 + Aps * Asr)
    mmatrix <- mmatrix/sqrt(1 + Aps * Aps)
    mmatrix <- mmatrix/sqrt(1 + Asr * Asr)
    mmatrix <- 1 - mmatrix

    # The Zamar performance criterion. 
    crit[j,12:15] <- rowSums(mmatrix)/numSim

    # The mean absolute difference as percentage of theoretical activity
    # A
    accuracy <- abs(Aps - Asr)/Asr
#    err[j,12:15] <- rowMeans(accuracy)
    err[j,12:15] <- apply(accuracy, 1, quantile, probs=0.9)
    precision[j,12:15] <- apply(Aps, 1, mad)/A

}

# combine into a single data frame and write out
results <- rbind(results, crit, err, precision)

# FOURTH SIMULATION: vary the noise in reporter expression. No
# contamination (sd12 = sd11). Fix transfection efficiency to 0.25, and
# activity to 10.

# Reset the number of samples to 10.
numSamples <- 10

# The standard deviations we want to test
sd11s <- seq(1, 10)

# Initialize variables to store results
crit <- data.frame("Noise", "Zamar", teff, numSamples, R, A,
                   sd11s, sd11s, sd11s * A, sd11s * A, e, NA, NA, NA, NA)

colnames(crit) <- results_colnames

err <- data.frame("Noise", "Error", teff, numSamples, R, A,
                   sd11s, sd11s, sd11s * A, sd11s * A, e, NA, NA, NA, NA)

colnames(err) <- results_colnames

precision <- data.frame("Noise", "Precision", teff, numSamples, R, A,
                   sd11s, sd11s, sd11s * A, sd11s * A, e, NA, NA, NA, NA)

colnames(precision) <- results_colnames


for (j in seq(1,10))
{

    # set the noise levels
    sd11 <- sd11s[j]
    sd12 <- sd11
    sd21 <- A * sd11
    sd22 <- A * sd12

    simParams <- matrix(
                        rep(c(numSamples, 
                               R, 
                               A,
                               beta_a,
                               beta_b,
                               sd11,
                               sd12,
                               sd21,
                               sd22,
                               e), 
                          numSim), 
                  nrow=numSim, byrow=TRUE)

    # carry out simulations and get activities
    Aps <- apply(simParams, 1, simulateReporterExperiment)

    # compute the performance criterion according to Zamar (1989)
    Asr <- t(replicate(4, rep(A, numSim)))

    mmatrix <- abs(1 + Aps * Asr)
    mmatrix <- mmatrix/sqrt(1 + Aps * Aps)
    mmatrix <- mmatrix/sqrt(1 + Asr * Asr)
    mmatrix <- 1 - mmatrix

    # The Zamar performance criterion. 
    crit[j,12:15] <- rowSums(mmatrix)/numSim

    # The mean absolute difference as percentage of theoretical activity
    # A
    accuracy <- abs(Aps - Asr)/Asr
#    err[j,12:15] <- rowMeans(accuracy)
    err[j,12:15] <- apply(accuracy, 1, quantile, probs=0.9)
    precision[j,12:15] <- apply(Aps, 1, mad)/A

}

# combine into a single data frame and write out
results <- rbind(results, crit, err, precision)

# FIFTH SIMULATION: vary the level of contaminating noise in reporter
# expression. Fix noise (sd11) to 3. Fix transfection efficiency to 0.25,
# and activity to 10.

# set A to 3
A <- 3

# Reset the noise level (sd11, sd21)
sd11 <- 3
sd21 <- A * sd11

# The standard deviations of the contaminating noise
sd12s <- seq(3, 30, by = 3)

# Larger sample size to increase chances of outliers
numSamples <- 30

# Initialize variables to store results
crit <- data.frame("Contam", "Zamar", teff, numSamples, R, A,
                   sd11, sd12s, sd21, sd12s * A, e, NA, NA, NA, NA)

colnames(crit) <- results_colnames

err <- data.frame("Contam", "Error", teff, numSamples, R, A,
                   sd11, sd12s, sd21, sd12s * A, e, NA, NA, NA, NA)

colnames(err) <- results_colnames

precision <- data.frame("Contam", "Precision", teff, numSamples, R, A,
                   sd11, sd12s, sd21, sd12s * A, e, NA, NA, NA, NA)

colnames(precision) <- results_colnames

for (j in seq(1,10))
{

    # set the contamination noise levels
    sd12 <- sd12s[j]
    sd22 <- A * sd12

    simParams <- matrix(
                        rep(c(numSamples, 
                               R, 
                               A,
                               beta_a,
                               beta_b,
                               sd11,
                               sd12,
                               sd21,
                               sd22,
                               e), 
                          numSim), 
                  nrow=numSim, byrow=TRUE)

    # carry out simulations and get activities
    Aps <- apply(simParams, 1, simulateReporterExperiment)

    # compute the performance criterion according to Zamar (1989)
    Asr <- t(replicate(4, rep(A, numSim)))

    mmatrix <- abs(1 + Aps * Asr)
    mmatrix <- mmatrix/sqrt(1 + Aps * Aps)
    mmatrix <- mmatrix/sqrt(1 + Asr * Asr)
    mmatrix <- 1 - mmatrix

    # The Zamar performance criterion. 
    crit[j,12:15] <- rowSums(mmatrix)/numSim

    # The mean absolute difference as percentage of theoretical activity
    # A
    accuracy <- abs(Aps - Asr)/Asr
#    err[j,12:15] <- rowMeans(accuracy)
    err[j,12:15] <- apply(accuracy, 1, quantile, probs=0.9)
    precision[j,12:15] <- apply(Aps, 1, mad)/A

}

# combine into a single data frame and write out
results <- rbind(results, crit, err, precision)




# Write out the results
write.csv(results, outfname)
