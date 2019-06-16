# simulateReporterExperiment()
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
#           A 2-element vector estA = c(A1, A2), where
#
#                   A1: activity inferred by the average ratio of Firefly
#                   to Renilla luminescence
#
#                   A2: activity inferred by Robust EIV method
#
#
# DESCRIPTION
# -----------
# 
# Given the parameters of a particular simulation (params) calls
# genSyntheticData() to generate Firefly/Renilla deviates according to the
# contaminated Gaussian model and then estimates Firefly activity accoring
# to the ratio or robust EIV method. Returns the estimated activity.

simulateReporterExperiment <- function(params)
{

    numSamples <- params[1]

    # generate the synthetic data
    rl <- genSyntheticData(params)


    # check if we got less than numSamples, if so call genSynthetic data
    # until we have at least numSamples
    while ((length(rl) == 0) | (nrow(matrix(rl,ncol=2)) < numSamples))
    {
        rln <- genSyntheticData(params)
        rl <- rbind(rl, rln)
    }

    # take first numSamples
    rl <- rl[1:numSamples,]

    # estimate the activity according to Firefly/Renilla ratio
    A1 <- mean(rl[,2]/rl[,1])

    # estimate the activity according to least squares regression
    lmfit <- lm(rl[,2] ~ rl[,1] - 1)
    A2 <- lmfit$coefficients[[1]]

    # estimate the activity according to EIV regression
    eivfit <- eiv(rl[,1], rl[,2])
    A3 <- eivfit$b

    # estimate the activity according to robust EIV regression
    robusteiv_result <- robusteiv(rl[,1], rl[,2], ci_opts=NULL)
    A4 <- robusteiv_result$b

    # return the inferred activities
    return(c(A1,A2,A3,A4))

}
