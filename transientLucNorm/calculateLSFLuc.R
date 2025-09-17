calculateLSFLuc <- function(ldata, rdata){

    numcrms <- ncol(ldata)

    crms <- colnames(ldata)

    slopes <- vector(mode="numeric", length=numcrms)
    stderrs <- vector(mode="numeric", length=numcrms)

    for (i in 1:numcrms)
    {
        fit <- lm(ldata[, crms[i]] ~ rdata[, crms[i]] - 1)
        m <- summary(fit)
        slopes[i] <- m$coefficients[, 1]
        stderrs[i] <- m$coefficients[, 2]

    }

    return(matrix(c(slopes, stderrs), nrow = numcrms, ncol = 2)) 

}
