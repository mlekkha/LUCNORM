eiv <- function(x, y, alpha = 0.05, calcIntercept = FALSE) {

    a = NA
    b = NA

    n = length(x)

    x = as.double(x)
    y = as.double(y)

    if (calcIntercept)
    {

        sxx = sum((x-mean(x))^2)
        syy = sum((y-mean(y))^2)
        sxy = sum((x-mean(x))*(y-mean(y)))

    } else {

        sxx = sum(x*x)
        syy = sum(y*y)
        sxy = sum(x*y)

    }

    b = 0.5 * (1/sxy) * ( (syy - sxx) + sqrt((sxx-syy)^2 + 4*sxy*2)) 

    sigmab_num = (1+b^2)^2 * (sxx*syy - sxy^2)
    sigmab_den = (sxx - syy)^2 + 4*sxy^2

    sigmab = sqrt(sigmab_num/sigmab_den)

    b_lower = b + (1/sqrt(n)) * sigmab * qnorm(alpha/2)
    b_upper = b + (1/sqrt(n)) * sigmab * qnorm(1-alpha/2)


    ci = c(b_lower, b_upper)

    # Gleser's CI

    b_lower_g = b + (1/sqrt(n-1)) * sigmab * qt(alpha/2, n-1)
    b_upper_g = b + (1/sqrt(n-1)) * sigmab * qt(1-alpha/2, n-1)

    gleser_ci = c(b_lower_g, b_upper_g)


    if (calcIntercept)
       a = mean(y) - b*mean(x) 

    return(list(a = a,b = b, ci = ci, gleser_ci = gleser_ci))

}
