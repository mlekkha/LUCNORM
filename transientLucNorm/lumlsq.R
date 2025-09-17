# lumlsq()
# Wrapper function for the least squares method of normalizing luciferase
# activity.  Takes an n x 2 matrix as input, and computes the mean activity
# by performing an ordinary least squares regression between Firefly
# luciferase luminescence and Renilla luciferase luminescence.


lumlsq <- function(lumdata)
{

    lmfit <- lm(lumdata$Luc ~ lumdata$Ren - 1)
    activity <- coefficients(lmfit)

    return(activity)

}
