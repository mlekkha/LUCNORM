# lumeiv()
# Wrapper function for the errors-in-variables method of normalizing
# luciferase activity.  Takes an n x 2 matrix as input, and computes the
# mean activity by performing an EIV regression between Firefly luciferase
# luminescence and Renilla luciferase luminescence.

lumeiv <- function(lumdata)
{

    eivfit <- eiv(lumdata$Ren, lumdata$Luc)

    return(eivfit$b)

}
