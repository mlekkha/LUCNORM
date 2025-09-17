# lumrobusteiv()

# Wrapper function for the robust errors-in-variables method of normalizing
# luciferase activity.  Takes an n x 2 matrix as input, and computes the
# mean activity by performing a robust errors-in-variables regression
# between Firefly luciferase luminescence and Renilla luciferase
# luminescence.


lumrobusteiv <- function(lumdata)
{

    reivfit <- robusteiv(lumdata$Ren, lumdata$Luc, ci_opts=NULL)

    return(reivfit$b)

}
