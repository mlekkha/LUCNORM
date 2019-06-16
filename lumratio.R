# ratio()
# Wrapper function for the Ratio method of normalizing luciferase activity.
# Takes an n x 2 matrix as input, and computes the mean activity as the
# ratio of Firefly luciferase luminescence and Renilla luciferase
# luminescence.


lumratio <- function(lumdata)
{
    activities <- lumdata$Luc/lumdata$Ren

    return(mean(activities))

}
