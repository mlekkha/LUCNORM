# TEST REIV ON DATA
# v0.1
# 
# Author: Manu
# 06/06/19

# A script to test the performace of different methods of normalizing
# Firefly Luciferase and Renilla Luciferase luminescence data. Firefly
# Luciferase expression was driven by the Cebpa proximal promoter, while
# Renilla Luciferase expression was driven by the CMV promoter in PUER
# cells. The dataset contains 85 replicate measurements spread over 8
# independent experiments. The script subsamples this dataset with varying
# sample sizes and then determines the
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
# The script also determines the activity according to each method for each
# individual experiment.

# source the functions to be used
source("robusteiv.R")
source("eiv.R")

# output file name
outfname1 <- "results_datatest1.csv"
outfname2 <- "results_datatest2.csv"

# read in the dataset
lucdat <- read.csv("testdata/puer_cebpa_prom_luc_data.csv")

# Make sure that "Experiment" is a factor variable
lucdat$Experiment <- as.factor(lucdat$Experiment)

# First determine the activity of each experiment according to the methods
experiments <- levels(lucdat$Experiment)
numexp <- length(experiments)

# Initialize data frame to store results
results_colnames <- c("Experiment", "Method", "Activity", "Lower", "Upper")
methods <- c(rep("Ratio", numexp),
             rep("LSQ", numexp),
             rep("EIV", numexp),
             rep("REIV", numexp))

results <- data.frame(rep(experiments,4), methods, NA, NA, NA) 
colnames(results) <- results_colnames

# loop over experiments and compute activity according to various methods
for (experiment in levels(lucdat$Experiment))
{

    lucsubset <- lucdat[lucdat$Experiment == experiment,]
    numSamples <- nrow(lucsubset)

    # Ratio method
    rowtomod <- results$Experiment == experiment & results$Method == "Ratio"

    activities <- lucsubset$Luc/lucsubset$Ren
    activity <- mean(activities)
    activity_stderr <- sd(activities)/sqrt(numSamples)
    results[rowtomod, 3:5] <- c(activity, 
                                activity - 1.96*activity_stderr,
                                activity + 1.96*activity_stderr)

    # Least squares (LSQ) method
    rowtomod <- results$Experiment == experiment & results$Method == "LSQ"

    lmfit <- lm(lucsubset$Luc ~ lucsubset$Ren - 1)
    activity <- coefficients(lmfit)
    ci <- confint(lmfit)
    results[rowtomod, 3:5] <- c(activity[[1]], 
                                ci[1],
                                ci[2])

    # Errors-in-variables (EIV) method
    rowtomod <- results$Experiment == experiment & results$Method == "EIV"

    eivfit <- eiv(lucsubset$Ren, lucsubset$Luc)
    results[rowtomod, 3:5] <- c(eivfit$b, 
                                eivfit$gleser_ci[1],
                                eivfit$gleser_ci[2])

    # Robust errors-in-variables (REIV) method
    rowtomod <- results$Experiment == experiment & results$Method == "REIV"

    reivfit <- robusteiv(lucsubset$Ren, lucsubset$Luc)
    results[rowtomod, 3:5] <- c(reivfit$b, 
                                reivfit$boot_positive_ci[1],
                                reivfit$boot_positive_ci[2])



}

# loop over the number of samples
for (numSamples in c(3, 6, 10, 20, 30))
{



}




# Write out the results
#write.csv(results, outfname1)
