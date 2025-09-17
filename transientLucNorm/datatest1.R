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

# source the functions to be used
source("robusteiv.R")
source("eiv.R")
source("lumratio.R")
source("lumlsq.R")
source("lumeiv.R")
source("lumrobusteiv.R")

# output file name
outfname1 <- "datatest_summary.csv"
outfname2 <- "datatest_activities.csv"

# read in the dataset
lucdat <- read.csv("testdata/puer_cebpa_prom_luc_data.csv")

# Make sure that "Experiment" is a factor variable
lucdat$Experiment <- as.factor(lucdat$Experiment)

# Uncomment for sampling high expression datapoints only.
#lucdat <- lucdat[lucdat$Ren > 8,]
#outfname1 <- "datatest_summary_high.csv"
#outfname2 <- "datatest_activities_high.csv"

# total number of samples
totalSamples <- nrow(lucdat)

# the number of times to subsample
numReplicates <- 1000

# The sample sizes to test
SampleSizes <- c(3, 6, 10, 20, 30)
numSampleSizes <- length(SampleSizes)

# Initialize data frames to store results

# data frame to store activity summary statistics
act_summary_colnames <- c("SampleSize", "Method", "Activity", 
                          "StdDev", "Lower", "Upper", "StdErr", 
                          "CILower", "CIUpper")
methods <- c(rep("Ratio", numSampleSizes),
             rep("OLS", numSampleSizes),
             rep("EIV", numSampleSizes),
             rep("REIV", numSampleSizes))

act_summary <- data.frame(rep(SampleSizes,4), methods, NA, NA, NA, NA, NA, NA, NA) 
colnames(act_summary) <- act_summary_colnames

# data frame to store activities
act_samples_colnames <- c("SampleSize", "Method", "Activity")
methods <- c(rep("Ratio", numSampleSizes*numReplicates),
             rep("OLS", numSampleSizes*numReplicates),
             rep("EIV", numSampleSizes*numReplicates),
             rep("REIV", numSampleSizes*numReplicates))

act_samples <- data.frame(rep(SampleSizes,4*numReplicates), methods, NA) 
colnames(act_samples) <- act_samples_colnames



# loop over the number of samples
for (SampleSize in SampleSizes)
{

    # sample SampleSize points numReplicates times
    samples <- replicate(numReplicates, 
                                lucdat[sample(seq(1:totalSamples), 
                                           SampleSize, replace=T),
                                    1:2], 
                        simplify=F)

    # Ratio method
    activities <- sapply(samples, lumratio)

    # save summary statistics
    rowtomod <- act_summary$SampleSize == SampleSize & act_summary$Method == "Ratio"
    act_summary[rowtomod, 3:9] <- as.numeric(
                                 c(mean(activities), 
                                   sd(activities),
                                   mean(activities) - sd(activities),
                                   mean(activities) + sd(activities),
                                   sd(activities)/sqrt(length(activities)),
                                quantile(activities, probs=c(0.025, 0.975)))
                                        )

    # save raw activities
    rowtomod <- act_samples$SampleSize == SampleSize & act_samples$Method == "Ratio"
    act_samples[rowtomod, 3] <- activities

    # Least squares (OLS) method
    activities <- sapply(samples, lumlsq)

    # save summary statistics
    rowtomod <- act_summary$SampleSize == SampleSize & act_summary$Method == "OLS"
    act_summary[rowtomod, 3:9] <- as.numeric(
                                 c(mean(activities), 
                                   sd(activities),
                                   mean(activities) - sd(activities),
                                   mean(activities) + sd(activities),
                                   sd(activities)/sqrt(length(activities)),
                                quantile(activities, probs=c(0.025, 0.975)))
                                        )

    # save raw activities
    rowtomod <- act_samples$SampleSize == SampleSize & act_samples$Method == "OLS"
    act_samples[rowtomod, 3] <- activities

    # Errors-in-variables (EIV) method
    activities <- sapply(samples, lumeiv)

    # save summary statistics
    rowtomod <- act_summary$SampleSize == SampleSize & act_summary$Method == "EIV"
    act_summary[rowtomod, 3:9] <- as.numeric(
                                 c(mean(activities), 
                                   sd(activities),
                                   mean(activities) - sd(activities),
                                   mean(activities) + sd(activities),
                                   sd(activities)/sqrt(length(activities)),
                                quantile(activities, probs=c(0.025, 0.975)))
                                        )
    # save raw activities
    rowtomod <- act_samples$SampleSize == SampleSize & act_samples$Method == "EIV"
    act_samples[rowtomod, 3] <- activities

    # Robust errors-in-variables (REIV) method
    activities <- sapply(samples, lumrobusteiv)

    # save summary statistics
    rowtomod <- act_summary$SampleSize == SampleSize & act_summary$Method == "REIV"
    act_summary[rowtomod, 3:9] <- as.numeric(
                                 c(mean(activities), 
                                   sd(activities),
                                   mean(activities) - sd(activities),
                                   mean(activities) + sd(activities),
                                   sd(activities)/sqrt(length(activities)),
                                quantile(activities, probs=c(0.025, 0.975)))
                                        )

    # save raw activities
    rowtomod <- act_samples$SampleSize == SampleSize & act_samples$Method == "REIV"
    act_samples[rowtomod, 3] <- activities

}

# Write out the results
write.csv(act_summary, outfname1)
write.csv(act_samples, outfname2)
