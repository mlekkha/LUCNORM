# PLOT_DATATEST
# v0.1
# 
# Author: Manu
# 05/28/19

# A script to plot the results produced by datatest1.R

# source robusteiv
source("robusteiv.R")

# require the necessary packages
require(reshape2)
require(ggplot2)

# filenames to be read
infname3 <- "testdata/puer_cebpa_prom_luc_data.csv"

if (!file.exists(infname3))
{

    stop("Results file ", infname3,  " not found!")

} 

lucdat <- read.csv(infname3)

# Make sure that "Experiment" is a factor variable
lucdat$Experiment <- as.factor(lucdat$Experiment)

# Get experiment 1
lucdat <- lucdat[lucdat$Experiment == 1,]

# get the REIV slope and confidence intervals
reivfit <- robusteiv(lucdat$Ren, lucdat$Luc)

# Add column with labels for outliers
lucdat <- cbind(lucdat, c("", "*", "", "", "", "*", "", "", "", ""))
colnames(lucdat) <- c(colnames(lucdat)[1:8], "Outliers")

pl <- ggplot(data=lucdat, aes(x=Ren, y=Luc, label=Outliers))

pl <- pl + geom_point() + 
           geom_text(nudge_x=1, nudge_y=-2, size=8)

pl <- pl + geom_abline(intercept = 0,
                       slope = reivfit$b)

pl <- pl + geom_abline(intercept = 0,
                       slope = reivfit$boot_positive_ci[1],
                       linetype = 2)

pl <- pl + geom_abline(intercept = 0,
                       slope = reivfit$boot_positive_ci[2],
                       linetype = 2)

pl <- pl + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(y="Firefly Luciferase luminescence") + 
            labs(x="Renilla Luciferase luminescence") + 
            theme(axis.text.x=element_text(size=8, family="Helvetica", colour="black")) + 
            theme(axis.text.y=element_text(size=8, family="Helvetica", colour="black")) 

# figure dimensions
figwidth <- 3.54
figheight <- 3.54

svg('lucdata.svg', width=figwidth, height=figheight)

print(pl)

# close the svg device
dev.off()
