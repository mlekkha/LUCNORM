# PLOT_DATATEST
# v0.1
# 
# Author: Manu
# 05/28/19

# A script to plot the results produced by datatest1.R

# require the necessary packages
require(reshape2)
require(ggplot2)
require(gridExtra)
require(grid)
require(gtable)

# filenames to be read
infname1 <- "datatest_activities.csv"
infname2 <- "datatest_activities_high.csv"
infname3 <- "testdata/puer_cebpa_prom_luc_data.csv"

if (!file.exists(infname1))
{

    stop("Results file ", infname1,  " not found!")

} 

if (!file.exists(infname2))
{

    stop("Results file ", infname2,  " not found!")

} 

if (!file.exists(infname3))
{

    stop("Results file ", infname3,  " not found!")

} 

# Read in the files 
print("Reading the files...")
activities <- read.csv(infname1)
activities$SampleSize <- as.factor(activities$SampleSize)
activities$Method <- factor(activities$Method, 
                               levels=c("Ratio", "OLS", "EIV", "REIV"))

activities_high <- read.csv(infname2)
activities_high$SampleSize <- as.factor(activities_high$SampleSize)
activities_high$Method <- factor(activities_high$Method, 
                               levels=c("Ratio", "OLS", "EIV", "REIV"))

lucdat <- read.csv(infname3)

# Initialize list of plots
pl <- vector("list", 3)

pl[[1]] <- ggplot(data=lucdat, aes(x=Ren, y=Luc))

pl[[1]] <- pl[[1]] + geom_point(size=0.5)

pl[[1]] <- pl[[1]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(y="Firefly Luciferase luminescence") + 
            labs(x="Renilla Luciferase luminescence") + 
            theme(axis.text.x=element_text(size=8, family="Helvetica", colour="black")) + 
            theme(axis.text.y=element_text(size=8, family="Helvetica", colour="black")) + 
            labs(tag = "(a)")

pl[[2]] <- ggplot(data=activities, 
          aes(x = SampleSize, y=Activity, color=Method))


pl[[2]] <- pl[[2]] + geom_boxplot(outlier.size=0.5)

pl[[2]] <- pl[[2]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            coord_cartesian(ylim = c(9, 32)) + 
            labs(y="Inferred activity") + 
            labs(x="Sample size (N)") + 
            theme(axis.text.x=element_text(size=8, family="Helvetica", colour="black")) + 
            theme(axis.text.y=element_text(size=8, family="Helvetica", colour="black")) + 
            labs(tag = "(b)") +
            guides(color=guide_legend(title=NULL)) + 
            theme(legend.position=c(0.95,1),
                  legend.background = element_blank(),
                  legend.justification = c("right", "top"),
                  legend.box.just = "right",
                  legend.margin = margin(2, 2, 2, 2))

# Distributions for activities inferred from experiments with higher
# Renilla expression
pl[[3]] <- ggplot(data=activities_high, 
          aes(x = SampleSize, y=Activity, color=Method))


pl[[3]] <- pl[[3]] + geom_boxplot(outlier.size=0.5)

pl[[3]] <- pl[[3]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            coord_cartesian(ylim = c(9, 32)) + 
            labs(y=element_blank()) + 
            labs(x="Sample size (N)") + 
            theme(axis.text.x=element_text(size=8, family="Helvetica", colour="black")) + 
            theme(axis.text.y=element_blank()) + 
            labs(tag = "(c)") +
            guides(color="none")


# figure dimensions
figwidth <- 7.08

svg('datatest1.svg', width=figwidth)

# make the multipanel grid of objects
mygrobs <- lapply(pl, ggplotGrob)

grid.arrange(
  grobs = mygrobs,
  widths = c(1,1),
  layout_matrix = rbind(c(1, 1),
                        c(2, 3))
)

# close the svg device
dev.off()
