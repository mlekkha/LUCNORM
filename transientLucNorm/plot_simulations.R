# PLOT_SIMULATIONS
# v0.1
# 
# Author: Manu
# 05/28/19

# A script to plot the results produced by simulations.R

# require the necessary packages
require(reshape2)
require(ggplot2)
#require(gridExtra)
require(grid)
require(gtable)

# filename to be read
infname <- "results.csv"

if (!file.exists(infname))
{

    stop("Results file ", infname,  " not found!")

} else {

    print("Reading the file...")
    results <- read.csv(infname)
    
}

# Make long form data frame
thefactors <- colnames(results)[!(colnames(results) 
                                %in% c("Ratio","OLS","EIV","REIV"))]

resultslong <- melt(results, id.vars=thefactors)

# Initialize list of plots
pl <- vector("list", 6)

# Plots for the transfection efficiency simulations
restoplot <- resultslong[resultslong$SimType == "Efficiency",]

pl[[1]] <- ggplot(data=restoplot[restoplot$Measure == "Zamar",], 
          aes(x = TxnEff, y = value, color = variable))


pl[[1]] <- pl[[1]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[1]] <- pl[[1]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(y="Zamar criterion") + 
            labs(x=element_blank()) + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_text(size=8, family="Helvetica", colour="black", margin=margin(t=0, r=0, b=0, l=10))) + 
            coord_cartesian(ylim = c(0, 0.003)) + 
            labs(tag = "(a)") +
            guides(color=guide_legend(title=NULL)) + 
            scale_color_discrete( breaks=c("Ratio", "OLS")) + 
            theme(legend.position=c(0.95,1),
                  legend.direction="horizontal",
                  legend.background = element_blank(),
                  legend.justification = c("right", "top"),
                  legend.box.just = "right",
                  legend.margin = margin(2, 2, 2, 2))


pl[[2]] <- ggplot(data=restoplot[restoplot$Measure == "Error",], 
          aes(x = TxnEff, y = value, color = variable))

pl[[2]] <- pl[[2]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[2]] <- pl[[2]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(y=bquote("Relative bias ("~P[90]~")")) + 
            labs(x=element_blank()) + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_text(size=8, family="Helvetica", colour="black", margin=margin(t=0, r=0, b=0, l=10))) + 
            coord_cartesian(ylim = c(0, 3)) + 
            guides(color="none")

pl[[3]] <- ggplot(data=restoplot[restoplot$Measure == "Precision",], 
          aes(x = TxnEff, y = value, color = variable))

pl[[3]] <- pl[[3]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[3]] <- pl[[3]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(y="Relative MAD") + 
            labs(x=bquote("Mean transfection efficiency ("~bar(t)~")")) + 
            theme(axis.text.x=element_text(size=8, family="Helvetica", colour="black")) + 
            theme(axis.text.y=element_text(size=8, family="Helvetica", colour="black", margin=margin(t=0, r=0, b=0, l=10))) + 
            coord_cartesian(ylim = c(0, 0.75)) + 
            guides(color="none")


# Plots for the activity level simulations
restoplot <- resultslong[resultslong$SimType == "Activity",]

pl[[4]] <- ggplot(data=restoplot[restoplot$Measure == "Zamar",], 
          aes(x = A, y = value, color = variable))


pl[[4]] <- pl[[4]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[4]] <- pl[[4]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=element_blank()) + 
            labs(y=element_blank()) + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_blank()) +
            coord_cartesian(ylim = c(0, 0.003)) + 
            labs(tag = "(b)") +
            guides(color=guide_legend(title=NULL)) + 
            scale_color_discrete( breaks=c("EIV", "REIV")) + 
            theme(legend.position=c(0.95,1),
                  legend.direction="horizontal",
                  legend.background = element_blank(),
                  legend.justification = c("right", "top"),
                  legend.box.just = "right",
                  legend.margin = margin(2, 2, 2, 2))

pl[[5]] <- ggplot(data=restoplot[restoplot$Measure == "Error",], 
          aes(x = A, y = value, color = variable))

pl[[5]] <- pl[[5]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[5]] <- pl[[5]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=element_blank()) + 
            labs(y=element_blank()) + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_blank()) +
            coord_cartesian(ylim = c(0, 3)) + 
            guides(color="none")

pl[[6]] <- ggplot(data=restoplot[restoplot$Measure == "Precision",], 
          aes(x = A, y = value, color = variable))

pl[[6]] <- pl[[6]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[6]] <- pl[[6]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=bquote("Relative activity ("~A~")")) + 
            labs(y=element_blank()) + 
            theme(axis.text.x=element_text(size=8, family="Helvetica", colour="black")) + 
            theme(axis.text.y=element_blank()) +
            coord_cartesian(ylim = c(0, 0.75)) + 
            guides(color="none")

# Plots for the sample size simulations
restoplot <- resultslong[resultslong$SimType == "SampleSize",]

pl[[7]] <- ggplot(data=restoplot[restoplot$Measure == "Zamar",], 
          aes(x = N, y = value, color = variable))


pl[[7]] <- pl[[7]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[7]] <- pl[[7]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=element_blank()) + 
            labs(y=element_blank()) + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_blank()) +
            coord_cartesian(ylim = c(0, 0.003)) + 
            labs(tag = "(c)") +
            guides(color="none")




pl[[8]] <- ggplot(data=restoplot[restoplot$Measure == "Error",], 
          aes(x = N, y = value, color = variable))

pl[[8]] <- pl[[8]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[8]] <- pl[[8]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=element_blank()) + 
            labs(y=element_blank()) + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_blank()) +
            coord_cartesian(ylim = c(0, 3)) + 
            guides(color="none")

pl[[9]] <- ggplot(data=restoplot[restoplot$Measure == "Precision",], 
          aes(x = N, y = value, color = variable))

pl[[9]] <- pl[[9]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[9]] <- pl[[9]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=bquote("Sample size ("~N~")")) + 
            labs(y=element_blank()) + 
            theme(axis.text.x=element_text(size=8, family="Helvetica", colour="black")) + 
            theme(axis.text.y=element_blank()) +
            coord_cartesian(ylim = c(0, 0.75)) + 
            guides(color="none")

# Plots for the noise level simulations
restoplot <- resultslong[resultslong$SimType == "Noise",]

pl[[10]] <- ggplot(data=restoplot[restoplot$Measure == "Zamar",], 
          aes(x = sd11, y = value, color = variable))


pl[[10]] <- pl[[10]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[10]] <- pl[[10]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=element_blank()) + 
            labs(y="Zamar criterion") + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_text(size=8, family="Helvetica", colour="black", margin=margin(t=0, r=0, b=0, l=10))) + 
            coord_cartesian(ylim = c(0, 0.003)) + 
            labs(tag = "(a)") +
            guides(color=guide_legend(title=NULL)) + 
            scale_color_discrete( breaks=c("Ratio", "OLS")) + 
            theme(legend.position=c(0,1),
                  legend.direction="horizontal",
                  legend.background = element_blank(),
                  legend.justification = c("left", "top"),
                  legend.box.just = "right",
                  legend.margin = margin(2, 2, 2, 2))


pl[[11]] <- ggplot(data=restoplot[restoplot$Measure == "Error",], 
          aes(x = sd11, y = value, color = variable))

pl[[11]] <- pl[[11]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[11]] <- pl[[11]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=element_blank()) + 
            labs(y=bquote("Relative bias ("~P[90]~")")) + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_text(size=8, family="Helvetica", colour="black", margin=margin(t=0, r=0, b=0, l=10))) + 
            coord_cartesian(ylim = c(0, 3)) + 
            guides(color="none")

pl[[12]] <- ggplot(data=restoplot[restoplot$Measure == "Precision",], 
          aes(x = sd11, y = value, color = variable))

pl[[12]] <- pl[[12]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[12]] <- pl[[12]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=bquote("Renilla error ("~sigma[11]~")")) + 
            labs(y="Relative MAD") + 
            theme(axis.text.x=element_text(size=8, family="Helvetica", colour="black")) + 
            theme(axis.text.y=element_text(size=8, family="Helvetica", colour="black", margin=margin(t=0, r=0, b=0, l=10))) + 
            coord_cartesian(ylim = c(0, 0.75)) + 
            guides(color="none")

# Plots for the contaminating noise level simulations
restoplot <- resultslong[resultslong$SimType == "Contam",]

pl[[13]] <- ggplot(data=restoplot[restoplot$Measure == "Zamar",], 
          aes(x = sd12, y = value, color = variable))


pl[[13]] <- pl[[13]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[13]] <- pl[[13]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=element_blank()) + 
            labs(y=element_blank()) + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_blank()) +
            coord_cartesian(ylim = c(0, 0.003)) + 
            labs(tag = "(b)") +
            guides(color=guide_legend(title=NULL)) + 
            scale_color_discrete( breaks=c("EIV", "REIV")) + 
            theme(legend.position=c(0,1),
                  legend.direction="horizontal",
                  legend.background = element_blank(),
                  legend.justification = c("left", "top"),
                  legend.box.just = "right",
                  legend.margin = margin(2, 2, 2, 2))

pl[[14]] <- ggplot(data=restoplot[restoplot$Measure == "Error",], 
          aes(x = sd12, y = value, color = variable))

pl[[14]] <- pl[[14]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[14]] <- pl[[14]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=element_blank()) + 
            labs(y=element_blank()) + 
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_blank()) +
            coord_cartesian(ylim = c(0, 3)) + 
            guides(color="none")

pl[[15]] <- ggplot(data=restoplot[restoplot$Measure == "Precision",], 
          aes(x = sd12, y = value, color = variable))

pl[[15]] <- pl[[15]] + geom_line(size=0.5) + geom_point(size=1.5)
pl[[15]] <- pl[[15]] + 
            theme_bw() + 
            theme(plot.background=element_blank()) +
            theme(panel.grid = element_blank()) +
            theme(text=element_text(size=10, family="Helvetica")) + 
            labs(x=bquote("Contaminating error ("~sigma[12]~")")) + 
            labs(y=element_blank()) + 
            theme(axis.text.x=element_text(size=8, family="Helvetica", colour="black")) + 
            theme(axis.text.y=element_blank()) +
            coord_cartesian(ylim = c(0, 0.75)) + 
            guides(color="none")

# make the multipanel grid of objects
mygrobs <- lapply(pl, ggplotGrob)

# first plot for txn efficiency, activity, and sample size
g <- rbind(
           cbind(mygrobs[[1]], mygrobs[[4]], mygrobs[[7]], size="first"),
           cbind(mygrobs[[2]], mygrobs[[5]], mygrobs[[8]], size="first"),
           cbind(mygrobs[[3]], mygrobs[[6]], mygrobs[[9]], size="first"),
           size="first")

# figure dimensions
figwidth <- 7.09
figheight <- 4.68

svg('syntheticdata_sims1.svg', width=figwidth, height=figheight)

# Draw the panels
grid.newpage()
grid.draw(g)


# close the svg device
dev.off()

# second plot for noise and contamination
g <- rbind(
           cbind(mygrobs[[10]], mygrobs[[13]], size="first"),
           cbind(mygrobs[[11]], mygrobs[[14]], size="first"),
           cbind(mygrobs[[12]], mygrobs[[15]], size="first"),
           size="first")

# figure dimensions
figwidth <- 4.99
figheight <- 4.68

svg('syntheticdata_sims2.svg', width=figwidth, height=figheight)


# Draw the panels
grid.newpage()
grid.draw(g)


# close the svg device
dev.off()
