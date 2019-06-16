getCRMComparisonPlot <- function(data, toCompare, conditions,
                                regmethod=robusteiv, cim="boot_positive_ci",
				legend=TRUE)
{

    require(ggplot2)
    require(RColorBrewer)

    drows <- dim(data)[1]
    dcols <- dim(data)[2]

    # subset conditions
    cnames <- names(conditions)
    numconds <- length(cnames)

    subset_rows <- rep(T, drows)
    for (cname in cnames)
        subset_rows <- subset_rows & (data[,cname] == get(cname,
                                                        conditions))

    subset <- subset(data, Construct %in% toCompare & subset_rows)

    toCompare <- sort(toCompare)
    numCRMs = length(toCompare)
    
    # calculate EIV fits
    fits <- list()
    for (crm in toCompare)
    {

        x <- subset[subset$Construct == crm,]$Ren
        y <- subset[subset$Construct == crm,]$Luc

        fits[[crm]] <- regmethod(x,y)

    }    

    # make the plot

    ylimits = c(0, max(subset$Luc)*1.05)
    xlimits = c(0, max(subset$Ren)*1.05)

    # Get color palette
    palnames <- c("Blues","Reds","Greens","Purples","Oranges","Greys")

    ccol <- NULL
    for (j in 1:numCRMs)
        ccol <- cbind(ccol, rev(brewer.pal(n=8, name = palnames[j])))

    crmlabels <- NULL
    for (crm in toCompare)
    {
        slope <- fits[[crm]]$b
        crmlabels <- c(crmlabels, 
                        paste("CRM", crm, format(round(slope, 2), nsmall=2)))
    }
    # Create plot object
    pl <- ggplot(subset, aes(x=Ren, y=Luc, color=Construct))

    pl <- pl + theme_bw()
    pl <- pl + theme(panel.grid.major = element_blank())
    pl <- pl + scale_x_continuous(limits = xlimits, expand = c(0,0))
    pl <- pl + scale_y_continuous(limits = ylimits, expand = c(0,0))

    pl <- pl + geom_point() 

    if (legend) 
	    pl <- pl + scale_color_manual(values=ccol[5,], labels=crmlabels)
    else
	    pl <- pl + scale_color_manual(values=ccol[5,], labels=NULL)


    # Add regression lines
    for (crm in toCompare)
    {

        pos <- match(crm, toCompare)

        if (is.na(fits[[crm]]$a)) 
            interc <- 0
        else
            interc <- fits[[crm]]$a

        pl <- pl + geom_abline(intercept = interc,
                               slope = fits[[crm]]$b,
                               color = ccol[2,pos])

        pl <- pl + geom_abline(intercept = interc,
                               slope = fits[[crm]][[cim]][1],
                               linetype = 2,
                               color = ccol[5,pos])

        pl <- pl + geom_abline(intercept = interc,
                               slope = fits[[crm]][[cim]][2],
                               linetype = 2,
                               color = ccol[5,pos])
    }                           


    return(pl)

}

getCondComparisonPlot <- function(data, construct, toCompare, 
                                otherconditions,
                                regmethod=robusteiv, cim="boot_positive_ci")
{

    require(ggplot2)
    require(RColorBrewer)

    drows <- dim(data)[1]
    dcols <- dim(data)[2]

    # subset conditions
    cnames <- names(otherconditions)
    numconds <- length(cnames)

    subset_rows <- rep(T, drows)
    for (cname in cnames)
        subset_rows <- subset_rows & (data[,cname] == get(cname,
                                                    otherconditions))

    subset <- subset(data, subset_rows)

    # subset CRM

    subset <- subset(subset, Construct == construct)

    toCompare <- sort(toCompare)
    numCRMs = length(toCompare)
    
    # calculate EIV fits
    fits <- list()
    for (crm in toCompare)
    {

        x <- subset[subset$Construct == crm,]$Ren
        y <- subset[subset$Construct == crm,]$Luc

        fits[[crm]] <- regmethod(x,y)

    }    

    # make the plot

    ylimits = c(0, max(subset$Luc)*1.05)
    xlimits = c(0, max(subset$Ren)*1.05)

    # Get color palette
    palnames <- c("Blues","Reds","Greens","Purples","Oranges","Greys")

    ccol <- NULL
    for (j in 1:numCRMs)
        ccol <- cbind(ccol, rev(brewer.pal(n=8, name = palnames[j])))

    crmlabels <- NULL
    for (crm in toCompare)
    {
        slope <- fits[[crm]]$b
        crmlabels <- c(crmlabels, 
                        paste(crm, format(round(slope, 2), nsmall=2)))
    }
    # Create plot object
    pl <- ggplot(subset, aes(x=Ren, y=Luc, color=Construct))

    pl <- pl + theme_bw()
    pl <- pl + theme(panel.grid.major = element_blank())
    pl <- pl + scale_x_continuous(limits = xlimits, expand = c(0,0))
    pl <- pl + scale_y_continuous(limits = ylimits, expand = c(0,0))

    pl <- pl + geom_point() 
    pl <- pl + scale_color_manual(values=ccol[5,], labels=crmlabels)


    # Add regression lines
    for (crm in toCompare)
    {

        pos <- match(crm, toCompare)

        if (is.na(fits[[crm]]$a)) 
            interc <- 0
        else
            interc <- fits[[crm]]$a

        pl <- pl + geom_abline(intercept = interc,
                               slope = fits[[crm]]$b,
                               color = ccol[2,pos])

        pl <- pl + geom_abline(intercept = interc,
                               slope = fits[[crm]][[cim]][1],
                               linetype = 2,
                               color = ccol[5,pos])

        pl <- pl + geom_abline(intercept = interc,
                               slope = fits[[crm]][[cim]][2],
                               linetype = 2,
                               color = ccol[5,pos])
    }                           


    return(pl)

}



plotCRMComparisonAllConds <- function(data, toCompare, 
                                regmethod=eiv, cim="gleser_ci",
                                ignore=c("Luc","Ren","Construct","Exp","Note"))
{

    # Get combinations of all conditions except Luc, Ren, and Construct
    cond_combinations <- getConditionCombinations(data, 
                                                ignore=ignore)
    
    # Call getCRMComparisonPlot to get plot object for each combination of
    # conditions present in the data
    plots <- list()
    for (j in 1:dim(cond_combinations)[1])
    {

        conditions = as.list(cond_combinations[j,])
        
        p <- getCRMComparisonPlot(data, toCompare, conditions, 
                                    regmethod, cim)

        title <- ""
        for (condition in names(conditions))
            title <- paste(title, condition,"=",
                                    conditions[[condition]])

        p <- p + ggtitle(title)
        p <- p + theme(plot.title = element_text(size=12))

        plots[[length(plots) + 1]] <- p                            


    }

    
    return(plots)

}

getConditionCombinations <- function(data, ignore = c("Luc", "Ren"))
{

    columns <- colnames(data)

    # Build a list of conditions and their possible values 
    conds = list()

    for (col in columns) 
    {
        
        if (col %in% ignore)
            next

        conds[[col]] <- levels(as.factor(data[[col]]))

    }

    # Formulate all possible combinations of condition values
    cond_combinations <- expand.grid(conds)
    numcombs <- dim(cond_combinations)[1]
    numconds <- dim(cond_combinations)[2]

    # Filter out non-existent combinations
    ccf <- cond_combinations
    for (j in 1:numcombs)
    {

       ss <- data

       for (col in names(conds))
       {
            
            ss <- subset(ss, ss[[col]] == cond_combinations[[col]][j])

       }

       if (dim(ss)[1] == 0)
            ccf[j,] <- rep(NA, n=numconds)


    }
   
    cond_combinations <- subset(ccf, !is.na(ccf[,1]))

    return(cond_combinations)

}

calcSlopesCIs <- function(data, alpha = 0.05, 
                                regmethod=eiv, cim="gleser_ci",
                                ignore=c("Luc","Ren","Exp","Note"))
{

    # Get combinations of all conditions except Luc, Ren, and Construct
    cond_combinations <- getConditionCombinations(data, 
                                ignore=ignore)
    
    numcombs = dim(cond_combinations)[1]
    slopes <- rep(NA, numcombs)
    ci_lower <- rep(NA, numcombs)
    ci_upper <- rep(NA, numcombs)

    # Calculate the slope and CI for 
    for (j in 1:numcombs)
    {
        
       ss <- data

       for (col in colnames(cond_combinations))
       {
            
            ss <- subset(ss, ss[[col]] == cond_combinations[[col]][j])

       }

       x <- ss$Ren
       y <- ss$Luc

       fit <- regmethod(x,y,alpha=alpha)

       slopes[j] <- fit$b
       ci <- fit[[cim]]
       ci_lower[j] <- ci[1]
       ci_upper[j] <- ci[2]

    }

    fits <- as.data.frame(c(cond_combinations,
                        data.frame(slope = slopes),
                        data.frame(ci_lower = ci_lower),
                        data.frame(ci_upper = ci_upper)))

    return(fits)                    
}

getActivityPlot <- function(slopecis, constructs = NULL, 
                          by = "Constructs",
                          conditionorder = list(Cytokine=c("IL3","GCSF"),
                                                OHT=c(FALSE, TRUE)),
			  ptype = "point", dv = 0.25, colorbars=FALSE, freescales=FALSE)
{

    require(reshape2)
    require(ggplot2)

    # Sort the data frame according to the order specified in the
    # conditionorder list as also the factor levels
    for (cond in names(conditionorder))
    {

        # Sort the data, doesn't do anything for plotting
        condsindata <- as.vector(slopecis[[cond]])
        condorder <- getSortByVectorOrder(condsindata, 
                                            conditionorder[[cond]])
        slopecis <- slopecis[condorder, ]

        # Sort the factor levels, essential for the conditions to be
        # plotted in the required order
        slopecis[[cond]] <- factor(slopecis[[cond]], 
                                    levels = conditionorder[[cond]])

    }

    # If a subset or order of constructs is specified then subset and order
    if (!is.null(constructs))
    {
    
	numcons <- length(constructs)
	slopecis <- subset(slopecis, slopecis$Construct %in% constructs)

	    if (by == "Constructs")
	    {
		consindata <- as.vector(slopecis$Construct)
		consorder <- getSortByVectorOrder(consindata, constructs)

		slopecis <- slopecis[consorder, ]

	    }
    }
    # Sort construct factor levels in desired order, required for the
    # constructs to appear in the correct order in the plot
    slopecis$Construct <- factor(slopecis$Construct, levels = constructs)

    # Construct all the combinations of all the 
    # conditions soecified in the conditionorder argument
    combinations <- slopecis[[names(conditionorder)[1]]]

    for (cond in names(conditionorder)[-1])
        combinations <- interaction(combinations, slopecis[[cond]])
    
    slopecis <- as.data.frame(c(slopecis, data.frame(Cond = combinations)))

    # Make long form data frame
    thefactors <- colnames(slopecis)[!(colnames(slopecis) 
                                %in% c("slope","ci_lower","ci_upper"))]

    slcislong <- melt(slopecis, id.vars=thefactors)
    limits <- aes(ymax = slopecis$ci_upper, 
                            ymin = slopecis$ci_lower)

    if (by == "Constructs")
    {
    
        # Create plot with Construct on the x axis
        pl <- ggplot(data=slcislong[slcislong$variable=="slope",],
                            aes(x = Construct, y = value, 
				color = Cond, fill = Cond))

	if (ptype == "point")
        {

		# plot with points
		pl <- pl + geom_point(position=position_dodge(width=dv))

	} else if (ptype == "bar")
	{

		# plot with bars
		pl <- pl + geom_bar(stat = "identity", 
				position=position_dodge(width=dv),
				width = 0.5)

	}

        # error bars
	if (colorbars==FALSE)        
		pl <- pl + geom_errorbar(limits, width = 0,
                            position=position_dodge(width=dv), colour="black")
        else 
		pl <- pl + geom_errorbar(limits, width = 0,
                            position=position_dodge(width=dv))

        # axis labels
        pl <- pl + labs(y = "Activity")
    }

    if (by == "Conditions")
    {

        pl <- ggplot(data=slcislong[slcislong$variable=="slope",],
                            aes(x = Construct, y = value, 
				    color = Cond, fill = Cond))

        # plot with points
        if (ptype == "point")
            pl <- pl + geom_point()
        else
            # plot with bars
            pl <- pl + geom_bar(stat = "identity", 
                    position=position_dodge(width=dv),
                    width = 0.5)

			
        # error bars
        if (colorbars==FALSE)
                pl <- pl + geom_errorbar(limits, width = 0, colour="black")
        else
            pl <- pl + geom_errorbar(limits, width = 0)		

        # facets for conditions
        if (!freescales)
            pl <- pl + facet_wrap(~ Cond, nrow = 1)
        else
            pl <- pl + facet_wrap(~ Cond, scales = "free", drop=TRUE)

    }

    # Change the looks

    # BW theme
    pl <- pl + theme_bw()

    # get rid of major grid lines
    pl <- pl + theme(panel.grid.major = element_blank())

}

getSortByVectorOrder <- function(tosort, sortby)
{

        numby <- length(sortby)

        neword <- NULL
        for (j in 1:numby)
           neword <- c(neword, which(match(tosort,sortby) %in% j))

        return(neword)
}
