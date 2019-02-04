#' Relative abundance calculation 
#'
#' Calculate the relative abundance (fraction) from a table of counts
#' @param counts counts table (rownames are features)
#' @return a data frame
#' @examples
#' relab(counts.table)

relab <- function(counts){

    relab <- (sweep(counts, 2, colSums(counts), "/"))
    return(relab)
}

#' Check the type of abundance table
#'
#' Check the nature of abundance (count, fraction or percentage)
#' @param x data frame of abundances
#' @return a string (fraction, percentage or count)
#' @examples
#' checkValue(counts.table)

checkValue <- function(x){

    # use column 1 to check for what the abundance
    # is measured in - could be counts, fraction or
    # percentage

    # use sum to check
    cat("checking abundance type")
    if (sum(x[,1]) == 1){
        abundance.type <- "fraction"}

    else if (sum(x[,1]) == 100){
	abundance.type <- "percentage"}
    else{
        abundance.type <- "count"}
    return(abundance.type)
}

#' Read in an abundance file
#'
#' Reads in data file
#' @param abundanceFile .tsv file with first column as feature. Columns are samples
#' @return a data frame
#' @examples
#' readAbundanceFile("abundances.tsv")
#' @export

readAbundanceFile <- function(abundanceFile){

    dat <- read.csv(abundanceFile,
                    header=T,
		    stringsAsFactors=F,
		    sep="\t",
		    row.names=1)

    abundance.type <- checkValue(dat)
    cat(paste0("; abundance type = ", abundance.type, "\n"))

    # convert if not fraction
    if (abundance.type == "percentage"){
        dat <- dat/100}
    else if (abundance.type == "count"){
        dat <- relab(dat)}

    nrows <- nrow(dat)
    nsamples <- ncol(dat)
    cat(paste0("read in data with ", nrows, " rows and ", nsamples, " samples\n"))

    return(dat)		     
}

#' Number of points to generate for a given feature
#'
#' Generate number of points to generate for a feature
#' @param value abundance in fraction from which to calculate the number of points
#' @return int
#' @examples
#' pointNumber(0.5)

pointNumber <- function(value){

    # out of a total of 5000 points
    n <- value*5000
    n <- round(n, 0)
    return(n)
}

#' Build points for each feature
#'
#' Generate a list of points for a given feature
#' @param taxon taxon name
#' @param n number of repititions as calculated using pointNumber()
#' @return a character vector
#' @examples
#' buildPointsPerTaxon("E.coli", 100)

buildPointsPerTaxon <- function(taxon, n){

    points.to.plot <- rep(taxon, n)
    return(points.to.plot)
}

#' Set the colours for ggplot
#'
#' Use rainbow() for generating colours to use for points
#' @param dat abundance data frame containing taxa for colouring
#' @return a character vector
#' @examples
#' abundances <- readAbundanceFile("test.tsv")
#' colours <- setColours(abundances)

setColours <- function(dat){

    # use rownames to set colours for
    # rainbow
    taxon.colours <- setNames(rainbow(nrow(dat)),
	                      rownames(dat))
    return(taxon.colours)
}

#' Subset a data frame by column number
#'
#' Subset a data frame for plotting each sample
#' @param dat abundance data frame
#' @param columnNumber column number to take
#' @return a data frame
#' @examples
#' subsetDataFrame(abundances, columnNumber=1)

subsetDataFrame <- function(dat, columnNumber=1){

    # keep rownames etc in subset of
    # data frame
    dat.subset <- data.frame(abundance = dat[,columnNumber])
    rownames(dat.subset) <- rownames(dat)
    return(dat.subset)
}

#' Build data frame with points that can be plotted over a sine wave
#'
#' Build data frame for downstream plotting
#' @param dat abundance data frame
#' @return a data frame
#' @examples
#' buildDataFrame(abundances)

buildDataFrame <- function(dat){

    taxa <- c()
	
    for (i in 1:nrow(dat)){
        v <- dat[i,]
	taxon <- rownames(dat)[i]
	n <- pointNumber(v)
	taxon.points <- buildPointsPerTaxon(taxon, n)
	taxa <- append(taxa, taxon.points)
    }
      
    # make sure its just 50000 points
    taxa <- taxa[1:5000]

    # randomoise taxa
    taxa <- sample(taxa)
	
    # create jittered wave
    t <- seq(0, 10, length=5000)
    t <- sin(t)

    # x variable
    x <- seq(1, 5000, 1)

    # randomise sizes of points
    # make it more likely to pick
    # smaller sizes

    pick.from <- c(rep(0.05, 75) , rep(0.2, 20), rep(1, 5))
    sizes <- sample(pick.from, 5000, replace=TRUE)

    df <- data.frame(taxon=taxa,
	             x=x,
                     y=t,
                     size=sizes)
    return(df)
}

#' Plot wave with or without epithelium using ggplot2
#'
#' Plot wave for each sample
#' @param dat data frame created with buildDataFrame
#' @param cols colours used for plotting (setColours())
#' @param g image file (epithelium)
#' @param with_image do you want to include the epithelium? Deafults to TRUE
#' @param background colour of background in the plots
#' @param text include text above the plot
#' @return a ggplot object (grob)
#' @examples
#' data(genus_abundances)
#' plotWave(genus_abundances, with_image=FALSE)
#' @import ggplot2

plotWave <- function(dat, cols, g="none", with_image=TRUE, background="white", text=""){

    p1 <- ggplot(dat, aes(x=x, y=jitter(y, amount=0.5), colour=taxon))

    if (with_image){
        if (g=="none"){
            stop("must provide image if with_image = TRUE")}
        p1 <- p1 + annotation_custom(g, xmin=1, xmax=5000, ymin=-3.1)}
    else{
	p1 <- p1}

    p2 <- p1 + geom_point(size=dat$size)
    p3 <- p2 + scale_color_manual(values=cols)
    p4 <- p3 + theme(panel.background=element_rect(fill=background, colour=background))
    p5 <- p4 + theme(panel.grid.major=element_line(colour=background))
    p6 <- p5 + theme(panel.grid.minor=element_line(colour=background))
    p7 <- p6 + theme(legend.position="none")
    p8 <- p7 + theme(axis.line=element_blank())
    p9 <- p8 + theme(axis.text.x=element_blank())
    p10 <- p9 + theme(axis.text.y=element_blank())
    p11 <- p10 + theme(axis.ticks=element_blank())
    p12 <- p11 + theme(axis.title.x=element_blank())
    p13 <- p12 + theme(axis.title.y=element_blank())
    p14 <- p13 + ylim(c(-1.8, 1.7))
    p15 <- p14 + annotate(geom="text", x=2500, y=1.7, family="Times", fontface="italic", size=6, label=text)
    return(p15)
}

#' Plot a wave based on microbial profiles
#'
#' Plot a wave representation of bacteria based on relative abundance data. Microwave
#' takes as input a data frame of microbial abundances and generates a set of points for
#' each taxon that are randomly asigned a position in a wave with an arbitrary colour and
#' size. This produces an image with an aesthetically pleasing view of an individual's
#' microbiome profile with the intended use for participant engagement. The aim is that
#' participants can receive an abstract view of their own microbiome. The wave can be produced
#' with or without a corresponding epithelial layer - maybe of use to folk that do not work
#' on gut samples but would like a visual representation of individual microbiomes.
#' @param abundances a data frame object 
#' @param with_image do you want to include the epithelium? Deafults to TRUE
#' @param background colour of background in the plots (colours supported in the
#' generic R colour palettes)
#' @param text include text above the plot. Defaults to no text.
#' @return .png files with plots per sample (columns in abundance data frame)
#' @examples
#' microwave(genus_abundances, with_image=FALSE)
#' @importFrom grid rasterGrob
#' @importFrom png readPNG
#' @export

microwave <- function(abundances, with_image=TRUE, animation=FALSE, background="white", text=""){

    nsamples <- ncol(abundances)
    img <- readPNG(system.file("data", "epithelium.png", package="micRowave"))
    g <- rasterGrob(img, interpolate=TRUE)

    dfs <- list()

    for (i in 1:ncol(abundances)){
        outfilename <- paste(colnames(abundances[i]), "png", sep=".")
	dat <- subsetDataFrame(abundances, columnNumber=i)
	colours <- setColours(dat)
	df <- buildDataFrame(dat)
        df$sample <- as.character(colnames(abundances)[i])
        dfs[[i]] <- df
	print(str(df))
        percent.done <- (i/nsamples)*100
	cat(paste0("creating plot ", i, "/", nsamples, " (", percent.done, "%)\n"))
	wave <- plotWave(df, colours, g=g, with_image=with_image, background=background, text=text)
	ggsave(outfilename, plot=wave, height=105, width=148, unit="mm")
    }
    if (animation==TRUE){
        library(gganimate)
	library(dplyr)
        to.animate <- bind_rows(dfs) 
        p <- plotWave(to.animate, colours, g=g, with_image=with_image, background=background, text=text)
	p <- p + transition_states(sample) + ease_aes("linear")
        animate(p, "allwaves.gif")
    }
}

