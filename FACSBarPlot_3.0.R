FACSplot <- function(dataset){
  
  ### dataset <- read.table("CJ exp9 Statistics_summary.csv", header=TRUE, sep=";", dec=",", stringsAsFactors=FALSE)
  ### dataset <- dataset[,2:10]
  
  # The dataset argument is a dataframe produced by the FACSSummary function
  
  ## Import packages
  
  library(ggplot2)
  library(grid)
  library(gridExtra)
  
  ## Choose conditions to be plotted (wavelength and protein concentration)
  
  wl <- "740"
  conc <- "s" 
  
  dataset <- dataset[which(dataset$LightCond == wl & (dataset$ProteinConc == "0" | dataset$ProteinConc == conc)),]
  
  ## Split dataframe into mean and error dataframes
  
  means <- dataset[,which(grepl("Median", colnames(dataset), fixed=TRUE) | grepl("Percentage", colnames(dataset), fixed=TRUE))]
  errors <- dataset[,which(grepl("SD", colnames(dataset), fixed=TRUE))]
  
  ## Convert data into dataframe usable by geom_bar
  
  x <- c(rep(dataset$SampleType, dim(means)[2]))   ### Sample type
  y <- c(means$FL1Median, means$FL3Median, means$FL6Median, means$GFPPercentage)
  groups <- c(rep(colnames(means), each=dim(means)[1]))      ### Type of data
  error <- c(errors$FL1SD, errors$FL3SD, errors$FL6SD, errors$GFPSD)
  
  df <- data.frame(x, as.numeric(y), groups, as.numeric(error), stringsAsFactors = TRUE)
  df <- setNames(df, c("x","y","groups", "error"))
  
  ## Separate df into medians and GFP+ percentage
  
  dfMedians <- df[which(df$groups == "FL1Median"  | df$groups == "FL3Median" | df$groups == "FL6Median"),]
  dfGFP <- df[which(df$groups == "GFPPercentage"),]
  
  ## Bar plot
  
  p <- c()
  ### Exp 9.1
  ### xlimits <- c("blank", "AAV pMH301", "AAV VP2-PIF6", "1", "2 pMH301", "2 VP2-PIF6", "3", "4 pMH301", "4 VP2-PIF6")
  ### xlabels <- c("Blank", "AAV pMH301", "AAV VP2-PIF6", "PhyB", "PhyB + pMH301", "PhyB + VP2-PIF6", "PhyB-DARPin", "PhyB-DARPin + pMH301", "PhyB-DARPin + VP2-PIF6")
  
  ### Exp 9.2 & 9.3
  ### xlimits <- c("blank", "1", "2", "3", "4", "5")
  ### xlabels <- c("Blank", "PhyB", "PhyB-DARPin", "PhyB + mVenus-PIF6", "PhyB-DARPin + mVenus-PIF6", "mVenus-PIF6")
  
  ### Exp 12
  xlimits <- c("blank", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
  xlabels <- c("Blank", "PhyB-DARPin c1", "AAV 303", "AAV 303 + PhyB-DARPin c1", "AAV 303 + PhyB-DARPin c2", "AAV 303 + PhyB-DARPin c3", "AAV 303 + PhyB-DARPin c4", 
                "AAV 303 + PhyB-DARPin -TCEP", "AAV 156 + FMD (ctrl)", "AAV 156 + FMD + AP21967 (ctrl)", "AAV 156 + FMD", "AAV 156 + FMD + AP21967")
  
  ymax <- 4
  
  p$Medians <- ggplot(data=dfMedians, aes(x=x,y=y, fill=groups)) + 
    geom_bar(stat="identity", position=position_dodge(), width=.8) + 
    geom_errorbar(aes(ymin = y - error, ymax = y + error, color = "red"), width=.1, position=position_dodge(.9), show.legend=FALSE) +
    scale_fill_manual(values = c("#D8D8D8", "#A4A4A4", "#848484")) + theme_classic() +
    labs(y="Median of fluorescence intensity", x="", fill="", subtitle="", title="", caption="") +
    scale_x_discrete(limits = xlimits, labels = xlabels) +
    theme(axis.text.x = element_text(hjust=1, vjust=0.9, angle=60), axis.title.x = element_text(vjust=-5), axis.title.y = element_text(vjust=2), legend.position = c(0.15,0.95)) + 
    scale_y_continuous(breaks=seq(0,ymax,0.2), limits=c(0,ymax))
  
  if (dataset$GFPPercentage[1] != "NA"){
    p$GFP <- ggplot(data=dfGFP, aes(x=x,y=y, fill=groups)) + 
      geom_bar(stat="identity", position=position_dodge(), width=.4,  show.legend=FALSE) + 
      geom_errorbar(aes(ymin = y - error, ymax = y + error, color = "red"), width=.1, position=position_dodge(.9), show.legend=FALSE) +
      scale_fill_manual(values = c("#D8D8D8", "#A4A4A4", "#848484")) + theme_classic() +
      labs(y="Percentage GFP+ cells", x = "", fill="", subtitle="", title="", caption="") +
      scale_x_discrete(limits = xlimits, labels = xlabels) +
      theme(axis.text.x = element_text(hjust=1, vjust=0.9, angle=60), axis.title.x = element_text(vjust=-5), axis.title.y = element_text(vjust=2)) +
      scale_y_continuous(breaks=seq(0,50,1), limits=c(0,7))
  } 
  
  ## Arrange plots into one grid
  
  if (dataset$GFPPercentage[1] != "NA"){
    grid.arrange(p$Medians, p$GFP, nrow = 1)
  }else{
    p$Medians
  }
  
}

