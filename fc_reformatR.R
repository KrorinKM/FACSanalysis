readFC <- function(fileName){
  
  ## Import packages
  
  library(stringr)
  library(plyr)
  library(dplyr)
  
  ## Read csv data
  
  df <- read.table(fileName, header=TRUE, sep=";", dec=",", stringsAsFactors=FALSE)
  df <- df[,c("Data.Set", "X.Parameter", "X.Gated", "X.Med")]
  df$X.Gated <- sub(",",".", df$X.Gated)
  df$X.Med <- sub(",",".", df$X.Med)
  
  ## Specify csv file parameters -> Change for every dataset
  
  rowsPerSample <- 11
  FL1FirstRow <- 5      ### Row at which FL1 Median appears for the first time
  FL6FirstRow <- 9      ### Row at which FL6 Median appears for the first time
  GFPFirstRow <- 6      ### Row at which the percentage of GFP+ cells appears for the first time
  FL3FirstRow <- 8
  
  ## Create vector with sample names
  
  sampleName <- c()
  for (i in seq(from = 1, to = dim(df)[1], by = rowsPerSample)){
    sampleName <- append(sampleName, toString(df[i,1]))
  }
  
  ## Create vectors with info extracted from sample names
  
  lightCond <- c()
  protConc <- c()
  sampleType <- c()
  
  for (i in sampleName){
    
  ## Extract protein concentration and sample type
    
    if (grepl("blank", i, fixed=TRUE)){
      
      conc <- 0
      sample <- "blank"
      
    }else{
      
      nameSplit <- strsplit(i, "_")
      
      conc <- nameSplit[[1]][2]
      conc <- str_sub(conc, str_length(conc), str_length(conc))
      
      sample <- nameSplit[[1]][3]
      sample <- strsplit(sample, " ")[[1]][1]
      
    }
    
    protConc <- append(protConc, conc)
    sampleType <- append(sampleType, sample)
    
  ## Extract illumination conditions  
    
    if (grepl("660-dark", i, fixed=TRUE)){
      wl <- "660-dark"
    }else if (grepl("740", i, fixed=TRUE)){
      wl <- "740"
    }else{
      wl <- "660"
    }
      
    lightCond <- append(lightCond, wl)
    
  }
  
  ## Create vector with FL1 median
  
  FL1Median <- c()
  for (i in seq(from = FL1FirstRow, to = dim(df)[1], by = rowsPerSample)){
    FL1Median <- append(FL1Median, as.numeric(df[i,4]))
  }

  ## Create vector with FL3 median
  
  FL3Median <- c()
  for (i in seq(from = FL3FirstRow, to = dim(df)[1], by = rowsPerSample)){
    FL3Median <- append(FL3Median, as.numeric(df[i,4]))
  }
  
  ## Create vector with FL6 median
  
  FL6Median <- c()
  for (i in seq(from = FL6FirstRow, to = dim(df)[1], by = rowsPerSample)){
    FL6Median <- append(FL6Median, as.numeric(df[i,4]))
  }
  
  ## Create vector with % GFP+
  
  GFPperc <- c()
  
  if (GFPFirstRow != 0){
    
    for (i in seq(from = GFPFirstRow, to = dim(df)[1], by = rowsPerSample)){
      GFPperc <- append(GFPperc, as.numeric(df[i,3]))
    }
  
  }else{
    
    GFPperc <- rep(NA, length(sampleName))
    
  }

  ## Create dataframe and save it as csv
  
  df_final <- data.frame(Sample = sampleName, ProteinConc = protConc, SampleType = sampleType, 
                  LightCond = lightCond, FL1Median = FL1Median, FL3Median, FL6Median = FL6Median, 
                  GFP_percentage = GFPperc, stringsAsFactors = FALSE)
  
  ## Remove rows with NAs
  
  df_final <- df_final[complete.cases(df_final[ , 5:6]),]
  
  ## Sort dataframe
  
  df_final <- arrange(df_final, Sample)
  
  ## Write output file
  
  write.table(df_final, file=paste(substr(fileName, 1, nchar(fileName)-4), "_re-formated.csv", sep=""), sep=";", dec=".", col.names=NA) 
  
  ## Return value
  
  df_final
}
