readFC <- function(fileName){
  
  ## Import packages
  
  library(stringr)
  
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
    
    if (grepl("cells", i, fixed=TRUE) | grepl("blank", i, fixed=TRUE)){
      
      conc <- 0
      sample <- "blank"
      
    }else if (grepl("aav", i, fixed=TRUE)){
      
      conc <- 0
      sample <- "AAV"
      
    }else{
      
      nameSplit <- strsplit(i, "_")
      
      conc <- nameSplit[[1]][2]
      conc <- str_sub(conc, str_length(conc), str_length(conc))
      
      sample <- nameSplit[[1]][3]
      
      ## Extract AAV type
      
    }
    
    if (sample == 2 | sample == 4 | sample == "AAV"){
      
      if (grepl("VP2-PIF6", i, fixed=TRUE)){
        sample <- paste(sample, "VP2-PIF6")
        
      }else if (grepl("pMH301", i, fixed=TRUE)){
        sample <- paste(sample, "pMH301")
      }
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

  ## Create vector with FL6 median
  
  FL6Median <- c()
  for (i in seq(from = FL6FirstRow, to = dim(df)[1], by = rowsPerSample)){
    FL6Median <- append(FL6Median, as.numeric(df[i,4]))
  }
  
  ## Create vector with % GFP+
  
  GFPperc <- c()
  for (i in seq(from = GFPFirstRow, to = dim(df)[1], by = rowsPerSample)){
    GFPperc <- append(GFPperc, as.numeric(df[i,3]))
  }

  ## Create dataframe and save it as csv
  
  df_final <- data.frame(Sample = sampleName, ProteinConc = protConc, SampleType = sampleType, 
                  LightCond = lightCond, FL1Median = FL1Median, FL6Median = FL6Median, 
                  GFP_percentage = GFPperc, stringsAsFactors = FALSE)
  
  ## Remove rows with NAs
  
  df_final <- df_final[complete.cases(df_final[ , 5:7]),]
  
  ## Sort dataframe
  
  arrange(df_final, LightCond, ProteinConc)
  
  ## Write output file
  
  write.table(df_final, file=paste(substr(fileName, 1, nchar(fileName)-4), "_re-formated.csv", sep=""), sep=";", 
              dec=".", col.names=NA) 
  
  ## Return value
  
  df_final
}
