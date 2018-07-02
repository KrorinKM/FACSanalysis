readFC <- function(fileName){
  
  ## fileName <- "C:/Users/Carolina/Google Drive (cjlongres@gmail.com)/Carolina Master thesis/Exp007 - Flow cytometry 2/CJ Exp7 Statistics.csv"
  ## fileName <- "C:/Users/Carolina/Google Drive (cjlongres@gmail.com)/Carolina Master thesis/Exp006 - Protein expression/CJ Exp6 FC data/CJ Exp6 Statistics.csv"
  
  ## Import packages
  
  library(stringr)
  
  ## Read csv data
  
  df <- read.table(fileName, header=TRUE, sep=";", dec=",")
  df <- df[,c("Data.Set", "X.Parameter", "X.Gated", "X.Med")]
  
  ## Create vector with sample names
  
  sampleName <- c()
  for (i in seq(from = 1, to = dim(df)[1], by = 12)){
    sampleName <- append(sampleName, toString(df[i,1]))
  }
  
  ## Create vectors with info extracted from sample names
  
  lightCond <- c()
  protConc <- c()
  sampleType <- c()
  aavType <- c()
  
  for (i in sampleName){
    
  ## Extract protein concentration and sample type
    
    if (grepl("cells", i, fixed=TRUE) | grepl("blank", i, fixed=TRUE)){
      
      conc <- 0
      type <- "blank"
      
    }else if (grepl("aav", i, fixed=TRUE)){
      
      conc <- 0
      type <- "AAV"
      
    }else{
      
      nameSplit <- strsplit(i, "_")
      
      conc <- nameSplit[[1]][2]
      type <- nameSplit[[1]][2]
      
      conc <- str_sub(conc, str_length(conc), str_length(conc))
      type <- str_sub(type, 1, 1)
    }
    
    protConc <- append(protConc, conc)
    sampleType <- append(sampleType, type)
    
  ## Extract illumination conditions  
    
    if (grepl("660-dark", i, fixed=TRUE)){
      wl <- "660-dark"
    }else if (grepl("740", i, fixed=TRUE)){
      wl <- "740"
    }else{
      wl <- "660"
    }
      
    lightCond <- append(lightCond, wl)
  
  ## Extract AAV type
    
    if (grepl("VP2-PIF6", i, fixed=TRUE)){
      aav <- "VP2-PIF6"
    }else if (grepl("pMH301", i, fixed=TRUE)){
      aav <- "pMH301"
    }else{
      aav <- NA
    } 
    
    aavType <- append(aavType, aav)
    
  }
  
  ## Create vector with FL1 median
  
  FL1Median <- c()
  for (i in seq(from = 1, to = dim(df)[1], by = 6)){
    FL1Median <- append(FL1Median, df[i,3])
  }

  ## Create vector with FL6 median
  
  FL6Median <- c()
  for (i in seq(from = 4, to = dim(df)[1], by = 6)){
    FL6Median <- append(FL6Median, df[i,3])
  }

  ## Create dataframe and save it as csv
  
  df_final <- data.frame(Sample=sampleName, ProteinConc=protConc, SampleType=sampleType, LightCond=lightCond, FL1Median=FL1Median, 
                   FL6Median=FL6Median, stringsAsFactors = FALSE)
  
  ## Write output file
  
  write.table(df_final, file=paste(substr(fileName, 1, nchar(fileName)-4), "_re-formated.csv", sep=""), sep=";", 
              dec=".", col.names=NA) 
  
  ## Return value
  
  df_final
}
