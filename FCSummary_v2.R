FCsummary <- function(fileName){
  
  ## Import data
  
  df <- readFC(fileName)
  
  ## Split data and for each combination of light condition, protein concentration and sample type, extract 
  ## mean fluorescence intensity in FL1 and FL6
  
  protConc <- c()
  sampleType <- c()
  lightCond <- c()
  FL1Median <- c()
  FL1SD <- c()
  FL3Median <- c()
  FL3SD <- c()
  FL6Median <- c()
  FL6SD <- c()
  GFPPercentage <- c()
  GFPSD <- c()
  
  df2 <- split(df, df$LightCond)
  
  for (i in df2){
    
    df3 <- split(i, i$ProteinConc)
    
    for (j in df3){
      
      df4 <- split(j, j$SampleType)
      
      for (k in df4){
        
        protConc <- append(protConc, toString(k[1,2]))
        sampleType <- append(sampleType, toString(k[1,3]))
        lightCond <- append(lightCond, toString(k[1,4]))
        FL1Median <- append(FL1Median, toString(round(mean(k$FL1Median), digits=2)))
        FL1SD <- append(FL1SD, toString(round(sd(k$FL1Median), digits=2)))
        FL3Median <- append(FL3Median, toString(round(mean(k$FL3Median), digits=2)))
        FL3SD <- append(FL3SD, toString(round(sd(k$FL3Median), digits=2)))
        FL6Median <- append(FL6Median, toString(round(mean(k$FL6Median), digits=2)))
        FL6SD <- append(FL6SD, toString(round(sd(k$FL6Median), digits=2)))
        GFPPercentage <- append(GFPPercentage, toString(round(mean(k$GFP_percentage), digits=2)))
        GFPSD <- append(GFPSD, toString(round(sd(k$GFP_percentage), digits=2)))
        
      }
    }
  }
  
  ## Put dataframe together
  
  df_final <- data.frame(LightCond=lightCond, ProteinConc=protConc, SampleType=sampleType, FL1Median=FL1Median, 
                         FL1SD=FL1SD, FL3Median=FL3Median, FL3SD=FL3SD, FL6Median=FL6Median, FL6SD=FL6SD, GFPPercentage=GFPPercentage, GFPSD=GFPSD, stringsAsFactors = FALSE)
  
  ## Write output file
  
  write.table(df_final, file=paste(substr(fileName, 1, nchar(fileName)-4), "_summary.csv", sep=""), sep=";", dec=".", col.names=NA) 
  
  ## Return value
  
  df_final
}