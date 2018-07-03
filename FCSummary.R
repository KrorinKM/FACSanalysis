FCsummary <- function(fileName){
  
  ## Import data
  
  df <- readFC(fileName)
  
  ## Split data and for each combination of light condition, protein concentration and sample type, extract 
  ## mean fluorescence intensity in FL1 and FL6
  
  protConc <- c()
  sampleType <- c()
  lightCond <- c()
  FL1Mean <- c()
  FL1SD <- c()
  FL6Mean <- c()
  FL6SD <- c()
  
  df2 <- split(df, df$LightCond)
  
  for (i in df2){
    
    df3 <- split(i, i$ProteinConc)
    
    for (j in df3){
      
      df4 <- split(j, j$SampleType)
      
      for (k in df4){
      
        protConc <- append(protConc, toString(k[1,2]))
        sampleType <- append(sampleType, toString(k[1,3]))
        lightCond <- append(lightCond, toString(k[1,4]))
        FL1Mean <- append(FL1Mean, toString(round(mean(k$FL1Median), digits=2)))
        FL1SD <- append(FL1SD, toString(round(sd(k$FL1Median), digits=2)))
        FL6Mean <- append(FL6Mean, toString(round(mean(k$FL6Median), digits=2)))
        FL6SD <- append(FL6SD, toString(round(sd(k$FL6Median), digits=2)))
        
      }
    }
  }
  
  ## Put dataframe together
  
  df_final <- data.frame(LightCond=lightCond, ProteinConc=protConc, SampleType=sampleType, FL1Mean=FL1Mean, 
                         FL1SD=FL1SD, FL6Mean=FL6Mean, FL6SD=FL6SD, stringsAsFactors = FALSE)
  
  ## Write output file
  
  write.table(df_final, file=paste(substr(fileName, 1, nchar(fileName)-4), "_summary.csv", sep=""), sep=";", 
              dec=".", col.names=NA) 
  
  ## Return value
  
  df_final
}