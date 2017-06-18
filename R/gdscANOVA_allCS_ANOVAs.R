
load('Fi_gdsc1000_DATA/gdscANOVA/R/CancerTypes.rdata')
load('Fi_gdsc1000_DATA/gdscANOVA/R/IncludeMSIfactor.rdata')

nCancerTypes<-length(CancerTypes)

for (i in 1:nCancerTypes){
  print(paste('- Executing ',CancerTypes[i],'-specific ANOVA',sep=''))
  
  gdscANOVA.settings.CELL_LINES<-CancerTypes[i]
  gdscANOVA.settings.includeMSI_Factor<-IncludeMSIfactor[i]
  
  source('Fi_gdsc1000_CODE/gdscANOVA/CELL LINES MULTIFACTOR/gdscANOVA_pipeline.R')
  
  vars<-ls()
  vars<-setdiff(vars,c('CancerTypes','IncludeMSIfactor','nCancerTypes','i'))
  rm(list=vars)
  
  print('+ DONE!')
}

