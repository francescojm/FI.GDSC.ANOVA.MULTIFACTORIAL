
load('Fi_gdsc1000_DATA/gdscANOVA/R/CancerTypes.rdata')
load('Fi_gdsc1000_DATA/gdscANOVA/R/IncludeMSIfactor.rdata')

gdscANOVA.settings.CELL_LINES<-'LUAD'
gdscANOVA.settings.analysisType<-'CS'
gdscANOVA.settings.includeMSI_Factor<-FALSE
gdscANOVA.settings.includeMEDIA_factor<-FALSE
source('Fi_gdsc1000_CODE/gdscANOVA/CELL LINES MULTIFACTOR/gdscANOVA_pipeline.R')
load('Fi_gdsc1000_DATA/ANNOTATIONS/DRUGS/DRUG_ANALYSIS_SET_20150123.rdata')

DRUG_ID<-'1019'

gdscANOVA_individualANOVA(DRUG_ID = DRUG_ID,"EGFR_mut",display = TRUE)
