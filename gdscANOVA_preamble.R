

## ANOVA paths and filenames

annotations.master_list.fn<-'../R_MAIN_PAPER_5.0/Fi_gdsc1000_DATA/ANNOTATIONS/CELL LINES/R/MASTER_LIST_EXTENDED_INFOS_Mar2016.Rdata'

gdscANOVA.screeMedium.fn<-'Fi_gdsc1000_DATA/DRUG SCREENING/R/ASSAY_PROPERTIES.rdata'
gdscANOVA.screening.fn<-'../R_MAIN_PAPER_5.0/Fi_gdsc1000_DATA/DRUG SCREENING/R/v17a_IC50s.Rdata'

#gdscANOVA.screening.DRcurves.fn<-'Fi_gdsc1000_DATA/DRUG SCREENING/R/DoseResponseCurves.rdata'
gdscANOVA.additionalFeatures.fn<-NULL
gdscANOVA.results.dir<-'../R_MAIN_PAPER_5.0/Fi_gdsc1000_RESULTS/gdscANOVA/MultiFactorial/'
gdscANOVA.cellLineFeatures.dir<-'Fi_gdsc1000_DATA/MultiOmics/CELL LINES/R/'

gdscANOVA.maxTestedConc.fn<-'Fi_gdsc1000_DATA/DRUG SCREENING/R/maximalTestedConcentrations.rdata'
gdscANOVA.drugProperties.fn<-'../R_MAIN_PAPER_5.0/Fi_gdsc1000_DATA/ANNOTATIONS/DRUGS/DRUG_ANALYSIS_SET_20160324.Rdata'


gdscANOVA.drugOwnership.fn<-'Fi_gdsc1000_DATA/ANNOTATIONS/DRUGS/DRUG_BY_COMPANIES_20150126.rdata'


## ANOVA file name suffixes
gdscANOVA.cellLineFeatures.fn_suffix<-'_simple_MOBEM.rdata'

## ANOVA settings
gdscANOVA.settings.SCREENING_VERSION<-paste('v17a - Mar2016',Sys.Date())
gdscANOVA.settings.DRUG_domain<-"GDSC1000_paper_set"
gdscANOVA.settings.CELL_LINES<-'PANCAN'
gdscANOVA.settings.includeMEDIA_factor<-TRUE
gdscANOVA.settings.includeMSI_Factor<-TRUE
gdscANOVA.settings.analysisType<-'PANCAN'
gdscANOVA.settings.additionalFeatures<-FALSE
gdsdANOVA.settings.additionalFeaturesOnly<-FALSE
gdscANOVA.settings.excludeHyperMetData<-TRUE
gescANOVA.settings.oneGeneOnly<-NULL
gdscANOVA.settings.resPackageHeader<-paste('ANOVA',gdscANOVA.settings.CELL_LINES,'result package')
gdscANOVA.settings.resPackageIncludeSangerLogo<-TRUE
gdscANOVA.settings.resPackageIncludeNKILogo<-TRUE
gdscANOVA.settings.resPackageIncludeEBIlogo<-TRUE
gdscANOVA.settings.featFactorPopulationTh<-3
gdscANOVA.settings.MSIfactorPopulationTh<-2

gdscANOVA.settings.pval_correction_method<-'qvalue'

if (gdscANOVA.settings.CELL_LINES == 'PANCAN'){
  gdscANOVA.settings.pval_TH<-10^-3  
}else{
  gdscANOVA.settings.pval_TH<-10^-3
}

gdscANOVA.settings.FDR_TH<-25

load(gdscANOVA.screening.fn)
load(paste(gdscANOVA.cellLineFeatures.dir,gdscANOVA.settings.CELL_LINES,gdscANOVA.cellLineFeatures.fn_suffix,sep=''))

TOTALBEM <-MoBEM

load(gdscANOVA.drugProperties.fn)
DRUG_PROPS<-
    DRUG_ANALYSIS_SET
load(gdscANOVA.drugOwnership.fn)
load(gdscANOVA.maxTestedConc.fn)
#load(gdscANOVA.screening.DRcurves.fn)

load(annotations.master_list.fn)


if(gdscANOVA.settings.additionalFeatures){
  print('Loading Additional Features...')
  load(gdscANOVA.additionalFeatures.fn)
  gdscANOVA.additionalFeatures=BEM
  print('+ Done!')
}else{
  gdscANOVA.additionalFeatures=NULL
}

