########################################################################
#                                                                      #
#     ANOVA_pipeline: Script performing ANOVA                     	   #	
#					 correlating drug response to cell line features,  #
#					 for the GDSC1000 project, including               #
#                    screen medium                                     #
#                    as co-factors                                     #
#                             v0.0.4 - 20151124                        # 
#                                                                      #
# Copyright (c) 2014 - 2019, EMBL - European Bioinformatics Institute  #
#       Author: Francesco Iorio (iorio@ebi.ac.uk)                      #
#       Distributed under the GPLv3 License.                           #
########################################################################
source('Fi_gdsc1000_CODE/Project_Preamble.R')
source('Fi_gdsc1000_CODE/MISC/my.compress_identical_patterns.R')

DATE<-Sys.Date()
TIME<-Sys.time()
SYSINFO<-Sys.info()

cat(paste('\n\n**** gdscANOVA analysis - ',TIME,' ****',sep=''))
cat(paste('\nExecuted by ',SYSINFO[7],' on ',SYSINFO[4],sep=''))
cat('\n*********************************************************')

cat('\n\n- Loading gdscANOVA packages...')
source('Fi_gdsc1000_CODE/gdscANOVA/CELL LINES MULTIFACTOR/gdscANOVA_packages.R')
cat('\n+ Done!\n\n')

cat('- Loading gdscANOVA R objects...')
source('Fi_gdsc1000_CODE/gdscANOVA/CELL LINES MULTIFACTOR//gdscANOVA_preamble.R')
cat('\n+ Done!\n\n')

cat('\n\n- Loading gdscANOVA functions and objects...')
source('Fi_gdsc1000_CODE/gdscANOVA/CELL LINES MULTIFACTOR//gdscANOVA_library.R')
cat('\n+ Done!\n\n')

TIME<-str_replace_all(TIME,'[:]','-')
TIME<-str_replace_all(TIME,'[ ]','_')

MAIN_DIR<-paste(gdscANOVA.results.dir,'/',
                gdscANOVA.settings.CELL_LINES,'_',
                gdscANOVA.settings.DRUG_domain,'_',
                TIME,'/',sep='')
INPUT_DIR<-paste(MAIN_DIR,'INPUT/',sep='')
OUTPUT_DIR<-paste(MAIN_DIR,'OUTPUT/',sep='')

if(!exists(MAIN_DIR)){
  dir.create(MAIN_DIR)
}

if(!exists(INPUT_DIR)){
  dir.create(INPUT_DIR)
}

if(!exists(OUTPUT_DIR)){
  dir.create(OUTPUT_DIR)
}


cat(paste('- Retrieving data for the selected drug domain (',gdscANOVA.settings.DRUG_domain,')...',sep=''))
cat(paste('\n\t\t'))
IC50s<-gdscANOVA_createDrugDataInput(gdscANOVA.settings.DRUG_domain)
cat(paste('\n\t\t(',gdscANOVA.settings.SCREENING_VERSION,')',sep=''))
cat('\n+ Done!\n\n')

cat('- Assembling input features')
cat(paste('\n\t\t'))
InputFeatures<-gdscANOVA_createInputFeatures(additional_features=gdscANOVA.additionalFeatures,
                                             oneGeneOnly = gescANOVA.settings.oneGeneOnly,
                                             excludeHyperMetData = gdscANOVA.settings.excludeHyperMetData,
                                             additional_features_only = gdsdANOVA.settings.additionalFeaturesOnly)
cat('\n+ Done!\n\n')

save(InputFeatures,file=paste(INPUT_DIR,'InputFeatures.rdata',sep=''))
save(IC50s,file=paste(INPUT_DIR,'IC50s.rdata',sep=''))

DIAGNOSTICS<-gdscANOVA_diagnostics()

gdscANOVA_create_systemInfos()

TOTRES<-gdscANOVA_totalANOVA(fn=paste(OUTPUT_DIR,'ANOVA_results',sep=''))
 
write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(OUTPUT_DIR,'ANOVA_results','.txt',sep=''))
save(TOTRES,file=paste(OUTPUT_DIR,'ANOVA_results','.rdata',sep=''))
 
source('Fi_gdsc1000_CODE/gdscANOVA/CELL LINES MULTIFACTOR/gdscANOVA_graphics.R')
source('Fi_gdsc1000_CODE/gdscANOVA/CELL LINES MULTIFACTOR/gdscANOVA_ind_pack_creation.R')

