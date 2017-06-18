load('../R_MAIN_PAPER_4.0/Fi_gdsc1000_DATA/DRUG SCREENING/R/simplyIC50s-Final-RMSE-Filtered-0.3.rdata')
load('../R_MAIN_PAPER_4.0/Fi_gdsc1000_DATA/SEQUENCING/CELL LINES/R/BEMs/PANCAN_SEQ_BEM.rdata')

BEM<-NGS_BEM$logical
commonSamples<-intersect(rownames(IC50s),colnames(BEM))
commonSamples<-commonSamples[which(!is.na(IC50s[commonSamples,'34']))]

DRUGresponseP<-IC50s[commonSamples,'34']
#34 is the drug id of Imatinib

MUTpattern<-BEM['BCR-ABL',commonSamples]

Tres<-t.test(DRUGresponseP~MUTpattern,var.equal = TRUE)
fit <- aov(DRUGresponseP~MUTpattern)
Y<-anova(fit)

print(Y[1,5])
print(Tres$p.value)