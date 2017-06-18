load('Fi_gdsc1000_RESULTS/gdscANOVA/PANCAN_GDSC1000_newModel_new_data/OUTPUT/ANOVA_results.rdata')
nnTOTRES<-TOTRES

load('Fi_gdsc1000_RESULTS/gdscANOVA/PANCAN_GDSC1000_newModel_old_data/OUTPUT/ANOVA_results.rdata')
noTOTRES<-TOTRES

# load('Fi_gdsc1000_RESULTS/gdscANOVA/PANCAN_GDSC1000_oldModel_old_data/OUTPUT/MANOVA_results.rdata')
# ooTOTRES<-TOTRES

commons<-rep(0,500)
no_only<-rep(0,500)
nn_only<-rep(0,500)

for (i in 1:500){
  print(i)
  commons[i]<-length(intersect(paste(noTOTRES[1:i,2],noTOTRES[1:i,3],noTOTRES[1:i,4]),paste(nnTOTRES[1:i,2],nnTOTRES[1:i,3],nnTOTRES[1:i,4])))
  no_only[i]<-length(setdiff(paste(noTOTRES[1:i,2],noTOTRES[1:i,3],noTOTRES[1:i,4]),paste(nnTOTRES[1:i,2],nnTOTRES[1:i,3],nnTOTRES[1:i,4])))
  nn_only[i]<-length(setdiff(paste(nnTOTRES[1:i,2],nnTOTRES[1:i,3],nnTOTRES[1:i,4]),paste(noTOTRES[1:i,2],noTOTRES[1:i,3],noTOTRES[1:i,4])))
}

plot(100*commons[1:300]/(1:300),type='b',ylim=c(0,100),pch=16,col='purple',ylab='',xlab='top n.significant hits')
par(new=TRUE)
plot(as.numeric(noTOTRES[1:300,"ANOVA FEATURE FDR %"]),type='l',ylim=c(0,100),pch=16,col='red',ylab='',xlab='',xaxt='n',yaxt='n',frame.plot = FALSE)
par(new=TRUE)
plot(as.numeric(nnTOTRES[1:300,"ANOVA FEATURE FDR %"]),type='l',ylim=c(0,100),pch=16,col='blue',ylab='',xlab='',xaxt='n',yaxt='n',frame.plot = FALSE)





