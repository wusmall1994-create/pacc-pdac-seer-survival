args<-commandArgs(trailingOnly=TRUE)
project<-if(length(args)) args[1] else '.'
outfile<-if(length(args)>1) args[2] else file.path(project,'04_results/revision/conditional_OS/riskset_composition.csv')
dir.create(dirname(outfile),recursive=TRUE,showWarnings=FALSE)
d<-readRDS(file.path(project,'02_processed_data/pacc_pdac_analysis_cohort.rds'))
stopifnot(nrow(d)==116516,sum(d$histology=='pACC')==580)
rows<-list();k<-0
for(cutoff in c(2023,2015,2013)) for(lm in c(0,1,2,3,5)) {
 z<-d[d$histology=='pACC' & d$year<=cutoff & d$os_time>lm,]
 qage<-quantile(z$age,c(.25,.5,.75),names=FALSE)
 qyear<-quantile(z$year,c(.25,.5,.75),names=FALSE)
 k<-k+1
 rows[[k]]<-data.frame(diagnosis_end=cutoff,landmark=lm,n=nrow(z),age_q1=qage[1],age_median=qage[2],age_q3=qage[3],year_q1=qyear[1],year_median=qyear[2],year_q3=qyear[3],localized=sum(z$stage=='Localized'),regional=sum(z$stage=='Regional'),distant=sum(z$stage=='Distant'),unknown=sum(z$stage=='Unknown'))
}
x<-do.call(rbind,rows)
stopifnot(all(rowSums(x[,c('localized','regional','distant','unknown')])==x$n))
write.csv(x,outfile,row.names=FALSE)
print(x)
