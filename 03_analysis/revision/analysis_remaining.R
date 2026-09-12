options(stringsAsFactors=FALSE,width=180)
suppressPackageStartupMessages(library(survival))
args<-commandArgs(trailingOnly=TRUE)
root <- if(length(args)>=1)args[1] else '.'
out <- if(length(args)>=2)args[2] else file.path(root,'04_results/revision/selection_and_calendar')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
savecsv<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='NA')
raw<-read.delim(file.path(root,'01_data_extraction/pacc_pdac_seer17_2004_2023_raw.txt'),sep='\t',quote='"',check.names=FALSE,colClasses='character',na.strings='NA')
raw<-raw[raw[['Histologic Type ICD-O-3']] %in% c('8140','8500','8550'),]
raw$histology<-ifelse(raw[['Histologic Type ICD-O-3']]=='8550','pACC','PDAC')
raw$confirmed<-raw[['Diagnostic Confirmation']] %in% c('Positive histology','Pos hist AND immunophenotyping AND/OR pos genetic studies')
raw$other_eligible<-!raw[['Type of Reporting Source']] %in% c('Autopsy only','Death certificate only') & !is.na(raw[['Survival Days']]) & !is.na(raw[['Vital status recode (study cutoff used)']])
raw$age<-as.numeric(sub('[^0-9].*$','',raw[['Age recode with single ages and 85+']]))
raw$year<-as.numeric(raw[['Year of diagnosis']])
sr<-raw[['Combined Summary Stage with Expanded Regional Codes (2004+)']]
raw$stage<-ifelse(grepl('^Localized',sr),'Localized',ifelse(grepl('^Regional',sr),'Regional',ifelse(grepl('^Distant',sr),'Distant','Unknown')))
raw$era<-ifelse(raw$year<=2009,'2004-2009',ifelse(raw$year<=2015,'2010-2015',ifelse(raw$year<=2017,'2016-2017','2018-2023')))
raw$sex<-raw[['Sex']]
raw$site<-ifelse(grepl('^C25.0',raw[['Primary Site - labeled']]),'Head',ifelse(grepl('^C25.[12]',raw[['Primary Site - labeled']]),'Body/tail','Other/unspecified'))
sc<-ifelse(raw$year<=2022,raw[['RX Summ--Surg Prim Site (1998-2022)']],raw[['RX Summ--Surg Prim Site 2023 (2023+)']])
raw$surgery<-ifelse(is.na(sc)|sc %in% c('Blank(s)','00','A000','99','A990'),'No/unknown','Yes')
raw$chemotherapy<-ifelse(raw[['Chemotherapy recode (yes, no/unk)']]=='Yes','Yes','No/unknown')
raw$confirmation<-raw[['Diagnostic Confirmation']]
flow<-do.call(rbind,lapply(c('PDAC','pACC'),function(h){z<-raw[raw$histology==h,]; data.frame(histology=h,target=nrow(z),confirmation_pass=sum(z$confirmed),confirmation_excluded=sum(!z$confirmed),later_excluded=sum(z$confirmed & !z$other_eligible),final=sum(z$confirmed & z$other_eligible),excluded_confirmation_only_otherwise_eligible=sum(!z$confirmed & z$other_eligible))}))
savecsv(flow,'selection_flow')
smd_binary<-function(p,q){den<-sqrt((p*(1-p)+q*(1-q))/2);if(den==0) if(p==q)0 else NA_real_ else abs(p-q)/den}
desc<-list();k<-0
for(universe in c('confirmation_step','otherwise_eligible')) for(h in c('PDAC','pACC')){
 z<-raw[raw$histology==h & (universe=='confirmation_step' | raw$other_eligible),]
 a<-z[z$confirmed,];b<-z[!z$confirmed,]
 for(v in c('age','year')){
  smd<-abs(mean(a[[v]])-mean(b[[v]]))/sqrt((var(a[[v]])+var(b[[v]]))/2)
  fmt<-function(x)sprintf('%.0f (%.0f-%.0f)',median(x),quantile(x,.25),quantile(x,.75))
  k<-k+1;desc[[k]]<-data.frame(universe=universe,histology=h,variable=v,level='Median (IQR)',included=fmt(a[[v]]),excluded=fmt(b[[v]]),included_n=nrow(a),excluded_n=nrow(b),included_percent=NA,excluded_percent=NA,SMD=smd)
 }
 for(v in c('sex','stage','era','site','surgery','chemotherapy','confirmation'))for(l in sort(unique(z[[v]]))){
  na<-sum(a[[v]]==l);nb<-sum(b[[v]]==l);pa<-na/nrow(a);pb<-nb/nrow(b)
  k<-k+1;desc[[k]]<-data.frame(universe=universe,histology=h,variable=v,level=l,included=sprintf('%d (%.1f%%)',na,100*pa),excluded=sprintf('%d (%.1f%%)',nb,100*pb),included_n=nrow(a),excluded_n=nrow(b),included_percent=100*pa,excluded_percent=100*pb,SMD=smd_binary(pa,pb))
 }
}
desc<-do.call(rbind,desc);savecsv(desc,'included_excluded_characteristics')
rm(raw)
d<-readRDS(file.path(root,'02_processed_data/pacc_pdac_analysis_cohort.rds'))
stopifnot(all(flow$final==c(sum(d$histology=='PDAC'),sum(d$histology=='pACC'))))
trend<-list();k<-0
for(group in c('pACC','PDAC pooled','PDAC 8500')){
 z<-d[d$stage!='Unknown' & (if(group=='pACC') d$histology=='pACC' else if(group=='PDAC pooled')d$histology=='PDAC' else d$hist_code=='8500'),]
 z<-droplevels(z);z$year_centered<-z$year-2004
 z$trend_era<-factor(ifelse(z$year<=2009,'2004-2009',ifelse(z$year<=2015,'2010-2015','2016-2023')))
 for(model in c('year','era')){
 f<-coxph(as.formula(paste('Surv(os_time,os_event)~',if(model=='year')'year_centered' else 'trend_era','+age+sex+race3+origin+stage+site+cluster(patient_id)')),data=z,ties='efron')
 sm<-summary(f);ii<-grep(if(model=='year')'^year_centered' else '^trend_era',rownames(sm$coefficients))
 for(i in ii){k<-k+1;trend[[k]]<-data.frame(group=group,model=model,term=rownames(sm$coefficients)[i],n=f$n,deaths=f$nevent,HR=sm$conf.int[i,1],lower=sm$conf.int[i,3],upper=sm$conf.int[i,4],p=sm$coefficients[i,'Pr(>|z|)'])}
 }
}
trend<-do.call(rbind,trend);savecsv(trend,'calendar_trend_models')
csfun<-function(z,lm,hz){
 n0<-nrow(z);z<-z[z$os_time>lm,];r<-z$os_time-lm;n<-nrow(z)
 e<-l<-u<-NA_real_; supported<-n>=2 && max(r)>=hz
 if(supported){s<-summary(survfit(Surv(r,z$os_event)~1,conf.type='log-log'),times=hz,extend=FALSE);e<-s$surv;l<-s$lower;u<-s$upper}
 data.frame(n_diagnosis=n0,n_landmark=n,n_end=sum(r>=hz),deaths=sum(z$os_event==1&r<=hz),censored=sum(z$os_event==0&r<hz),estimate=e,lower=l,upper=u,supported=supported)
}
# Era labels are fixed, but horizon eligibility can truncate their represented years.
eras<-data.frame(era=c('2004-2009','2010-2015','2016-2023'),start=c(2004,2010,2016),end=c(2009,2015,2023))
rows<-list();k<-0
for(group in c('pACC','PDAC pooled','PDAC 8500'))for(hz in c(3,5))for(lm in c(0,1,2,3,5))for(j in 1:3){
 cutoff<-min(eras$end[j],2023-lm-hz)
 z<-d[d$year>=eras$start[j]&d$year<=cutoff & (if(group=='pACC')d$histology=='pACC' else if(group=='PDAC pooled')d$histology=='PDAC' else d$hist_code=='8500'),]
 k<-k+1;rows[[k]]<-cbind(data.frame(group=group,era=eras$era[j],diagnosis_start=eras$start[j],diagnosis_end=cutoff,landmark=lm,horizon=hz),csfun(z,lm,hz))
}
era_cs<-do.call(rbind,rows);savecsv(era_cs,'era_conditional_OS')
savecsv(as.data.frame(table(hist_code=d$hist_code,year=d$year)),'code_year_counts')
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
cat('FLOW\n');print(flow,row.names=FALSE)
cat('CHARACTERISTICS AT FILTER STEP\n');print(subset(desc,universe=='confirmation_step' & (variable=='age' | (variable=='stage' & level=='Distant') | (variable=='surgery' & level=='Yes') | variable=='confirmation')),row.names=FALSE)
cat('TRENDS\n');print(trend,row.names=FALSE)
cat('ERA CONDITIONAL OS HORIZON3\n');print(subset(era_cs,horizon==3 & landmark %in% c(0,3,5)),row.names=FALSE)
