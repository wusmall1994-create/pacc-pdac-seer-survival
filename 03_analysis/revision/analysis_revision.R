# Targeted Pancreatology major-revision analyses. No original files are modified.
# Run: Rscript --vanilla analysis_revision.R [project_root] [output_directory]
options(stringsAsFactors=FALSE, width=180)
suppressPackageStartupMessages(library(survival))
args <- commandArgs(trailingOnly=TRUE)
project <- if(length(args)>=1) args[1] else '.'
out <- if(length(args)>=2) args[2] else file.path(project,'04_results/revision/conditional_OS')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
save_table <- function(x,name) write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='NA',fileEncoding='UTF-8')
d <- readRDS(file.path(project,'02_processed_data/pacc_pdac_analysis_cohort.rds'))
stopifnot(nrow(d)==116516, sum(d$histology=='pACC')==580, all(is.finite(d$os_time)), all(d$os_time>0), all(d$year %in% 2004:2023))
input_paths <- file.path(project,c('02_processed_data/pacc_pdac_analysis_cohort.rds','01_data_extraction/pacc_pdac_seer17_2004_2023_raw.txt','01_data_extraction/pacc_pdac_seer17_2004_2023_raw.dic','01_data_extraction/pacc_pdac_seer17_2004_2023_case_listing.sl','03_analysis/01_core_analysis.R','03_analysis/02_extended_analysis.R'))
present <- file.exists(input_paths)
save_table(data.frame(file=substring(input_paths[present],nchar(project)+2),md5=unname(tools::md5sum(input_paths[present]))),'input_provenance')

# Independently reconstruct eligibility from the untouched export and compare to RDS.
raw <- read.delim(input_paths[2],sep='\t',quote='"',na.strings='NA',check.names=FALSE,colClasses='character')
keep <- raw[['Histologic Type ICD-O-3']] %in% c('8140','8500','8550') & raw[['Diagnostic Confirmation']] %in% c('Positive histology','Pos hist AND immunophenotyping AND/OR pos genetic studies') & !raw[['Type of Reporting Source']] %in% c('Autopsy only','Death certificate only') & !is.na(raw[['Survival Days']]) & !is.na(raw[['Vital status recode (study cutoff used)']])
raw_d <- raw[keep,,drop=FALSE]
stopifnot(nrow(raw_d)==nrow(d))
for(nm in names(raw_d)) stopifnot(identical(raw_d[[nm]],d[[nm]]))
stopifnot(all(d$os_time==pmax(as.numeric(raw_d[['Survival Days']]),0.5)/365.25), all(d$os_event==as.integer(raw_d[['Vital status recode (study cutoff used)']]=='Dead')))
raw_n <- nrow(raw); rm(raw,raw_d)

# Definition is alive AND still observed at landmark; > and >= coincide here.
# No extrapolation past last observed follow-up. CI: Greenwood/log-log.
estimate_cs <- function(z,lm,hz) {
  boundary_n <- sum(z$os_time==lm)
  z <- z[z$os_time>lm,,drop=FALSE]
  r <- z$os_time-lm
  n <- nrow(z); n_end <- sum(r>=hz)
  ev <- sum(z$os_event==1 & r<=hz)
  cen <- sum(z$os_event==0 & r<hz)
  deaths_all <- sum(z$os_event)
  max_r <- if(n) max(r) else NA_real_
  est <- lo <- hi <- ratio_error <- NA_real_
  supported <- n>=2 && is.finite(max_r) && max_r>=hz
  if(supported) {
    fit <- survfit(Surv(r,z$os_event)~1,conf.type='log-log')
    s <- summary(fit,times=hz,extend=FALSE)
    stopifnot(length(s$surv)==1)
    est <- unname(s$surv); lo <- unname(s$lower); hi <- unname(s$upper)
  }
  data.frame(n_landmark=n,n_unique_patients=length(unique(z$patient_id)),n_at_horizon=n_end,
    deaths_within_horizon=ev,censored_before_horizon=cen,deaths_all_remaining_followup=deaths_all,
    max_remaining_followup_years=max_r,estimate=est,lower95=lo,upper95=hi,
    ci_width_pp=100*(hi-lo),supported=supported,
    sparse_landmark=n<20,sparse_horizon=n_end<20,exact_landmark_boundary=boundary_n)
}

# Validate original outputs using their exact original algorithm (including >=).
old <- read.csv(file.path(project,'04_results/tables/05_conditional_survival.csv'))
old <- old[old$endpoint=='OS',]
validation <- lapply(seq_len(nrow(old)),function(i){
  q <- old[i,]; z <- d[d$histology==q$histology & (q$stage=='All' | d$stage==q$stage),]
  z <- z[z$os_time>=q$landmark,]
  if(nrow(z)<2) return(NULL)
  r <- z$os_time-q$landmark
  f <- survfit(Surv(r,z$os_event)~1,conf.type='log-log')
  s <- summary(f,times=q$horizon,extend=TRUE)
  data.frame(histology=q$histology,stage=q$stage,landmark=q$landmark,horizon=q$horizon,
    old_estimate=q$estimate,reproduced_estimate=s$surv,absolute_error=abs(q$estimate-s$surv),
    lower_error=abs(q$lower-s$lower),upper_error=abs(q$upper-s$upper),
    n_matches=nrow(z)==q$n_risk,max_remaining=max(r),extrapolated=q$horizon>max(r))
})
validation <- do.call(rbind,validation)
stopifnot(max(validation$absolute_error,na.rm=TRUE)<1e-10,all(validation$n_matches),max(validation$lower_error,na.rm=TRUE)<1e-10,max(validation$upper_error,na.rm=TRUE)<1e-10)
save_table(validation,'original_OS_reproduction')

rows <- list(); k <- 0
for(comp in c('8140+8500','8500-only')) for(hz in c(3,5)) for(lm in c(0,1,2,3,5)) {
  # matched_window varies by landmark; common_window holds diagnosis cohort fixed.
  for(cohort in c('full','matched_window','common_window')) {
    cutoff <- switch(cohort,full=2023,matched_window=2023-lm-hz,common_window=2023-5-hz)
    dc <- d[d$year<=cutoff & (comp=='8140+8500' | d$hist_code %in% c('8500','8550')),]
    for(st in c('All','Localized','Regional','Distant','Unknown')) for(h in c('pACC','PDAC')) {
      z <- dc[dc$histology==h & (st=='All' | dc$stage==st),]
      e <- estimate_cs(z,lm,hz)
      # Independent probability identity check S(lm+hz)/S(lm).
      if(e$supported) {
        f0 <- survfit(Surv(z$os_time,z$os_event)~1,conf.type='log-log')
        ss <- summary(f0,times=c(lm,lm+hz),extend=FALSE)$surv
        ratio_error <- abs(e$estimate-ss[2]/ss[1])
        stopifnot(is.finite(ratio_error),ratio_error<1e-10)
      } else ratio_error <- NA_real_
      k<-k+1
      rows[[k]] <- cbind(data.frame(comparator=comp,cohort=cohort,diagnosis_start=2004,diagnosis_end=cutoff,
        histology=h,stage=st,landmark=lm,horizon=hz,n_diagnosis=nrow(z)),e,ratio_check_error=ratio_error)
    }
  }
}
cs <- do.call(rbind,rows)
stopifnot(all(cs$estimate[cs$supported]>=0 & cs$estimate[cs$supported]<=1), all(cs$exact_landmark_boundary==0))
save_table(cs,'conditional_OS_all_analyses')
save_table(cs[cs$stage=='All',],'conditional_OS_overall')
save_table(cs[cs$landmark==5 & cs$stage!='Unknown',],'conditional_OS_landmark5')

keys <- c('comparator','cohort','diagnosis_start','diagnosis_end','stage','landmark','horizon')
fields <- c('n_diagnosis','n_landmark','n_at_horizon','estimate','lower95','upper95','supported')
pa <- cs[cs$histology=='pACC',c(keys,fields)]; pd <- cs[cs$histology=='PDAC',c(keys,fields)]
names(pa)[match(fields,names(pa))] <- paste0('pACC_',fields)
names(pd)[match(fields,names(pd))] <- paste0('PDAC_',fields)
contrast <- merge(pa,pd,by=keys)
contrast$gap_pp <- 100*(contrast$pACC_estimate-contrast$PDAC_estimate)
save_table(contrast,'conditional_OS_descriptive_gaps')

# Membership/composition and support audit; aggregate only, no patient identifiers.
counts <- as.data.frame(table(hist_code=d$hist_code,diagnosis_year=d$year,stage=d$stage))
save_table(counts,'cohort_counts_by_code_year_stage')
landmark_year <- do.call(rbind,lapply(c(0,1,2,3,5),function(lm){
  z<-d[d$os_time>lm,]; a<-as.data.frame(table(hist_code=z$hist_code,diagnosis_year=z$year)); a$landmark<-lm; a
}))
save_table(landmark_year,'landmark_diagnosis_year_composition')
qc <- data.frame(check=c('raw_export_records','raw_to_processed_all_columns_identical','final_records','pACC_records','PDAC_8500_records','followup_beyond_2023_year_upper_bound','original_OS_max_absolute_error','original_OS_extrapolated_rows','original_OS_extrapolated_known_stage_or_overall_rows','new_CS_ratio_max_error','exact_landmark_boundary_records_across_runs'),
  value=c(raw_n,TRUE,nrow(d),sum(d$histology=='pACC'),sum(d$hist_code=='8500'),sum(d$os_time>(2024-d$year)+0.01),max(validation$absolute_error),sum(validation$extrapolated),sum(validation$extrapolated & validation$stage!='Unknown'),max(cs$ratio_check_error,na.rm=TRUE),sum(cs$exact_landmark_boundary)))
save_table(qc,'validation_summary')
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
cat('\nVALIDATION\n'); print(qc,row.names=FALSE)
cat('\nOVERALL LANDMARK 5\n'); print(cs[cs$stage=='All' & cs$landmark==5,c('comparator','cohort','diagnosis_end','histology','horizon','n_diagnosis','n_landmark','n_at_horizon','deaths_within_horizon','censored_before_horizon','estimate','lower95','upper95')],row.names=FALSE)
cat('\nCOMMON COHORT 0 TO 5 LANDMARKS\n'); print(contrast[contrast$stage=='All' & contrast$cohort!='matched_window' & contrast$landmark %in% c(0,5),c('comparator','cohort','diagnosis_end','landmark','horizon','pACC_estimate','PDAC_estimate','gap_pp')],row.names=FALSE)
message('Targeted revision analyses and validation complete.')
