# Run after the original pipeline; does not modify its results.
# Rscript --vanilla 03_analysis/run_revision.R [project_root] [revision_output]
args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(if(length(args)) args[1] else '.',mustWork=TRUE)
out <- if(length(args)>=2) args[2] else file.path(project,'04_results/revision')
file_arg <- grep('^--file=',commandArgs(),value=TRUE)
script_dir <- dirname(normalizePath(sub('^--file=','',file_arg[1]),mustWork=TRUE))
rscript <- file.path(R.home('bin'),if(.Platform$OS.type=='windows') 'Rscript.exe' else 'Rscript')
run <- function(script,arguments) {
 status <- system2(rscript,c('--vanilla',shQuote(file.path(script_dir,'revision',script)),shQuote(arguments)))
 if(status!=0) stop(script,' failed with exit status ',status)
}
run('analysis_revision.R',c(project,file.path(out,'conditional_OS')))
run('analysis_remaining.R',c(project,file.path(out,'selection_and_calendar')))
run('riskset_composition.R',c(project,file.path(out,'conditional_OS','riskset_composition.csv')))
run('flow_diagram.R',file.path(out,'figures'))
message('Revision analyses completed. Run the two independent Python validators as documented.')
