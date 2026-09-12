# Figure S1 contract: transparent sequential selection; two histology branches;
# every exclusion is shown at its actual filtering step, without inferred causes.
# R-only schematic; 183 x 165 mm; vector PDF/SVG and 600-dpi TIFF.
args<-commandArgs(trailingOnly=TRUE)
out<-if(length(args))args[1] else 'figures'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
draw<-function(){
 par(mar=c(0,0,0,0),family='sans');plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
 boxtext<-function(x,y,w,h,label,size=.85){rect(x-w/2,y-h/2,x+w/2,y+h/2,col='#F3F5F7',border='#53616D');text(x,y,label,cex=size)}
 arr<-function(x,y,x2,y2)arrows(x,y,x2,y2,length=.065,col='#53616D',lwd=1)
 boxtext(.5,.91,.84,.10,'Pancreatic tumor records exported\nn = 160,472')
 arr(.5,.86,.5,.76)
 text(.53,.81,'Other exported histology codes excluded: 212',adj=0,cex=.70)
 boxtext(.5,.71,.84,.10,'Target histology codes 8140, 8500 and 8550\nn = 160,260')
 arr(.35,.66,.25,.58);arr(.65,.66,.75,.58)
 boxtext(.25,.53,.40,.09,'PDAC codes 8140 + 8500\nn = 159,595')
 boxtext(.75,.53,.40,.09,'pACC code 8550\nn = 665')
 segments(.25,.485,.25,.455,col='#53616D');segments(.75,.485,.75,.455,col='#53616D')
 arr(.25,.375,.25,.345);arr(.75,.375,.75,.345)
 text(.25,.415,'Histologic confirmation filter\nExcluded: 43,482 (27.2%)',cex=.78)
 text(.75,.415,'Histologic confirmation filter\nExcluded: 84 (12.6%)',cex=.78)
 boxtext(.25,.30,.40,.09,'Histologic confirmation recorded\nn = 116,113',.80)
 boxtext(.75,.30,.40,.09,'Histologic confirmation recorded\nn = 581',.80)
 segments(.25,.255,.25,.235,col='#53616D');segments(.75,.255,.75,.235,col='#53616D')
 arr(.25,.155,.25,.125);arr(.75,.155,.75,.125)
 text(.25,.195,'Autopsy / death certificate only: 177\nMissing survival / vital status: 0',cex=.72)
 text(.75,.195,'Autopsy / death certificate only: 1\nMissing survival / vital status: 0',cex=.72)
 boxtext(.25,.08,.40,.09,'Final PDAC cohort\nn = 115,936')
 boxtext(.75,.08,.40,.09,'Final pACC cohort\nn = 580')
}
pdf(file.path(out,'Figure_S1_flow.pdf'),width=7.2,height=6.5,useDingbats=FALSE);draw();dev.off()
svg(file.path(out,'Figure_S1_flow.svg'),width=7.2,height=6.5);draw();dev.off()
tiff(file.path(out,'Figure_S1_flow.tiff'),width=7.2,height=6.5,units='in',res=600,compression='lzw');draw();dev.off()
png(file.path(out,'Figure_S1_flow.png'),width=7.2,height=6.5,units='in',res=180);draw();dev.off()
