suppressMessages(library("facets", quietly = TRUE))
suppressMessages(library("argparse", quietly = TRUE))

parser <- ArgumentParser(description="Run or Refit Facets")
parser$add_argument(
  "--outdir",
  required=TRUE,
  help="Output Directory")

parser$add_argument(
  "--output_prefix",
  required=FALSE,
  default="facets",
  help="Prefix for output files")

parser$add_argument(
  "--dbsnp",
  required=TRUE,
  help="dbsnp reference file")

parser$add_argument(
  "--tumor",
  required=TRUE,
  help="Tumor Bam")

parser$add_argument(
  "--normal",
  required=TRUE,
  help="Normal Bam")

parser$add_argument(
  "--unmatched",
  required=FALSE,
  help="Run Unmatched",
  action="store_true", default=FALSE)

parser$add_argument(
  "--cval",
  required=FALSE,
  help="Critical value for segmentation, recommend to use 150 for exome and 300 for genome. The lower number, the higher sensitivity.",
  default=150)

parser$add_argument(
  "--diplogR",
  required=FALSE,
  default=FALSE,
  help="Refit with Alternate DiplogR")

opt <- parser$parse_args()
setwd(opt$outdir) # run in outdir
set.seed(1234) # for reproducibility
snpout <- paste(opt$output_prefix, "snp_pileup_output.txt", sep="_")

# set up output file names
output.clonal.tab<-paste(opt$output_prefix,".clonal.tsv",sep="")
output.clonal.fig<-paste(opt$output_prefix,".clonal.pdf",sep="")
output.subclonal.tab<-paste(opt$output_prefix,".subclonal.tsv",sep="")
output.subclonal.fig<-paste(opt$output_prefix,".subclonal.pdf",sep="")

if (opt$diplogR == FALSE){
  system(paste(
    "snp-pileup -g -q30 -Q30 -r10,10",
    opt$dbsnp, snpout, opt$normal, opt$tumor))
}

rcmat <- readSnpMatrix(paste(snpout, ".gz", sep=''))
xx = preProcSample(rcmat, unmatched=opt$unmatched)

if(opt$diplogR != FALSE) {
  oo <- procSample(xx,cval=opt$cval, dipLogR = opt$d)
} else {
  oo <- procSample(xx,cval=opt$cval)
}

# clone
fit=emcncf(oo)
write.table(fit$cncf,file=output.clonal.tab,sep="\t",quote=F,row.names=F)
pdf(output.clonal.fig)
plotSample(x=oo,emfit=fit)
dev.off()

# subclone
fit2=emcncf2(oo)
write.table(fit2$cncf,file=output.subclonal.tab,sep="\t",quote=F,row.names=F)
pdf(output.subclonal.fig)
plotSample(x=oo,emfit=fit2)
dev.off()

file=paste(opt$output_prefix,"_estimates.txt",sep="")
sink(file=file,type="output")
paste("Clonal purity=",fit$purity,sep="")
paste("Clonal ploidy=",fit$ploidy,sep="")
paste("Subclonal purity=",fit2$purity,sep="")
paste("Subclonal ploidy=",fit2$ploidy,sep="")
paste("Original alBalLogR=")
oo$alBalLogR
sink()

rawFile=paste(opt$output_prefix,"_rawData.tsv",sep="")
jseg<-oo$jointseg
write.table(jseg,file=rawFile,quote=F,sep="\t",row.names=F)
