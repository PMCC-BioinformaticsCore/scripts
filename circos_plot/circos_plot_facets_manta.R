suppressMessages(library(OmicCircos))
suppressMessages(library(dplyr))
suppressMessages(library(VariantAnnotation))
suppressMessages(library(StructuralVariantAnnotation))
suppressMessages(library(rtracklayer))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
options(stringsAsFactors=F)

# Script works for Manta version 1.5 
# And facets-suite:2.0.8
# For human genome hg19 and hg38
# Example inputs: facets: *_hisens.rds; manta: somaticSV.vcf.gz
# Usage: circos_plot_facets_manta.R $tumour_id $normal_id $facet_rds_file $manta_sv_file $output_dir $manta_filter_qual[optional, use "PASS" by default]

gen_link_dat = function(svdat, keep_deletions=T, valid_chroms = c(1:22,"X","Y")){
  svdat = svdat[svdat$FILTER == "PASS",]
  ids = svdat$ID
  pairs = as.list(rep(NA,length(ids)))
  pairs2 = as.list(rep(NA,length(ids)))
  i=0
  j=1
  for(id in ids){
    i=i+1
    info=strsplit(svdat$INFO[i],";",fixed=T)[[1]]
    info_split = lapply(info, function(x) strsplit(x,"=")[[1]])
    info_dict = sapply(info_split, function(x) unlist(tail(x,n=1)))
    names(info_dict) = sapply(info_split, function(x) x[1])
    gene = paste("Gene",i,sep="")
    pos = svdat$POS[i]
    chrom = svdat$CHROM[i]
    if(substr(id,1,8)=="MantaDEL" & keep_deletions){
      len=as.numeric(info_dict["SVLEN"])*-1
      pairs2[[i]] = c(chrom, pos, gene, chrom, pos+len, gene, "del")
    }else if(substr(id,1,8)=="MantaBND"){
      mate = match(info_dict["MATEID"], ids)
      if(mate > i){
        info=strsplit(svdat$INFO[mate],";",fixed=T)[[1]]
        info_split = lapply(info, function(x) strsplit(x,"=")[[1]])
        info_dict = sapply(info_split, function(x) unlist(tail(x,n=1)))
        names(info_dict) = sapply(info_split, function(x) x[1])
        gene2 = paste("Gene",mate,sep="")
        chrom2 = svdat$CHROM[mate]
        pos2 = svdat$POS[mate]
        if(chrom == chrom2){
          # inversion
          alt = svdat$ALT[i]
          pos2 = gsub("\\[.*$|\\].*$", "", gsub(".*:", "", alt))
          pairs2[[i]] = c(chrom, pos, gene, chrom, pos2, gene, "inv")
        } else {
          # breakend
          pairs[[i]] = c(chrom, pos, gene, chrom2, pos2, gene2)
        }
        
      }
    } else if(substr(id,1,8)=="MantaDUP"){
      len=as.numeric(info_dict["SVLEN"])
      pairs2[[i]] = c(chrom, pos, gene, chrom, pos+len, gene, "dup")
    }
  }
  pairs=pairs[!is.na(pairs)]
  pairs2=pairs2[!is.na(pairs2)]

  # Check if pairs is empty
  if (length(pairs) == 0) {
    pairsdf <- data.frame()
  } else {
    pairsdf = data.frame(chr1=rep(NA,length(pairs)),po1=NA,gene1=NA,chr2=NA,po2=NA,gene2=NA)
    for(i in 1:length(pairs)){
      for(j in 1:6){pairsdf[i,j] = pairs[[i]][j]}
    }
    pairsdf = pairsdf[pairsdf$chr1%in%valid_chroms & pairsdf$chr2%in%valid_chroms,]
  }
  # Check if pairs is empty
  if (length(pairs2) == 0) {
    samechrompairsdf <- data.frame()
  } else {
    samechrompairsdf = data.frame(chr1=rep(NA,length(pairs2)),po1=NA,gene1=NA,chr2=NA,po2=NA,gene2=NA, type=NA)
    for(i in 1:length(pairs2)){
      for(j in 1:7){samechrompairsdf[i,j] = pairs2[[i]][j]}
    }
    samechrompairsdf = samechrompairsdf[samechrompairsdf$chr1%in%valid_chroms & samechrompairsdf$chr2%in%valid_chroms,]
  }

  return(list(pairsdf,samechrompairsdf))
}

# Agruments
args = commandArgs(trailingOnly = T)
samplename=args[1]
samplename_normal=args[2]
cnvf=args[3]
svf=args[4]
outdir=args[5]
genome=args[6]
if ( length(args) == 7  ) {
  manta_somaticsore=args[7]
  if ( is.integer(manta_somaticsore) ) { apply_filter=TRUE }
  else {
    print("Error: Please specify integer for manta filter score")
    quit()
   }
} else {
  apply_filter=FALSE
}

# MultiOmics anchor
data(UCSC.hg19.chr)
UCSC.hg19.chr$chrom = gsub("chr","",UCSC.hg19.chr$chrom)
seg.num=length(unique(UCSC.hg19.chr$chrom));
seg.name = unique(UCSC.hg19.chr$chrom)

if (genome=="hg19") {
  # hg19
  db=segAnglePo(seg.dat = UCSC.hg19.chr, seg = seg.name)
} else if (genome=="hg38") {
  # hg38
  seg_end <- vector()
  for(i in 1:seg.num) {
    seg_end <- c(seg_end, seqlengths(BSgenome.Hsapiens.UCSC.hg38)[i])
  }
  seg.dat <- data.frame(chr=seg.name, start=rep(1, seg.num), end=seg_end+1, name=seg.name, value=seg_end)
  db=segAnglePo(seg.dat = seg.dat, seg = seg.name)
} else {
  print("Error: unrecognised genome version.")
  print("Please specify the correct genome version: hg19, or hg38")
  quit()
}
colors = rainbow(seg.num, alpha=0.5)

print(paste("FACETS file:", cnvf))
mapdat_rd <- readRDS(cnvf)
mapdat <- data.frame(chr=c(mapdat_rd$segs[,"chrom"],1,1), 
                     start=c(mapdat_rd$segs[, "start"],0,0), 
                     end=c(mapdat_rd$segs[, "end"],0,0),
                     CN=pmin(c(mapdat_rd$segs[, "tcn.em"],0,4), 5))

print(paste("MANTA file:", svf))
svdat2 = read.delim(gzfile(svf), sep="\t",header=F, comment.char = "#",
                    col.names = c("CHROM", "POS", "ID", "REF", "ALT",
                                  "QUAL","FILTER","INFO","FORMAT",
                                  samplename_normal,samplename))
if (apply_filter==TRUE) {
  svdat = readVcf(svf)
  info_data <- info(svdat)
  keep <- rownames(info_data)[info_data$SOMATICSCORE >=manta_somaticsore]
  svdat <- svdat[keep, ]
  svdat2 = svdat2[svdat2$ID %in% keep,]
}

linkdat = gen_link_dat(svdat2)
cnvcols = sapply(mapdat$CN, function(x) {cn=ifelse(x>3, 4, x); c("blue","blue","black","red",rep("red",20))[cn+1]})

# Draw plot
print(paste("Creating output directory:", outdir, sep=" "))
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
outf=paste(outdir, "/", paste(samplename, samplename_normal, sep="--"),".pdf",sep="")
print(paste("Writing to file", outf, sep=" "))
pdf(outf,width=7,height=7)
par(mar=c(2,2,2,2))
plot(c(1,800), c(1,800), type="n", axes = F, xlab="", ylab="", main="", cex=2);
circos(R=400, cir=db, type="chr", col=colors, print.chr.lab=T, W=4, scale=T, cex=3)
circos(R=260, cir=db, W=120, mapping=mapdat, col.v = 4, type="arc", B=T, cutoff=2, lwd=4, col=cnvcols, scale=T,cex=10)
if ( ! ( length(linkdat[[1]]) == 0 ) ) {
  linkdat_breakends = linkdat[[1]]
  circos(R=260, cir=db, W=40, mapping=linkdat_breakends, type="link", lwd=2, col= rainbow(nrow(linkdat_breakends), alpha=0.5))
}
if ( ! ( length(linkdat[[2]]) == 0 ) ) {
  linkdat_samechrom = linkdat[[2]]
  circos(R=260, cir=db, W=20, mapping=linkdat_samechrom, type="link2", lwd=1,
         col=c("darkblue","green", 'orange')[match(linkdat_samechrom$type, c("del", "dup", "inv"))])
}
dev.off()
