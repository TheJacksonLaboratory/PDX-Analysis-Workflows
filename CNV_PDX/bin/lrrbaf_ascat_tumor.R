options(scipen = 999)

args=(commandArgs(TRUE))
snp_pos <- args[1]
gcfile <- args[2]

lrrbaf = read.table("lrr_baf1.txt", header = T, sep = "\t", row.names=1)

SNPpos = read.table(snp_pos,header=T,sep="\t",row.names=1)

firstline = read.table("lrr_baf1.txt", nrows=1, sep = "\t")
sample = sub(".CEL.Log.R.Ratio","",firstline[1,4])
#sample = sub(".CEL.Log.R.Ratio","",colnames(lrrbaf)[3])

Tumor_LogR = lrrbaf[rownames(SNPpos),3,drop=F]
colnames(Tumor_LogR) = sample

Tumor_BAF = lrrbaf[rownames(SNPpos),4,drop=F]
colnames(Tumor_BAF) = sample

#Normal_LogR = lrrbaf[rownames(SNPpos),5,drop=F]
#colnames(Normal_LogR) = sample

#Normal_BAF = lrrbaf[rownames(SNPpos),6,drop=F]
#colnames(Normal_BAF) = sample

#replace 2's by NA
Tumor_BAF[Tumor_BAF==2]=NA
#Normal_BAF[Normal_BAF==2]=NA

# Tumor_LogR: correct difference between copy number only probes and other probes
CNprobes = substring(rownames(SNPpos),1,2)=="CN"

Tumor_LogR[CNprobes,1] = Tumor_LogR[CNprobes,1]-mean(Tumor_LogR[CNprobes,1],na.rm=T)
Tumor_LogR[!CNprobes,1] = Tumor_LogR[!CNprobes,1]-mean(Tumor_LogR[!CNprobes,1],na.rm=T)

#Normal_LogR[CNprobes,1] = Normal_LogR[CNprobes,1]-mean(Normal_LogR[CNprobes,1],na.rm=T)
#Normal_LogR[!CNprobes,1] = Normal_LogR[!CNprobes,1]-mean(Normal_LogR[!CNprobes,1],na.rm=T)

# limit the number of digits:
Tumor_LogR = round(Tumor_LogR,4)
#Normal_LogR = round(Normal_LogR,4)

write.table(cbind(SNPpos,Tumor_BAF),paste(sample, ".tumor.BAF.txt", sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
#write.table(cbind(SNPpos,Normal_BAF),paste(sample, ".normal.BAF.txt", sep=""),sep="\t",row.names=T,col.names=NA,quote=F)

write.table(cbind(SNPpos,Tumor_LogR),paste(sample, ".tumor.LogR.txt", sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
#write.table(cbind(SNPpos,Normal_LogR),paste(sample, ".normal.LogR.txt", sep=""),sep="\t",row.names=T,col.names=NA,quote=F)

#run ASCAT functions

library(ASCAT)
file.tumor.LogR <- dir(pattern="tumor.LogR")
file.tumor.BAF <- dir(pattern="tumor.BAF")
#file.normal.LogR <- dir(pattern="normal.LogR")
#file.normal.BAF <- dir(pattern="normal.BAF")

gender <- read.table("gender.txt", sep="\t")
sex <- as.vector(gender[1,1])
sex[sex == "female"] <- "XX"
sex[sex == "male"] <- "XY"
sex[sex == "unknown"] <- "XX"

#samplename <- sub(".tumor.LogR.txt", "", file.tumor.LogR)

if (sex == "XX") {
    
    ascat.bc <- ascat.loadData(file.tumor.LogR, file.tumor.BAF, chrs=c(1:22, "X", "Y"), gender=sex)
    
} else if (sex == "XY") {
    
    ascat.bc <- ascat.loadData(file.tumor.LogR, file.tumor.BAF, chrs=c(1:22, "X","Y"), gender=sex)

}
#ascat.bc <- ascat.loadData(file.tumor.LogR, file.tumor.BAF, file.normal.LogR, file.normal.BAF, chrs=c(1:22, "X"), gender=sex)

#GC correction for SNP6 data
ascat.bc <- ascat.GCcorrect(ascat.bc, gcfile)

ascat.plotRawData(ascat.bc)

gg<-ascat.predictGermlineGenotypes(ascat.bc, platform = "AffySNP6")

ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)

ascat.plotSegmentedData(ascat.bc)

ascat.output = ascat.runAscat(ascat.bc)

#save ASCAT results

save.image(paste(sample,".RData",sep=""))

if ( length(ascat.output$failedarrays) == 0 ) {
    
    num_probes <- vector(mode="numeric", length=nrow(ascat.output$segments_raw))
    for (i in 1:nrow(ascat.output$segments_raw)) {
        
        #print(i)
        L1 = which(SNPpos$Chromosome == ascat.output$segments_raw$chr[i] & SNPpos$Physical.Position == ascat.output$segments_raw$startpos[i])
        L2 = which(SNPpos$Chromosome ==  ascat.output$segments_raw$chr[i] & SNPpos$Physical.Position == ascat.output$segments_raw$endpos[i])
        num_probes[i] = L2[length(L2)] - L1[1] + 1
        
    }
    seg_raw = cbind(ascat.output$segments_raw,num_probes)
    
    num_probes <- vector(mode="numeric", length=nrow(ascat.output$segments))
    for (i in 1:nrow(ascat.output$segments)) {
        
        #print(i)
        L1 = which(SNPpos$Chromosome == ascat.output$segments$chr[i] & SNPpos$Physical.Position == ascat.output$segments$startpos[i])
        L2 = which(SNPpos$Chromosome ==  ascat.output$segments$chr[i] & SNPpos$Physical.Position == ascat.output$segments$endpos[i])
        num_probes[i] = L2[length(L2)] - L1[1] + 1
        
    }
    seg = cbind(ascat.output$segments,num_probes)

    write.table(seg_raw, file=paste(sample,".segments_raw.txt",sep=""), sep="\t", quote=F, row.names=F)
    write.table(seg, file=paste(sample,".segments.txt",sep=""), sep="\t", quote=F, row.names=F)
    write.table(as.data.frame(ascat.output$aberrantcellfraction), file=paste(sample,".aberrantcellfraction.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
    write.table(as.data.frame(ascat.output$ploidy), file=paste(sample,".ploidy.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)

} else {
    
    write.table(as.data.frame(ascat.output$failedarrays), file=paste(sample,".failedarrays.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)

}

if ( !is.null(ascat.output$nonaberrantarrays) ) {

    write.table(as.data.frame(ascat.output$nonaberrantarrays), file=paste(sample,".nonaberrantarrays.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)

}
