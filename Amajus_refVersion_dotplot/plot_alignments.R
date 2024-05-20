#Read in PAF:
#v2_v3.5.paf
#paf<-read.table("../v2_v3.5.paf", sep="\t", header=F)
#paf<-read.table("../Amajus_v3.5_SLocus.fasta_GWHBJVT00000000.genome.Chr.fa.paf", sep="\t", header=F)
scaf_v3.5<-read.table("../GWHBJVT00000000.genome.Chr.gaps", sep="\t", header=F)
scaf_v3.5<-split(scaf_v3.5, scaf_v3.5[,1])
scaf_v3<-read.table("../Amajus_v3.5_SLocus.gaps", sep="\t", header=F)
scaf_v3<-split(scaf_v3, scaf_v3[,1])


pdf("dotplot.minimap.pdf")
for (i in 1:8) {
	rows<-which(paf[,1] == paste0("Chr",i) & paf[,6] == paste0("Chr",i) & paf[,10]/paf[,11] > 0.5 & paf[,11] >= 10000)
	plot(range(paf[rows,3:4])/1e6, range(paf[rows,8:9])/1e6, type="n", xlab="v3 positions (Mbp)", ylab="v3.5 positions (Mbp)", main=paste0("Chr",i, " by minimap2"));#"Zebrafinch ordering");
	#abline(v=scaf_v2[[paste0("Chr",i)]][,2]/1e6,col=rgb(192,12,40,50,maxColorValue=255))
	abline(v=scaf_v3[[paste0("Chr",i)]][,2]/1e6,col=rgb(192,12,40,50,maxColorValue=255))
	abline(h=scaf_v3.5[[paste0("Chr",i)]][,2]/1e6,col=rgb(192,12,40,50,maxColorValue=255))
	abline(a=0,b=1,lty=2);
	#=scaf_v3.5[[paste0("Chr",i)]][,2]/1e6,col=rgb(192,12,40,50,maxColorValue=255))
	apply(paf[rows,],MARGIN=1, FUN=function(i) {if(i[5] == "+") {lines(as.numeric(i[3:4])/1e6,as.numeric(i[8:9])/1e6)} else {lines(as.numeric(i[4:3])/1e6, as.numeric(i[8:9])/1e6)}}); 
}
dev.off()
