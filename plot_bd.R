#!/usr/bin/env Rscript

library(grid)


genetic.colours<-c(rgb(244,165,130,max=255),rgb(5,113,176,max=255)) # from plot_genome_scaffolds.R

geneticcolour<-function(x) ifelse(x == "blue", genetic.colours[2], ifelse(x=="red", genetic.colours[1], ifelse(x=="magenta", "red", "black")))


bd.length<-602276

rad.bd.snps.full<-read.delim("rad_bd_snps_full.tsv",stringsAsFactors=FALSE)

rad.bd.snps.orig<-read.delim("rad_bd_snps_orig.tsv",stringsAsFactors=FALSE)

wgs.bd.snps<-read.delim("wgs_bd_snps.tsv",stringsAsFactors=FALSE)

type.pos<-data.frame(Type=c("Maternal","Paternal","Intercross"),Pos=c(0.25,0.5,0.75))

rad.bd.ranges.orig<-data.frame(Start=c(29907,144657,332491,397442,474324),End=c(111817,214649,332491,397500,474390),Colour=c("blue","red","blue","red","blue"))

rad.bd.ranges.full<-data.frame(Start=c(23544,144544,237134,397267,474094,587513),End=c(111817,228639,376562,397685,474560,602037),Colour=c("blue","red","blue","red","blue","magenta"))

features<-data.frame(Name=c("Kinesin","Dennis","Rays","Optix"),Start=c(279741,323151,360978,438423),End=c(291617,330538,371228,439107))

wgs.bd.ranges<-data.frame(Start=c(236,131152,229607,377135,377992,394890,431078,586182),End=c(130929,228639,376983,377527,394160,430469,580977,602037),Colour=c("blue","red","blue","red","blue","red","blue","magenta"))

plot.snps<-function(snplist,ranges,min=0.1) {
    snpset.view<-viewport(0,min,just=c("left","bottom"),width=1,height=0.3,xscale=c(1,bd.length))
    pushViewport(snpset.view)
    snp.pos<-sapply(snplist$Type, function(x){type.pos[type.pos$Type==x,]$Pos})

    grid.polyline(
        unit(rep(snplist$Position,2),"native"),
        unit(c(snp.pos-0.05,snp.pos+0.05),"native"),
        id=rep(1:length(snplist$Position),2),
        gp=gpar(col=as.character(geneticcolour(snplist$Colour)))
    )
    
    grid.polyline(
        unit(c(ranges$Start,ranges$End),"native"),
        unit(rep(0.9,length(ranges$Start)*2),"native"),
        id=rep(1:length(ranges$Start),2),
        gp=gpar(col=as.character(geneticcolour(ranges$Colour)))
    )
    popViewport() #snpset.view
}

pdf("bd_plot_boston.pdf",width=10,height=7);
grid.newpage()
pushViewport(viewport(x=0.5,y=0.49,width=0.95,height=0.95))
scaffold.line<-viewport(0,0,just=c("left","bottom"),width=1,height=0.1,xscale=c(1,bd.length))
pushViewport(scaffold.line)

grid.lines(unit(c(1,bd.length),"native"),c(1,1),gp=gpar(lwd=1,lineend="round"))

ticks<-seq(0,bd.length,by=50000)
grid.text(sprintf("%7d",ticks), unit(ticks, "native"), 0.7, just="right",rot=90)
grid.polyline(
	unit(c(ticks,ticks),"native"),
	c(rep(0.8,length(ticks)),rep(1,length(ticks))),
	id=rep(1:length(ticks),2)
)
popViewport() #scaffold.line

snp.view<-viewport(0,0.1,just=c("left","bottom"),width=1,height=0.9,xscale=c(1,bd.length))
pushViewport(snp.view)
plot.snps(rad.bd.snps.orig,rad.bd.ranges.orig,0)
plot.snps(wgs.bd.snps,wgs.bd.ranges,0.35)
plot.snps(rad.bd.snps.full,rad.bd.ranges.full,0.7)

grid.polyline(
    unit(c(features$Start,features$End),"native"),
    unit(rep(0.04,length(features$Start)*2),"native"),
    id=rep(1:length(features$Start),2)
)
grid.text(
    features$Name,
    unit(features$Start+(features$End-features$Start)/2,"native"),
    unit(rep(0.02,length(features$Name)),"native")
)

popViewport() #snp.view
popViewport()