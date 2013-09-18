#!/usr/bin/env Rscript

# Author: John Davey johnomics@gmail.com
# Begun 18/09/2013

#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################


# Initialise
library(grid)
library(RColorBrewer)

genetic.colours<-c(rgb(244,165,130,max=255),rgb(5,113,176,max=255))

args<-commandArgs(trailingOnly=T)

read.delim(args[2])->chrmap
read.delim(args[3])->scfmap

read.delim(args[4],header=F)->scflen
names(scflen)=c("Scaffold","Length")
scflen<-scflen[rev(order(scflen$Length)),]
rownames(scflen)<-1:nrow(scflen)

chrmap.col<-data.frame(chrmap,Colour=rep(genetic.colours,nrow(chrmap)/2))

chrnum<-max(chrmap$Chromosome)

chrcol<-rainbow(chrnum)

maxcm<-max(chrmap$cM)
maxbp<-max(scfmap$ChrStart+scfmap$Length)

pdf(args[1], width=2*chrnum, height=10)
grid.newpage()

plotchrom<-function(chr,cmmap,scfcmmap) {
    
    print(chr)
    
    chromvp<-viewport(x=(chr*2-2)/(chrnum*2),y=0,width=1/chrnum,height=1,just=c("left","bottom"))
    pushViewport(chromvp)
    
    
    headervp<-viewport(0,0.9,width=1,height=0.1,just=c("left","bottom"))
    pushViewport(headervp)
    grid.text(chr,x=0,y=1,just=c("left","bottom"),gp=gpar(col=chrcol[chr]))
    popViewport() #headervp
    
    
    geneticvp<-viewport(0,0,width=0.5,height=0.9,yscale=c(maxcm,0),just=c("left","bottom"))
    pushViewport(geneticvp)
    grid.lines(c(0.6,0.6),unit(c(max(cmmap$cM),0),"native"),gp=gpar(col=rgb(141,160,203,max=255),lwd=3,lineend="round"))
    grid.text(cmmap$cM,0.35,unit(cmmap$cM,"native"),just=c("right","centre"),gp=gpar(fontsize=8))
    grid.polyline(
        c(rep(0.4,nrow(cmmap)),rep(0.8,nrow(cmmap))),
        unit(c(cmmap$cM,cmmap$cM),"native"),
        id=rep(1:nrow(cmmap),2),
        gp=gpar(col=genetic.colours,lwd=2,lineend="round")
    )
    popViewport() #geneticvp
    
    
    physicalvp<-viewport(0.5,0,width=0.5,height=0.9,yscale=c(maxbp,0),just=c("left","bottom"))
    pushViewport(physicalvp)

    apply (cmmap,1,
    	function(x) {
    		popViewport()
    		pushViewport(geneticvp)
    		grid.move.to(0.8,unit(x[2],"native"))
    		popViewport()
    		pushViewport(physicalvp)
    		midpoint = as.numeric(x[3])+as.numeric(x[4])/2
    		grid.line.to(0.1,unit(midpoint,"native"),gp=gpar(col=x[5],lwd=1,lty="dashed"))
    		grid.lines(c(0.1,0.1),unit(c(as.numeric(x[3])+20000,as.numeric(x[3])+as.numeric(x[4])-20000),"native"),gp=gpar(col=x[5],lwd=1,lineend="round"))
    	}
    )

    # physical scale
    chrmaxbp<-max(scfcmmap$ChrStart+scfcmmap$Length)
    mb.onemil.tick<-seq(8000,chrmaxbp,1000000)
    chrmb<-ceiling(chrmaxbp/1000000)
    grid.text(sprintf("%2d", 0:chrmb), 0.29, unit(mb.onemil.tick, "native"), just="right", gp=gpar(fontsize=4))
    grid.text("Mb",0.29, unit(chrmaxbp,"native"),just="right",gp=gpar(fontsize=4))
    grid.lines(c(0.4,0.4),unit(c(0,chrmaxbp),"native"),gp=gpar(col="grey",lineend="round"))
    grid.polyline(
    	c(rep(0.35,chrmb),rep(0.4,chrmb)),
    	unit(c(mb.onemil.tick,mb.onemil.tick),"native"),
    	id=rep(1:chrmb,2),
    	gp=gpar(col="grey",lwd=1,lineend="round")
    )



    scffontface<-sapply(scfcmmap$ScfOriented, function(x){if (x>0) "bold" else "plain"})
    scfxpos<-sapply(scfcmmap$ScfOriented, function(x){if (x>0) 0.5 else 0.75})
    scfmisassembled<-scfcmmap$ScfChroms > 1 | scfcmmap$ScfGaps > 0
    scfnamecol<-sapply(scfmisassembled, function(x){if (x==TRUE) "red" else "black"})
    
    grid.polyline(
        c(scfxpos,scfxpos),
        unit(c(scfcmmap$ChrStart,scfcmmap$ChrStart+scfcmmap$Length),"native"),
        id=rep(1:nrow(scfcmmap),2),
        gp=gpar(col=chrcol[chr],alpha=c(1,0.3),lwd=1,lineend="butt")
    )
    grid.text(scfcmmap$Scaffold,scfxpos+0.05,unit(scfcmmap$ChrStart+scfcmmap$Length/2,"native"),just=c("left","centre"),gp=gpar(col=scfnamecol, fontsize=2,fontface=scffontface))
    popViewport() #physicalvp
    
    
    popViewport() #chromvp
}

# Main viewport: set margins, common parameters
pushViewport(viewport(x=0.5,y=0.49,width=0.95,height=0.95,gp=gpar(lineend="butt")))

scaffoldvp<-viewport(x=0,y=0,width=1,height=0.2,just=c("left","bottom"),xscale=c(0,nrow(scflen)*3),yscale=c(0,max(scflen$Length)))
pushViewport(scaffoldvp)
grid.polyline(
    unit(c(1:nrow(scflen)*3-3,1:nrow(scflen)*3-3),"native"),
    unit(c(rep(0,nrow(scflen)),scflen$Length),"native"),
    id=rep(1:nrow(scflen),2),
    gp=gpar(lwd=0.1)
)

scfmapnum<-sapply(scfmap$Scaffold,function(x) which(scflen$Scaffold==as.character(x)))

grid.polyline(
    unit(c(scfmapnum*3-3,scfmapnum*3-3),"native"),
    unit(c(scfmap$ScfStart,scfmap$ScfEnd),"native"),
    id=rep(1:nrow(scfmap),2),
    gp=gpar(col=chrcol[scfmap$Chromosome],lwd=0.1)
)


popViewport() #scaffoldvp

genomevp<-viewport(x=0,y=0.2,width=1,height=0.8,just=c("left","bottom"))
pushViewport(genomevp)
for (chr in 1:21) {
    plotchrom(chr, chrmap.col[chrmap$Chromosome==chr,], scfmap[scfmap$Chromosome==chr,])
}
popViewport() #genomevp

dev.off()
quit()


# Draw scaffolds
popViewport()
pushViewport(physicalvp)

mb.onemil.tick<-seq(8000,16000000,1000000)
mb.tenmil.tick<-seq(8000,16000000,10000000)
mb.fivemil.tick<-seq(8000,16000000,5000000)
grid.text(sprintf("%2d", 0:15), 0.29, unit(mb.onemil.tick, "native"), just="right", gp=gpar(fontsize=8))
grid.polyline(
	c(rep(0.35,16),rep(0.4,16)),
	unit(c(mb.onemil.tick,mb.onemil.tick),"native"),
	id=rep(1:16,2),
	gp=gpar(col="grey",lwd=1,lineend="round")
)
grid.polyline(
	c(rep(0.35,4),rep(0.4,4)),
	unit(c(mb.fivemil.tick,mb.fivemil.tick),"native"),
	id=rep(1:4,2),
	gp=gpar(col="darkgrey",lwd=2,lineend="round")
)
grid.polyline(
	c(rep(0.35,2),rep(0.4,2)),
	unit(c(mb.tenmil.tick,mb.tenmil.tick),"native"),
	id=rep(1:2,2),
	gp=gpar(col="dimgrey",lwd=3,lineend="round")
)


chr.width=0.49
chr.start=0.4
grid.rect(
	y=unit(c(chr18.scf.pos$Start),"native"),
	height=unit(c(chr18.scf.pos$End-chr18.scf.pos$Start),"native"),
	x=rep(chr.start,length(chr18.scf.pos$Scaffold)),
	width=rep(chr.width,length(chr18.scf.pos$Scaffold)),
	just=c("left","bottom"),
    gp=gpar(fill=rgb(141,160,203,max=255),col="black")
)


ordered.scfs<-chr18.scf.pos[chr18.scf.pos$Link=="black",]
grid.polyline(
    y=unit(c(ordered.scfs$Start,ordered.scfs$End), "native"),
    x=c(rep(chr.start,length(ordered.scfs$Scaffold)),rep(chr.start,length(ordered.scfs$Scaffold))),
    id=rep(1:length(ordered.scfs$Scaffold),2),
    gp=gpar(lwd=4,col=as.character(ordered.scfs$Link),lineend="round")
)
grid.polyline(
    y=unit(c(ordered.scfs$Start,ordered.scfs$End), "native"),
    x=c(rep(chr.start+chr.width,length(ordered.scfs$Scaffold)),rep(chr.start+chr.width,length(ordered.scfs$Scaffold))),
    id=rep(1:length(ordered.scfs$Scaffold),2),
    gp=gpar(lwd=4,col=as.character(ordered.scfs$Link),lineend="round")
)

dev.off()