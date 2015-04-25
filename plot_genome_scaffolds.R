#!/usr/bin/env Rscript

# Author: John Davey johnomics@gmail.com
# Begun 25/04/2015 based on plot_genome_map.R

#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################


# Initialise
library(grid)
library(plyr)
library(RColorBrewer)

args<-commandArgs(trailingOnly=T)

plotchrom<-function(chr,scfmap) {

    
    grid.newpage()

    chrscf <- scfmap[scfmap$Chromosome==chr,]
    maxlengths <- tapply(chrscf$Length, chrscf$Pool, sum)
    parts <- tapply(chrscf$ID, chrscf$Pool, max)
    chrmaxparts <- max(parts)
    chrsize <- sum(maxlengths)
    cumlen <- cumsum(maxlengths)
    names(cumlen) <- NULL
    starts <- c(0, cumlen[1:(length(cumlen)-1)])
    print(chr)

    pushViewport(viewport(x=0.5,y=0.49,width=0.95,height=0.95,gp=gpar(lineend="butt")))

    headervp<-viewport(0,0.9,width=1,height=0.1,just=c("left","bottom"))
    pushViewport(headervp)
    grid.text(chr,x=0,y=1,just=c("left","bottom"))
    popViewport() #headervp
    popViewport()

    physicalvp<-viewport(0,0,width=1,height=0.9,yscale=c(chrsize+1000000,0),just=c("left","bottom"))
    pushViewport(physicalvp)

    # physical scale
    mb.onemil.tick<-seq(8000,chrsize+8000,1000000)
    chrmb<-ceiling(chrsize/1000000)
    grid.text(sprintf("%2d", 0:chrmb), 0.1, unit(mb.onemil.tick, "native"), just="right")
    grid.text("Mb",0.1, unit(chrsize,"native"),just="right")
    grid.lines(c(0.2,0.2),unit(c(0,chrsize),"native"),gp=gpar(lwd=4, col="grey",lineend="round"))
    grid.polyline(
        c(rep(0.15,chrmb),rep(0.2,chrmb)),
        unit(c(mb.onemil.tick,mb.onemil.tick),"native"),
        id=rep(1:chrmb,2),
        gp=gpar(col="grey",lwd=3,lineend="round")
    )

    scfcol<-sapply(scfmap$Type, function(x){if (x=="ok") "green4" else if (x=="orient") "orange" else "red"})
    scfxpos<-sapply(scfmap$Type, function(x){if (x=="ok") 0.25 else if (x=="orient") 0.3 else 0.35})

    print(length(scfmap$Type))
    print(length(scfxpos))
    print(scfmap$Pool)
    print(starts)
    print(length(starts[scfmap$Pool]))

    offset = (scfmap$ID-1)*(0.3/chrmaxparts)
    grid.polyline(
        unit(c(scfxpos+offset,scfxpos+offset),"native"),
        unit(c(starts[scfmap$Pool],starts[scfmap$Pool]+scfmap$Length),"native"),
        id=rep(1:nrow(scfmap),2),
        gp=gpar(col=scfcol,lwd=2,lineend="butt")
    )
    popViewport() #physicalvp
}

read.delim(args[2])->scfmap

pdf(args[1], width=8.27, height=11.69)

for (chr in 1:21) {
    plotchrom(chr, scfmap[scfmap$Chromosome==chr,])
}

dev.off()
quit()

