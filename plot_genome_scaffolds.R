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

genetic.colours<-c(rgb(244,165,130,max=255),rgb(5,113,176,max=255))


draw.header<-function(chr) {
    headervp<-viewport(0,0.9,width=1,height=0.1,just=c("left","bottom"))
    pushViewport(headervp)
    grid.text(chr,x=0,y=1,just=c("left","bottom"))
    popViewport() #headervp    
}

draw.linkage.map<-function(chrmap, scale) {
    maxcm <- max(chrmap$cM)
    geneticvp<-viewport(0.05, 0, width=0.89*scale, height=0.2, xscale=c(0, maxcm), just=c("left","bottom"))
    pushViewport(geneticvp)
    
    grid.lines(unit(c(0,maxcm),"native"),c(0.6,0.6),gp=gpar(col=rgb(141,160,203,max=255),lwd=3,lineend="round"))
    grid.text(chrmap$cM,unit(chrmap$cM,"native"),0.35,just=c("right","centre"),rot=90,gp=gpar(fontsize=8))
    grid.polyline(
        unit(c(chrmap$cM,chrmap$cM),"native"),
        c(rep(0.4,nrow(chrmap)),rep(0.8,nrow(chrmap))),
        id=rep(1:nrow(chrmap),2),
        gp=gpar(col=genetic.colours,lwd=2,lineend="round")
    )
    
    popViewport()
    
    return(geneticvp)
}

draw.physical.map<-function(scfmap, chr, baseline, scale=1, plotnames=FALSE) {
    chrscf <- scfmap[scfmap$Chromosome==chr,]
    poollengths <- tapply(chrscf$Length, chrscf$Pool, sum)
    parts <- tapply(chrscf$ID, chrscf$Pool, max)
    chrmaxparts <- max(parts)
    chrsize <- sum(poollengths)
    cumlen <- cumsum(poollengths)
    names(cumlen) <- NULL
    starts <- c(0, cumlen[1:(length(cumlen)-1)])

    physicalvp<-viewport(0.05,baseline,width=0.95*scale,height=0.35,xscale=c(0, chrsize+1000000),just=c("left","bottom"))
    pushViewport(physicalvp)

    # physical scale
    mb.onemil.tick<-seq(8000,chrsize+8000,1000000)
    chrmb<-ceiling(chrsize/1000000)
    grid.text(sprintf("%2d", 0:chrmb), unit(mb.onemil.tick, "native"), 0.15, just="right")
    grid.text("Mb", unit(chrsize+500000,"native"), 0.15, just="right")
    grid.lines(unit(c(0,chrsize), "native"), c(0.25,0.25), gp=gpar(lwd=4, col="grey",lineend="round"))
    grid.polyline(
        unit(c(mb.onemil.tick,mb.onemil.tick),"native"),
        c(rep(0.2,chrmb),rep(0.25,chrmb)),
        id=rep(1:chrmb,2),
        gp=gpar(col="grey",lwd=3,lineend="round")
    )

    scfmapbylength <- scfmap[with(scfmap, order(Length)),]
    scfcol<-sapply(scfmapbylength$Type, function(x){if (x=="ok") "green4" else if (x=="orient") "orange" else "red"})
    scfxpos<-sapply(scfmapbylength$Type, function(x){if (x=="ok") 0.3 else if (x=="orient") 0.35 else 0.4})
    offset = (scfmapbylength$ID-1)*((baseline+0.2)/chrmaxparts)

    grid.polyline(
        unit(c(starts[scfmapbylength$Pool],starts[scfmapbylength$Pool]+scfmapbylength$Length),"native"),
        unit(c(scfxpos+offset,scfxpos+offset),"native"),
        id=rep(1:nrow(scfmap),2),
        gp=gpar(col=scfcol,lwd=2,lineend="butt")
    )
    
    grid.text(
        scfmapbylength$Scaffold,
        unit(c(starts[scfmapbylength$Pool]),"native"),
        unit(c(scfxpos+offset+0.005),"native"),
        just=c("left","bottom"),
        gp=gpar(fontsize=2)
    )

    popViewport() #physicalvp
    
    return(physicalvp)
}

connect.maps<-function(chrmap, geneticvp, physicalvp) {
    colour.rows<-nrow(chrmap)/2
    chr.colours<-rep(genetic.colours, colour.rows)
    if (colour.rows != floor(colour.rows)) chr.colours<-c(chr.colours, genetic.colours[1])
    chrmap.col<-data.frame(chrmap,Colour=chr.colours)

    apply (chrmap.col,1,
        function(x) {
            pushViewport(geneticvp)
            grid.move.to(unit(x[2],"native"),0.8)
            popViewport()
            pushViewport(physicalvp)
            midpoint = (as.numeric(x[3])+as.numeric(x[4]))/2
            grid.line.to(unit(midpoint,"native"),0.043,gp=gpar(col=x[6],lwd=1,lty="dotted"))
            grid.lines(unit(c(as.numeric(x[3]),as.numeric(x[4])),"native"),c(0.05,0.05),gp=gpar(col=x[6],lwd=2,lineend="butt"))
            popViewport()
        }
    )
}

draw.breakpoints<-function(chr, physicalvp, breakpoints, chrsize) {
    if (is.null(breakpoints)) return()
    if (!(any(breakpoints$Chromosome==chr))) return()
    
    chrbreak<-breakpoints[breakpoints$Chromosome==chr,]
    pushViewport(physicalvp)
    grid.lines(unit(c(1,chrbreak$Start),"native"),c(0.07,0.07),gp=gpar(col='darkmagenta',lwd=2))
    grid.lines(unit(c(chrbreak$End,chrsize),"native"),c(0.07,0.07),gp=gpar(col='sienna',lwd=2))
    popViewport()
}

plotchrom<-function(chr,chrmap, revisedmap, draftmap, breakpoints) {

    grid.newpage()

    pushViewport(viewport(x=0.5,y=0.49,width=0.95,height=0.95,gp=gpar(lineend="butt")))
    
    draw.header(chr)

    draft_size <- sum(draftmap[draftmap$Chromosome==chr,]$Length)+1000000
    revised_size <- sum(revisedmap[revisedmap$Chromosome==chr,]$Length)+1000000
    if (draft_size > revised_size) {
        draft_scale <- 1
        revised_scale <- revised_size / draft_size
    }
    else {
        draft_scale <- draft_size / revised_size
        revised_scale <- 1
    }
    
    geneticvp<-draw.linkage.map(chrmap, revised_scale)

    draft<-draw.physical.map(draftmap, chr, 0.55, draft_scale)

    revisedvp<-draw.physical.map(revisedmap, chr, 0.2, revised_scale, plotnames=TRUE)
        
    connect.maps(chrmap, geneticvp, revisedvp)
    
    draw.breakpoints(chr, revisedvp, breakpoints, revised_size-1000000)
    
    popViewport()
}

read.delim(args[2])->scfmap
read.delim(args[3])->chrmap
read.delim(args[4])->draftmap

breakpoints<-NULL
if (length(args) > 4) {
    read.delim(args[5])->breakpoints
}
chrmap<-chrmap[chrmap$cM != -1,]

pdf(args[1], width=11.69, height=8.27)

for (chr in 1:max(chrmap$Chromosome)) {
    plotchrom(chr, chrmap[chrmap$Chromosome==chr,], scfmap[scfmap$Chromosome==chr,], draftmap[draftmap$Chromosome==chr,], breakpoints)
}

dev.off()
quit()

