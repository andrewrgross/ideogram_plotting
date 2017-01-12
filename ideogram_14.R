### Ideogram generator -- Andrew R Gross -- 2015-10-02
### Second generation ideogram generator

### Header ####################################################################################
library(ggplot2)
library(GenomeGraphs)
require(IdeoViz)
require(RColorBrewer) ### nice colours

### Parse ideo d.f. for gieStain data #########################################################
###############################################################################################

ideo <- getIdeo("hg18")                     # downloads chromosome band data from the UCSC genome browser. 
ideoNames=ideo$chrom
ideoRanges = IRanges(ideo$chromStart,end=ideo$chromEnd,names=ideo$name)
#ideoStrand = rep("+",862)

### Functions #########################################################
### Subset a dataframe by specified rows #########################################################

subset<-function(dataframe,start,end) {
  newDF <- data.frame(as.matrix(dataframe[start:end,]))
  return(newDF)
}

### Subset a dataframe by specifying the chromosome position #########################################################
subsetByChrom<-function(dataframe,chromBegin,chromEnd) {
  #chrom <- chromosomes[chromNo]
  chromBegin <- paste0("chr",chromBegin)
  chromEnd <- paste0("chr",chromEnd)
  start <- which(dataframe[1] == chromBegin)[1]
  end  <-  which(dataframe[1] == chromEnd)[length(which(dataframe[1] == chromEnd))]
  newDF <- data.frame(as.matrix(dataframe[start:end,]))
  return(newDF)
}

### Format a genomic range with bin sizes corresponding to g-bands #########################################################
grFormat<-function(chromosomes) {
  chromLength <- c()
  Range <- c()
  Strand <- c()
  for (chrom in chromosomes) {
    print(chrom)
    chromStart <- which(ideo$chrom == chrom)[1]
    print(chromStart)
    chromEnd  <-  which(ideo$chrom == chrom)[length(which(ideo$chrom == chrom))]
    print(chromEnd)
    chromLength <- c(chromLength,chromEnd-chromStart+1)
    newRange = window(ideoRanges,chromStart,chromEnd) # Defines range for GRO by copying subsection of full ideo d.f. range
    Range <- append(Range,newRange)
  }
  Name <- Rle(chromosomes,chromLength)
  Strand = Rle(c("*"),length(Range))         # creates an RLE object of the appropriate length
  gr <- GRanges(seqnames = Name,             # Creates a GRO based on the given names, ranges, and strands
                ranges = Range,   
                strand = Strand)
  return(gr)
}

### formats abnormalities for plotting within a single chromosome #########################################################
chromFormat <- function(chromStart,chromEnd,chromLength,grSub,abnormSub) {
  for (i in 1:nrow(abnormSub)) {                      ### runs through the each abnormality listed within a chromosome
    band1 <- abnormSub[[2]][i]               # reads dataframe and defines current abnormality to plot
    band2 <- abnormSub[[3]][i]
    abnormType <- abnormSub[[4]][i]
    barPos1 <- grep(band1, names(grSub))
    barPos2 <- grep(band2, names(grSub))
    barX <- i*scaleFactor                                         # defines width of bars.  Can add scaling factor later
    #print(paste0(c("## band:","abnormType:"),c(as.character(band1),as.character(abnormType))))
    #print(paste("barPos:",as.character(barPos1),as.character(barPos2)))
    colorList <- c(colorList,"white")                 # assigns a color to the relative position on the color list
    ##### For full length abnormaltities, the color is set and then the range is defined as all positions
    if (abnormType == "tri" | abnormType == "mono" | abnormType == "der" | abnormType == "blank") { 
      if (abnormType == "tri") {
        colorList <- c(colorList,10000)
      } else if (abnormType == "mono") {
        colorList <- c(colorList, 20000)  
      } else if (abnormType == "der") {
        colorList <- c(colorList, 60000) 
      } else if (abnormType == "blank") {
        colorList <- c(colorList, "white")
      }
      beforeBar <- chromStart
      barLength <- chromEnd-chromStart
      afterBar <- length(gr)-chromEnd
      #print(barLength)
      ##### For greater/less than abnorms, color is set, then the bar is defined and the range is set to include all above or below
    } else if (abnormType == "trans" | abnormType == "del" | abnormType == "add" | abnormType == "inv" | abnormType == "dup") {                  # If the abnormality is a translocation or a deletion, process as such
      # Colors are set
      if        (abnormType == "trans") {colorList <- c(colorList, 50000)
      } else if (abnormType == "del")   {colorList <- c(colorList, 40000)             
      } else if (abnormType == "inv")   {colorList <- c(colorList, 60000)
      } else if (abnormType == "dup")   {colorList <- c(colorList, 30000)}
      if (grepl("p",band1) == TRUE) {        # generates bar if on p arm
        barPos1 <- max(barPos1)
        arm <- "p"
      } else if (grepl("q",band1) == TRUE) {
        barPos1 <- min(barPos1)
        arm <- "q"
      }
      if (band2 != "-") {
        if (grepl("p",band2) == TRUE) {
          barPos2 <- min(barPos2)
        } else if (grepl("q",band2) == TRUE) {
          barPos2 <- max(barPos2)
        }
      } else {
        if (arm == "p") {
          barPos2 <- 1
        } else if (arm == "q") {
          barPos2 <- chromLength
        }
      }
      barLength <- abs(barPos1-barPos2) +1
      preBar <- min(c(barPos1,barPos2)) -1
      postBar <- chromLength - max(c(barPos1,barPos2))
      beforeBar <- chromStart + preBar
      afterBar <- length(gr) - chromEnd + postBar 
      print(paste(c("barPos1:",barPos1,"barPos2:",barPos2,beforeBar,barLength,afterBar)))
    }
    barRange <- Rle(c(0,barX,0),c(beforeBar,barLength,afterBar))
    #print(c(beforeBar,barLength,afterBar))
    blankBar <- Rle(c(0,barX-scaleFactor+0.2,0),c(beforeBar,barLength,afterBar))
    dataList[[length(dataList)+1]] <- blankBar
    dataList[[length(dataList)+1]] <- barRange
    dataNames <- c(dataNames,paste0(c("","Blank_"),as.character(chromosomes[h]),"-",as.character(i)))
    
  }
  return(list(dataList,dataNames,colorList))
}

### Make a Genomic Range Object #############################################################
##### This GRO specifies the chromosomes to display, the positions of the bands, and the data

chromLength <- c()
Range <- c()
Strand <- c()

setwd("C:\\Users/grossar/Bioinform//DATA/ideogram_data/")

file <- "Fib_abnormalitites.txt" ; scaleFactor <- 0.75
file <- "PBMC_abnormalities.txt" ; scaleFactor <- 0.948
file <- "LCL_abnormalities.txt" ; scaleFactor <- 3
file <- "Epithelial_abnormalities.txt" ; scaleFactor <- 7.385
file <- "iPSCs_abnormalities.txt" ; scaleFactor <- 7.385
file <- "Legend2.txt" ; scaleFactor <- 1

abnormFull<-read.table(file)
colorSet<-c("#006837", "#A50026", "#1A9850", "#D73027", "#313695", "#4575B4")  # Final choice

maxY <- 12
chromosomes <- levels(abnormFull[,1])
abnorm<-abnormFull
chromosomes <- levels(abnorm[,1])
gr<-grFormat(chromosomes)

### Formatting data columns ###################################################################
### Specifying chromosomal positiosn by band ##################################################

dataList <- list()
dataNames <- c()
colorList <- c()

for (h in 1:length(chromosomes)) {                    ### runs through the number of chromosomes in the list
  chromStart <- start(seqnames(gr))[h] -1               # defines starting value of the current chr within the full GRO
  chromEnd <- end(seqnames(gr))[h]                    # defines ending value of the current chr within the full GRO
  chromLength <- chromEnd-chromStart
  grSub <- window(gr,chromStart+1,chromEnd)             # generates subset of the ideo data from the starting and ending values
  abnormSub <- data.frame(split(abnorm,abnorm$V1)[h]) # creates a dataframe made by splitting the original data by chrom
  print(paste("# abnormSub:",as.character(levels(seqnames(gr))[h])))
  print(paste(c("chromLength:",chromLength)))
  output <- chromFormat(chromStart,chromEnd,chromLength,grSub,abnormSub)
  dataList <- output[[1]]
  dataNames <- output[[2]]
  colorList <- output[[3]]
}

names(dataList) <- dataNames                        # assigns names to what will become columns
dataList <- rev(dataList)                           # reverses list so right-most column will be behind other bands
colorList <- rev(colorList)

### Color assigner #######################################################################
##### Unspecified values get specific colors assigned to them

for (j in 1:length(colorSet)) {
  print(paste("color #",as.character(j)))
  #print(colorSet[j])
  print(grep(j*10000,colorList))
  print(colorSet[j])
  colorList<-replace(colorList,grep(j*10000,colorList),colorSet[j])
}

mcols(gr) <- data.frame(dataList)                   # appends finished columns to GRO object of chromosome of interest

### Plotting GRO + data #######################################################################
##### The following commands plot the data specified in the GRO and format them as desired

plotOnIdeo(chrom=seqlevels(gr),            # which chrom to plot?
           ideoTable=ideo,                  # ideogram name
           values_GR=gr,                   # data goes here
           value_cols=colnames(mcols(gr)), # col to plot
           col= colorList,
           addScale = F,
           val_range=c(0,maxY),                # set y-axis range
           plot_title=paste(file),
           plotType='rect',  # plot as bars
           vertical=T)

#g + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())

### Output instructions ###################################################################

# Once generated, output as a png of 5400 px wide, 800 tall.
