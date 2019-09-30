library(shiny)
library(shinyalert)
library(shinythemes)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyTree)
library(shinyWidgets)
library(data.table)
library(TnT)
library(biovizBase)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
options(scipen = 999)

simpleDebug = TRUE
clickGene <- FALSE
clickPos <- FALSE


annot <- new.env()
load("/shares/CIBIO-Storage/CIBIO/sharedRL/Projects/MotifsData/MasterMotifs_20181017.RData",envir = annot)
# DEFAULT
# encodeFolder = "/shares/CIBIO-Storage/CO/elaborazioni/sharedCO/Characterization_SNP1376350/ENCODE/"
encodeFolder = "/shares/CIBIO-Storage/BCGLAB/ddalfovo/consensus/"
TBA_folder = "/shares/CIBIO-Storage/BCGLAB/ddalfovo/"
dasServer="http://genome.ucsc.edu/cgi-bin/das/hg19/dna"
defaultGene = "BRCA2"
elements <- c("chr7","139,424,940","141,784,100")
hg19 = fread(paste0(TBA_folder,"hg19.chrom.bed"),data.table = F)

# Ensembl data
seqlevelsStyle(EnsDb.Hsapiens.v75) <- "UCSC"
genesEns = as.data.table(genes(EnsDb.Hsapiens.v75)$symbol)
colnames(genesEns) = c("gene_symbol")
chrFilter = AnnotationFilterList(SeqNameFilter("chr1"),SeqNameFilter("chr2"),SeqNameFilter("chr3"),SeqNameFilter("chr4"),
                                 SeqNameFilter("chr5"),SeqNameFilter("chr6"),SeqNameFilter("chr7"),SeqNameFilter("chr8"),
                                 SeqNameFilter("chr9"),SeqNameFilter("chr10"),SeqNameFilter("chr11"),SeqNameFilter("chr12"),
                                 SeqNameFilter("chr13"),SeqNameFilter("chr14"),SeqNameFilter("chr15"),SeqNameFilter("chr16"),
                                 SeqNameFilter("chr17"),SeqNameFilter("chr18"),SeqNameFilter("chr19"),SeqNameFilter("chr20"),
                                 SeqNameFilter("chr21"),SeqNameFilter("chr22"),SeqNameFilter("chrX"),SeqNameFilter("chrY"),
                                 SeqNameFilter("chrM"),logicOp = rep(c("|"),24))

# Increase the windows by 1M before end after the selection (if it is possible, or max at the chr limits)
increaseWindow = 1000000
window2Load <- function(elements){
  c(if(as.numeric(elements[2])-increaseWindow>0){as.numeric(elements[2])-increaseWindow}else{1},
    if(as.numeric(elements[3])+increaseWindow>as.numeric(hg19[hg19$V1==elements[1],3])){as.numeric(hg19[hg19$V1==elements[1],3])}else{as.numeric(elements[3])+increaseWindow})
}

# data.frame for tooltip. It is empty at the beginning beacause otherwise an error will be print
df_tooltip = NULL
trackOut = TnT::BlockTrack(GRanges("chr1",IRanges(0,1)),
                           color = "#EEEEEE",background = "#EEEEEE",
                           height = 15,label=NULL)
# It converts the input regions in numeric removing the commas
convertPosition <- function(input,v=TRUE){
  if(v){
    input <- strsplit(input,":|-")[[1]]
  }
  input[2] <- as.numeric(gsub(",","",input[2]))
  input[3] <- as.numeric(gsub(",","",input[3]))
  return(input)
}

source("scripts/validate.R")
source("scripts/consensus.R")

generateEmptyTrack <- function(chr){
  return(TnT::BlockTrack(GRanges(chr,IRanges(0,1)),
                               color = "#EEEEEE",background = "#EEEEEE",
                               height = 15,label=NULL))
}

consensusColor <- function(consensus,peak){
  if(consensus=='promoter'){
    if(peak=='narrow'){
      return('darkorange')
    } else {
      return('darkorange4')
    } 
  } else if(consensus=='enhancer'){
    if(peak=='narrow'){
      return('royalblue')
    } else {
      return('royalblue4')
    } 
  } else if(consensus=='active'){
    if(peak=='narrow'){
      return('springgreen')
    } else {
      return('springgreen4')
    } 
  } else {
    return('black')
  }
}

list.depth <- function(this, thisdepth = 0) {
  # http://stackoverflow.com/a/13433689/1270695
  if(!is.list(this)) {
    return(thisdepth)
  } else {
    return(max(unlist(lapply(this, list.depth, thisdepth = thisdepth+1))))    
  }
}

### COUNT PFM ###
load("/shares/CIBIO-Storage/BCGLAB/ddalfovo/countPFMS.RData")
minCount = 100
# Load all the tissue consensus available. NOT USED, instead I decided to use the cell line consensus mapping
# tissueDF = tissueAvailability(tissue,c('narrow','broad'),c('promoter','active','enhancer'))
# tissueDF$tissue = unlist(lapply(tissueDF$listFiles,function(x)strsplit(as.character(x),"\\.")[[1]][1]))

#######
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# elements=c("chr7","139424940","141784100")

# Useful for searching by genes
# TnTGenome(ensGeneTrack, view.range = gene[gene$symbol == "BRCA2"][1] * .7)

loadSNPs <- function(elements,window_load){
  snps.filtered = system(paste0("tabix ",TBA_folder,"dbSNP_v151_hg19/dbSNP_151_hg19_complete_TOPMED_ALL.bed.gz ",gsub("chr","",elements[1]),":",window_load[1],"-",window_load[2]),intern = T)
  snps.filtered = data.frame(do.call(rbind, strsplit(snps.filtered, "\t", fixed=TRUE)))
  
  ir <- IRanges(start = as.numeric(as.character(snps.filtered$X2)), width = 1)
  gpos <- GRanges(elements[1], ir)
  gpos$score <- rep(1,nrow(snps.filtered))
  snps.filtered$alleles = paste(snps.filtered$X4,snps.filtered$X5,sep=">")
  gtooltip <- snps.filtered[,c(3,6,7,8)]
  colnames(gtooltip) = c("rs ID", "1000G MAF", "TOPMED MAF", "Alleles")
  pt = PinTrack(gpos,value = gpos$score,color = "black",label = "dbSNP v151", domain = c(0,1), tooltip = gtooltip,background = "#eeeeee")
  emptyTrack = TnT::BlockTrack(GRanges(elements[1],IRanges(0,1)),
                               color = "#EEEEEE",background = "#EEEEEE",
                               height = 15,label=NULL)
  return(append(pt,emptyTrack))
}


# output$messageMenu <- renderMenu({
#   dropdownMenu(type = "messages",
#                .list = lapply(validateInput(input$peaks,input$regElement,input$consensus,input$tissue),
#                               function(x)messageItem(from="Admin",message = x)))
# })