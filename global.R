library(shiny)
library(shinyalert)
library(shinythemes)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyTree)
library(shinyWidgets)
library(data.table)
library(DT)
library(dplyr)

library(TnT)

library(biovizBase)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
options(scipen = 999)
options(shinyTree.refresh = TRUE)

simpleDebug = TRUE
clickGene <- FALSE
clickPos <- FALSE

# In order to improve the shiny load time, the data and functions are precomputed and stored into a RData file
load("/shares/CIBIO-Storage/BCGLAB/CONREL/shinyApp/CONREL/global.RData")

# Ensembl data
seqlevelsStyle(EnsDb.Hsapiens.v75) <- "UCSC"
genesEns = as.data.table(genes(EnsDb.Hsapiens.v75)$symbol)
allEns = genes(EnsDb.Hsapiens.v75,columns="gene_biotype",return.type="DataFrame")$gene_biotype
levelsColor = unique(allEns)
#mapColor <- c(grDevices::rainbow(10),grDevices::rainbow(10,start = 0.2),grDevices::rainbow(length(levelsColor)-20,start = 0.4))
mapColor <- grDevices::rainbow(length(levelsColor))
c25 <- c(
  "#1c86ee", "#E31A1C", "#008b00", "#6A3D9A","#FF7F00", "#d4d4d4", "#ffd700", "#7ec0ee", "#FB9A99", "#90ee90", "#CAB2D6","#FDBF6F",
  "#b3b3b3", "#eee685", "#b03060", "#ff83fa", "#ff1493", "#0000ff", "#36648b", "#00ced1", "#00ff00", "#8b8b00", "#cdcd00", "#8b4500", "#a52a2a"
)
mapColor = c(c25,c25)[1:length(levelsColor)]
table(allEns)
tabEns = table(allEns)[match(levelsColor,names(table(allEns)))]
legend_idx = order(tabEns,decreasing = T)

colnames(genesEns) = c("gene_symbol")
chrFilter = AnnotationFilterList(SeqNameFilter("chr1"),SeqNameFilter("chr2"),SeqNameFilter("chr3"),SeqNameFilter("chr4"),
                                 SeqNameFilter("chr5"),SeqNameFilter("chr6"),SeqNameFilter("chr7"),SeqNameFilter("chr8"),
                                 SeqNameFilter("chr9"),SeqNameFilter("chr10"),SeqNameFilter("chr11"),SeqNameFilter("chr12"),
                                 SeqNameFilter("chr13"),SeqNameFilter("chr14"),SeqNameFilter("chr15"),SeqNameFilter("chr16"),
                                 SeqNameFilter("chr17"),SeqNameFilter("chr18"),SeqNameFilter("chr19"),SeqNameFilter("chr20"),
                                 SeqNameFilter("chr21"),SeqNameFilter("chr22"),SeqNameFilter("chrX"),SeqNameFilter("chrY"),
                                 SeqNameFilter("chrM"),logicOp = rep(c("|"),24))
