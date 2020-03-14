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
colnames(genesEns) = c("gene_symbol")
chrFilter = AnnotationFilterList(SeqNameFilter("chr1"),SeqNameFilter("chr2"),SeqNameFilter("chr3"),SeqNameFilter("chr4"),
                                 SeqNameFilter("chr5"),SeqNameFilter("chr6"),SeqNameFilter("chr7"),SeqNameFilter("chr8"),
                                 SeqNameFilter("chr9"),SeqNameFilter("chr10"),SeqNameFilter("chr11"),SeqNameFilter("chr12"),
                                 SeqNameFilter("chr13"),SeqNameFilter("chr14"),SeqNameFilter("chr15"),SeqNameFilter("chr16"),
                                 SeqNameFilter("chr17"),SeqNameFilter("chr18"),SeqNameFilter("chr19"),SeqNameFilter("chr20"),
                                 SeqNameFilter("chr21"),SeqNameFilter("chr22"),SeqNameFilter("chrX"),SeqNameFilter("chrY"),
                                 SeqNameFilter("chrM"),logicOp = rep(c("|"),24))
