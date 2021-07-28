library(shiny)
library(shinyalert)
library(shinythemes)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(waiter)

library(shinyTree)
library(shinyWidgets)
library(data.table)
library(DT)
library(dplyr)

library(TnT)

# library(biovizBase)
library(AnnotationFilter)
library(GenomicFeatures)

options(scipen = 999)
options(shinyTree.refresh = TRUE)

simpleDebug = TRUE
clickGene <- FALSE
clickPos <- FALSE

# In order to improve the shiny load time, the data and functions are precomputed and stored into a RData file
### moved in chooseORG
# load("./global_v2021.RData")
# data.frame for tooltip. It is empty at the beginning because otherwise an error will be print
df_tooltip = NULL
df_tba = NULL
df_cell = NULL

trackOut = TnT::BlockTrack(GRanges("chr1",IRanges(0,1)),
                           color = "#EEEEEE",background = "#EEEEEE",
                           height = 15,label=NULL)