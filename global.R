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

# library(TnT)

# library(biovizBase)
# library(AnnotationFilter)
# library(GenomicFeatures)

options(scipen = 999)
options(shinyTree.refresh = TRUE)

simpleDebug = TRUE
clickGene <- FALSE
clickPos <- FALSE

if(simpleDebug){cat(file=stderr(), paste0('start globalFile\n'))}


# In order to improve the shiny load time, the data and functions are precomputed and stored into a RData file
### moved in chooseORG
# load("./global_v2021.RData")
# data.frame for tooltip. It is empty at the beginning because otherwise an error will be print
df_tooltip = NULL
df_tba = NULL
df_cell = NULL

trackOut = TnT::BlockTrack(GenomicRanges::GRanges("chr1",IRanges::IRanges(0,1)),
                           color = "#EEEEEE",background = "#EEEEEE",
                           height = 15,label=NULL)

if(simpleDebug){cat(file=stderr(), paste0('step1 globalFile\n'))}

# a function to check the presence of the folder with data, otherwise the buttons are disabled
checkFolderAssembly <- function() {
  if(!dir.exists('/CONREL_Data/hg19')){disable('hg19')}
  if(!dir.exists('/CONREL_Data/hg38')){disable('hg38')}
  if(!dir.exists('/CONREL_Data/mm10')){disable('mm10')}
}

if(simpleDebug){cat(file=stderr(), paste0('step2 globalFile\n'))}

dataModal <- function(failed = FALSE) {
  modalDialog(title = "",
              div(
                id = "main_content",
                fluidRow(
                  column(12,align='center',
                         p('Welcome in CONREL - CONsensus Regulatory ELement.',style = 'font-size: 21px;color: dimgrey;font-family: "Georgia", Times, "Times New Roman", serif;font-weight: bold;'),
                         p('',style = 'font-size: 21px;color: dimgrey;font-family: "Georgia", Times, "Times New Roman", serif;font-weight: bold;')
                  )
                ),
                # fluidRow(
                #   column(8,offset=2,
                #          fluidRow(
                #            div(class='borderd-content',style='height: 300px;',
                #                div(class='title','Human (Homo sapiens)'),
                #                div(class='content',
                #                    column(12,actionButton('hg38','',width = '100%',
                #                                           style = "color:red;font-weight: bold;width: 120px; height: 120px;background: url('human38.png');  background-size: cover; background-position: center;")),
                #                    column(12,actionButton('hg19','',width='100%',
                #                                           style = "color:red;font-weight: bold;width: 120px; height: 120px;background: url('human19.png');  background-size: cover; background-position: center;"))
                #                )),
                #            div(class='borderd-content',style='height: 150px;',
                #                div(class='title','Mouse (Mus musculus)'),
                #                div(class='content',
                #                    column(12,actionButton('mm10','',width='100%',
                #                                           style = "color:red;font-weight: bold;width: 120px; height: 120px;background: url('mouse.png');  background-size: cover; background-position: center;")),
                #                ))
                #          ))
                # )
                column(10,offset=1,userBox(width = 12,
                                           id = "userbox",
                                           title = userDescription(
                                             title = "Human",
                                             subtitle = "Homo sapiens",
                                             type = 2,
                                             image = "david.png",
                                           ),
                                           status = "navy",
                                           gradient = TRUE,
                                           background = "olive",
                                           boxToolSize = "sm",collapsible = F,closable = F,collapsed = F,
                                           footer=fluidRow(column(10,offset=1,
                                                                  actionButton('hg38','GRCh38/hg38',width = '100%',style='text-transform: none !important;'),
                                                                  actionButton('hg19','GRCh37/hg19',width = '100%',style='text-transform: none !important;')))
                )),
                column(10,offset=1,userBox(width = 12,
                                           id = "userbox",
                                           title = userDescription(
                                             title = "Mouse",
                                             subtitle = "Mus musculus",
                                             type = 2,
                                             image = "mus.png",
                                           ),
                                           status = "navy",
                                           gradient = TRUE,
                                           background = "orange",
                                           boxToolSize = "sm",collapsible = F,closable = F,collapsed = F,
                                           footer=fluidRow(column(10,offset=1,
                                                                  actionButton('mm10','GRCm38/mm10',width = '100%',style='text-transform: none !important;')))
                ))

              ),
              div(
                id = "loading_page",
                fluidRow(align='center',
                         span('Please wait, loading data....'),hr(),
                         spin_2()
                ),br(),
                fluidRow(column(2,offset=2,
                                div(id='ok1',icon("check"))),
                         column(8,span('Preparation parameters and functions')
                         )
                ),br(),
                fluidRow(column(2,offset=2,
                                div(id='ok2',icon("check"))),
                         column(8,span('Load files mapping cell-lines, tissues, etc...')
                         )
                ),br(),
                fluidRow(column(2,offset=2,
                                div(id='ok3',icon("check"))),
                         column(8,span('Load SNPs, TSSs, TBAs, etc..')
                         )
                ),br(),
                fluidRow(column(2,offset=2,
                                div(id='ok4',icon("check"))),
                         column(8,span('Loading Ensembl data')
                         )
                )
              ),
              if (failed)
                div(tags$b("An error has occured", style = "color: red;")),
              
              footer = tagList(
              )
  )
}

if(simpleDebug){cat(file=stderr(), paste0('step3 globalFile\n'))}

###################################################################################
# With global <<- assignment
loadAssemblyData_old <- function(assembly,session,mouse=F) {
  if(simpleDebug){cat(file=stderr(), paste0("Load...",assembly,'\n'))}
  
  removeUI(
    ## pass in appropriate div id
    selector = '#logo-lg'
  )
  removeUI(
    ## pass in appropriate div id
    selector = '#logo-sq'
  )
  
  if(simpleDebug){cat(file=stderr(), paste0('Preparation parameters and functions\n'))}
  shinyjs::show("ok1")
  
  base_path <<- paste0("/CONREL_Data/",assembly,"/")
  encodeFolder <<- paste0("/CONREL_Data/",assembly,"/consensus/")
  inputFolder <<- paste0("/CONREL_Data/",assembly,"/")
  TBA_folder <<- paste0("/CONREL_Data/",assembly,"/")
  hg19 <<- fread(paste0(inputFolder,"data/",assembly,".chrom.bed"),data.table = F)
  hg_assembly <<- assembly
  
  if(assembly=='hg19') {
    EnsDb <<- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    regionEx1 <<- 'chrX:66,764,465-66,950,461'
    regionEx2 <<- 'chr19:51,353,417-51,367,543'
  } else if(assembly=='hg38') {
    EnsDb <<- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    regionEx1 <<- 'chrX:67,545,147-67,694,722'
    regionEx2 <<- 'chr19:50,849,918-50,863,080'
  } else if(assembly=='mm10') {
    EnsDb <<- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
    regionEx1 <<- 'chrX:98,149,782-98,316,880'
    regionEx2 <<- 'chr6:125,157,240-125,171,081'
    
  }
  
  
  # In order to improve the shiny load time, the data and functions are precomputed and stored into a RData file
  # load("./global_v2021.RData")
  ###############
  encodeFolder <<- paste0(base_path,"consensus/")
  TBA_folder <<- base_path
  defaultGene <<- "BRCA2"
  elements <<- c("chr7","139,424,940","141,784,100")
  
  annot <<- new.env()
  load(paste0(inputFolder,"data/MasterMotifs_20181017.RData"),envir = annot)
  
  # Increase the windows by 1M before end after the selection (if it is possible, or max at the chr limits)
  increaseWindow <<- 1000000
  window2Load <<- function(elements){
    c(if(as.numeric(elements[2])-increaseWindow>0){as.numeric(elements[2])-increaseWindow}else{1},
      if(as.numeric(elements[3])+increaseWindow>as.numeric(hg19[hg19$V1==elements[1],3])){as.numeric(hg19[hg19$V1==elements[1],3])}else{as.numeric(elements[3])+increaseWindow})
  }
  
  
  # It converts the input regions in numeric removing the commas
  convertPosition <<- function(input,v=TRUE){
    if(v){
      input <- strsplit(input,":|-")[[1]]
    }
    input[2] <- as.numeric(gsub(",","",input[2]))
    input[3] <- as.numeric(gsub(",","",input[3]))
    return(input)
  }
  
  
  # Check the availability of tissue consensus for narrow and broad peak
  tissueAvailability <<- function(tissue,peaks,regElement) {
    r = lapply(peaks,function(p){
      regElement=gsub("active enhancer","active",regElement)
      listFiles = list.files(paste0(encodeFolder,p,"Peaks/bgz_tissueConsensus"),
                             pattern = paste(paste0("(",paste(tissue,collapse = "|"),")"),
                                             paste0("(",paste(regElement,collapse = "|"),")"),
                                             "union",
                                             sep = "."))
      listFiles = listFiles[-grep(".tbi",listFiles)]
      if(length(listFiles)>0){
        data.frame(p,listFiles)
      }
    })
    return(as.data.frame(do.call(rbind,r)))
  }
  
  # The function returns the name of the cell-line selected in the tree
  extractCellline <<- function(cellList){
    unlist(lapply(cellList,function(x){
      if(list.depth(x)==3) {
        names(x[[1]][[1]])
      }
    }))
  }
  
  # Check the availability of tissue consensus for narrow and broad peak
  celllineAvailability <<- function(cellList,peaks,regElement) {
    cellline = extractCellline(cellList)
    r = lapply(peaks,function(p){
      regElement=gsub("active enhancer","active",regElement)
      listFiles = list.files(paste0(encodeFolder,p,"Peaks/bgz_consensus"),
                             pattern = paste(paste0("(",paste(cellline,collapse = "|"),")"),
                                             "intersect",
                                             paste0("(",paste(regElement,collapse = "|"),")"),
                                             sep = "."))
      listFiles = listFiles[-grep(".tbi",listFiles)]
      if(length(listFiles)>0){
        data.frame(p,listFiles)
      }
    })
    return(as.data.frame(do.call(rbind,r)))
  }
  
  # The function returns the files (combination of, in a DF) for global consensus selected by user
  mappingInput <<- function(consensus,peaks) {
    consensus = gsub("active enhancer","active",consensus)
    combFolder = expand.grid(peaks,consensus)
    return(combFolder)
  }
  
  if(simpleDebug){cat(file=stderr(), paste0('Load files mapping cell-lines, tissues, etc...\n'))}
  shinyjs::show("ok2")
  # TISSUES consensus for all ENCODE (the list in the bcglab DB)
  mapping <<- fread(paste0(inputFolder,"consensus/mapping_cellLines.csv"),data.table = F)
  mapping$name <<- gsub("'","",mapping$name)
  mapping$name <<- gsub(" ","_",mapping$name)
  mapping$name <<- gsub("/","_",mapping$name)
  mapping$tissue <<- gsub(" ","_",mapping$tissue)
  
  tissue <<- sort(unique(mapping$tissue))
  tissueOG <<- tissue[-grep("/",tissue)]
  
  # CELL-LINES consensus (list of cell-lines available) for narrow and broad peaks
  narrow_cell <<- fread(paste0(inputFolder,"consensus/narrowPeaks/dataTableNoReplicates.tsv"),data.table = F)
  narrow_cell$Biosample.term.name <<- gsub(" ","_",narrow_cell$Biosample.term.name)
  narrow_cell$tissueMapped <<- mapping[match(narrow_cell$Biosample.term.name,mapping$name),"tissue"]
  narrow_cell <<- unique(narrow_cell[,c("Biosample.term.name","tissueMapped")])
  idx <<- grep("/",narrow_cell$tissueMapped)
  split_multiTissue <<- strsplit(narrow_cell[idx,2],"/")
  res <<- lapply(1:length(split_multiTissue),function(x){
    data.frame(Biosample.term.name=narrow_cell[idx[x],1],tissueMapped=split_multiTissue[[x]])
  })
  narrow_cell <<- rbind(narrow_cell[-idx,],do.call(rbind,res))
  
  broad_cell <<- fread(paste0(inputFolder,"consensus/broadPeaks/dataTableNoReplicates.tsv"),data.table = F)
  broad_cell$Biosample.term.name <<- gsub(" ","_",broad_cell$Biosample.term.name)
  broad_cell$tissueMapped <<- mapping[match(broad_cell$Biosample.term.name,mapping$name),"tissue"]
  broad_cell <<- unique(broad_cell[,c("Biosample.term.name","tissueMapped")])
  idx <<- grep("/",broad_cell$tissueMapped)
  split_multiTissue <<- strsplit(broad_cell[idx,2],"/")
  res <<- lapply(1:length(split_multiTissue),function(x){
    data.frame(Biosample.term.name=broad_cell[idx[x],1],tissueMapped=split_multiTissue[[x]])
  })
  broad_cell <<- rbind(broad_cell[-idx,],do.call(rbind,res))
  
  
  
  generateEmptyTrack <<- function(chr){
    return(TnT::BlockTrack(GenomicRanges::GRanges(chr,IRanges::IRanges(0,1)),
                           color = "#EEEEEE",background = "#EEEEEE",
                           height = 15,label=NULL))
  }
  
  
  
  checkValidRegion <<- function(region) {
    elements <<- convertPosition(region,TRUE)
    window_load <<- window2Load(elements)
    start_geneView <<- window_load[1]
    end_geneView <<- window_load[2]
    chr_geneView <<- elements[1]
    gene <<- ensembldb::genes(EnsDb,
                              filter=AnnotationFilter::AnnotationFilterList(AnnotationFilter::GeneStartFilter(start_geneView,condition = '>'),
                                                                            AnnotationFilter::GeneEndFilter(end_geneView,condition = '<'),
                                                                            AnnotationFilter::SeqNameFilter(chr_geneView),
                                                                            logicOp = c('&','&')))
    return(length(gene)==0)
  }
  
  consensusColor <<- function(consensus,peak){
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
  
  list.depth <<- function(this, thisdepth = 0) {
    # http://stackoverflow.com/a/13433689/1270695
    if(!is.list(this)) {
      return(thisdepth)
    } else {
      return(max(unlist(lapply(this, list.depth, thisdepth = thisdepth+1))))
    }
  }
  
  if(simpleDebug){cat(file=stderr(), paste0('Load SNPs, TSSs, TBAs, etc..\n'))}
  shinyjs::show("ok3")
  ### COUNT PFM ###
  load(paste0(inputFolder,"data/countPFMS_v2.RData"))
  pfms <<- pfms
  countPFMS <<- countPFMS
  minCount <<- 50
  
  
  loadSNPs <<- function(elements,window_load){
    snps.filtered = system(paste0("tabix ",inputFolder,"dbSNP_v151_hg19/dbSNP_151_hg19_complete_TOPMED_ALL.bed.gz ",gsub("chr","",elements[1]),":",window_load[1],"-",window_load[2]),intern = T)
    snps.filtered = data.frame(do.call(rbind, strsplit(snps.filtered, "\t", fixed=TRUE)))
    
    ir <- IRanges::IRanges(start = as.numeric(as.character(snps.filtered$X2)), width = 1)
    gpos <- GenomicRanges::GRanges(elements[1], ir)
    gpos$score <- rep(1,nrow(snps.filtered))
    snps.filtered$alleles = paste(snps.filtered$X4,snps.filtered$X5,sep=">")
    gtooltip <- snps.filtered[,c(3,6,7,8)]
    colnames(gtooltip) = c("rs ID", "1000G MAF", "TOPMED MAF", "Alleles")
    pt = TnT::PinTrack(gpos,value = gpos$score,color = "black",label = "dbSNP v151", domain = c(0,1), tooltip = gtooltip,background = "#eeeeee")
    emptyTrack = TnT::BlockTrack(GenomicRanges::GRanges(elements[1],IRanges::IRanges(0,1)),
                                 color = "#EEEEEE",background = "#EEEEEE",
                                 height = 15,label=NULL)
    return(append(pt,emptyTrack))
  }
  
  
  
  loadTSS <<- function(elements,window_load){
    tss.filtered = system(paste0("tabix ",inputFolder,"data/UCSC_TSS.bed.gz ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
    tss.filtered = data.frame(do.call(rbind, strsplit(tss.filtered, "\t", fixed=TRUE)))
    
    ir <- IRanges::IRanges(start = as.numeric(as.character(tss.filtered$X2)), end = as.numeric(as.character(tss.filtered$X3)))
    gpos <- GenomicRanges::GRanges(elements[1], ir)
    # gpos$score <- rep(1,nrow(tss.filtered))
    # snps.filtered$alleles = paste(snps.filtered$X4,snps.filtered$X5,sep=">")
    gtooltip <- tss.filtered[,c(4,5)]
    colnames(gtooltip) = c("strand","confScore")
    pt = TnT::BlockTrack(gpos,color = "black",label = "SwitchGear TSS", tooltip = gtooltip,background = "#eeeeee",
                    height = 15)
    emptyTrack = TnT::BlockTrack(GenomicRanges::GRanges(elements[1],IRanges::IRanges(0,1)),
                                 color = "#EEEEEE",background = "#EEEEEE",
                                 height = 15,label=NULL)
    return(append(pt,emptyTrack))
  }
  
  
  map_Ncells <<- fread("/CONREL_Data/hg19/data/mapping_Ncells.csv")
  map_Ncells_Global <<- data.table(
    rbind(data.table("peaks"="narrow",
                     "cre"="promoter",
                     "count"=sum(map_Ncells$`Narrow peaks`>=1),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=1,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=1,3])),
          data.table("peaks"="narrow",
                     "cre"="enhancer",
                     "count"=sum(map_Ncells$`Narrow peaks`>=2),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=2,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=2,3])),
          data.table("peaks"="narrow",
                     "cre"="active",
                     "count"=sum(map_Ncells$`Narrow peaks`>=3),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=3,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=3,3])),
          data.table("peaks"="broad",
                     "cre"="promoter",
                     "count"=sum(map_Ncells$`Broad peaks`>=1),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=1,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=1,5])),
          data.table("peaks"="broad",
                     "cre"="enhancer",
                     "count"=sum(map_Ncells$`Broad peaks`>=2),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=2,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=2,5])),
          data.table("peaks"="broad",
                     "cre"="active",
                     "count"=sum(map_Ncells$`Broad peaks`>=3),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=3,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=3,5]))
    ))
  
  map_Ncells_Global <<- map_Ncells_Global[, list(cells=list(.SD)), by = list(peaks,cre,count)]
  for(i in 1:length(map_Ncells_Global$cells)) {
    colnames(map_Ncells_Global$cells[[i]]) = c("cell line","Number of experiments")
  }
  
  map_Ncells_Tissue <<- lapply(unique(map_Ncells$Tissue),function(x){
    rbind(data.table("peaks"="narrow",
                     "cre"="promoter",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=1&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=1&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=1&map_Ncells$Tissue==x,3])),
          data.table("peaks"="narrow",
                     "cre"="enhancer",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=2&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=2&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=2&map_Ncells$Tissue==x,3])),
          data.table("peaks"="narrow",
                     "cre"="active",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=3&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=3&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=3&map_Ncells$Tissue==x,3])),
          data.table("peaks"="broad",
                     "cre"="promoter",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=1&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=1&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=1&map_Ncells$Tissue==x,5])),
          data.table("peaks"="broad",
                     "cre"="enhancer",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=2&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=2&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=2&map_Ncells$Tissue==x,5])),
          data.table("peaks"="broad",
                     "cre"="active",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=3&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=3&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=3&map_Ncells$Tissue==x,5]))
    )})
  map_Ncells_Tissue_bind <<- do.call(rbind,map_Ncells_Tissue)
  map_Ncells_Tissue_bind <<- map_Ncells_Tissue_bind[map_Ncells_Tissue_bind$count!=0,]
  map_Ncells_Tissue_table <<- map_Ncells_Tissue_bind[, list(cells=list(.SD)), by = list(peaks,cre,count,tissue)]
  for(i in 1:length(map_Ncells_Tissue_table$cells)) {
    colnames(map_Ncells_Tissue_table$cells[[i]]) <<- c("cell line","Number of experiments")
  }
  
  map_Ncells_Cell <<- lapply(unique(map_Ncells$`Cell line`),function(x){
    rbind(data.table("peaks"="narrow",
                     "cre"="promoter",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=1&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=1&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=1&map_Ncells$`Cell line`==x,3])),
          data.table("peaks"="narrow",
                     "cre"="enhancer",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=2&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=2&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=2&map_Ncells$`Cell line`==x,3])),
          data.table("peaks"="narrow",
                     "cre"="active",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=3&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=3&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=3&map_Ncells$`Cell line`==x,3])),
          data.table("peaks"="broad",
                     "cre"="promoter",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=1&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=1&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=1&map_Ncells$`Cell line`==x,5])),
          data.table("peaks"="broad",
                     "cre"="enhancer",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=2&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=2&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=2&map_Ncells$`Cell line`==x,5])),
          data.table("peaks"="broad",
                     "cre"="active",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=3&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=3&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=3&map_Ncells$`Cell line`==x,5]))
    )})
  map_Ncells_Cell_bind <<- do.call(rbind,map_Ncells_Cell)
  map_Ncells_Cell_bind <<- map_Ncells_Cell_bind[map_Ncells_Cell_bind$count!=0,]
  map_Ncells_Cell_table <<- map_Ncells_Cell_bind[, list(cells=list(.SD)), by = list(peaks,cre,count,cell)]
  for(i in 1:length(map_Ncells_Cell_table$cells)) {
    colnames(map_Ncells_Cell_table$cells[[i]]) <<- c("cell line","Number of experiments")
  }
  
  ###############
  if(simpleDebug){cat(file=stderr(), paste0('Loading Ensembl data!\n'))}
  shinyjs::show("ok4")
  ensembldb::seqlevelsStyle(EnsDb) <<- "UCSC"
  genesEns <<- suppressWarnings(as.data.table(ensembldb::genes(EnsDb,return.type='DataFrame')))
  genesEns <<- genesEns[genesEns$seq_name%in%paste0('chr',c(1:22,'X','Y','M')),'symbol']
  allEns <<- ensembldb::genes(EnsDb,columns="gene_biotype",return.type="DataFrame")$gene_biotype
  allEns_tx <<- ensembldb::genes(EnsDb,columns="tx_biotype",return.type="DataFrame")$tx_biotype
  levelsColor <<- unique(allEns)
  levelsColor_tx <<- unique(allEns_tx)
  #mapColor <- c(grDevices::rainbow(10),grDevices::rainbow(10,start = 0.2),grDevices::rainbow(length(levelsColor)-20,start = 0.4))
  mapColor <<- grDevices::rainbow(length(levelsColor))
  c25 <<- c(
    "#1c86ee", "#E31A1C", "#008b00", "#6A3D9A","#FF7F00", "#d4d4d4", "#ffd700", "#7ec0ee", "#FB9A99", "#90ee90", "#CAB2D6","#FDBF6F",
    "#b3b3b3", "#eee685", "#b03060", "#ff83fa", "#ff1493", "#0000ff", "#36648b", "#00ced1", "#00ff00", "#8b8b00", "#cdcd00", "#8b4500", "#a52a2a"
  )
  mapColor <<- c(c25,c25)[1:length(levelsColor)]
  mapColor_tx <<- c(c25,c25)[1:length(levelsColor_tx)]
  
  tabEns <<- table(allEns)[match(levelsColor,names(table(allEns)))]
  tabEns_tx <<- table(allEns_tx)[match(levelsColor_tx,names(table(allEns_tx)))]
  legend_idx <<- order(tabEns,decreasing = T)
  legend_idx_tx <<- order(tabEns_tx,decreasing = T)
  
  colnames(genesEns) <<- c("gene_symbol")
  chrFilter <<- AnnotationFilter::AnnotationFilterList(AnnotationFilter::SeqNameFilter("chr1"),AnnotationFilter::SeqNameFilter("chr2"),AnnotationFilter::SeqNameFilter("chr3"),AnnotationFilter::SeqNameFilter("chr4"),
                                     AnnotationFilter::SeqNameFilter("chr5"),AnnotationFilter::SeqNameFilter("chr6"),AnnotationFilter::SeqNameFilter("chr7"),AnnotationFilter::SeqNameFilter("chr8"),
                                     AnnotationFilter::SeqNameFilter("chr9"),AnnotationFilter::SeqNameFilter("chr10"),AnnotationFilter::SeqNameFilter("chr11"),AnnotationFilter::SeqNameFilter("chr12"),
                                     AnnotationFilter::SeqNameFilter("chr13"),AnnotationFilter::SeqNameFilter("chr14"),AnnotationFilter::SeqNameFilter("chr15"),AnnotationFilter::SeqNameFilter("chr16"),
                                     AnnotationFilter::SeqNameFilter("chr17"),AnnotationFilter::SeqNameFilter("chr18"),AnnotationFilter::SeqNameFilter("chr19"),AnnotationFilter::SeqNameFilter("chr20"),
                                     AnnotationFilter::SeqNameFilter("chr21"),AnnotationFilter::SeqNameFilter("chr22"),AnnotationFilter::SeqNameFilter("chrX"),AnnotationFilter::SeqNameFilter("chrY"),
                                     AnnotationFilter::SeqNameFilter("chrM"),logicOp = rep(c("|"),24))
  
  updateSelectizeInput(session, "genes",
                       choices = as.character(genesEns$gene_symbol),
                       server = TRUE,selected = NULL)
  if(mouse){
    square_imm = paste0("mus.svg")
    corner_imm = paste0("mus",assembly,"_corner.png")
  }else{
    square_imm = paste0("david.svg")
    corner_imm = paste0("david",assembly,"_corner.png")
  }
  insertUI(
    selector = ".logo-lg",
    where = "afterEnd",
    ui = tagList(div(img(class = "logo-lg", 
                         id='logo-lg', 
                         src=corner_imm,
                         width = "100%", 
                         height = "auto",
                         style="float: left;display:inline-block;margin-left:-25px;")),
                 img(class = "logo-sq", 
                     id='logo-sq', 
                     src = square_imm)
    )
  )
  # removeModal()
  updateTabItems(session, "sideBar", "home")
  
  if(simpleDebug){cat(file=stderr(), paste0('End LoadAssembly\n'))}
}


######### FUNCTIONS
# Increase the windows by 1M before end after the selection (if it is possible, or max at the chr limits)
window2Load <- function(elements){
  c(if(as.numeric(elements[2])-increaseWindow>0){as.numeric(elements[2])-increaseWindow}else{1},
    if(as.numeric(elements[3])+increaseWindow>as.numeric(hg19[hg19$V1==elements[1],3])){as.numeric(hg19[hg19$V1==elements[1],3])}else{as.numeric(elements[3])+increaseWindow})
}
# It converts the input regions in numeric removing the commas
convertPosition <- function(input,v=TRUE){
  if(v){
    input <- strsplit(input,":|-")[[1]]
  }
  input[2] <- as.numeric(gsub(",","",input[2]))
  input[3] <- as.numeric(gsub(",","",input[3]))
  return(input)
}
# Check the availability of tissue consensus for narrow and broad peak
tissueAvailability <- function(tissue,peaks,regElement) {
  r = lapply(peaks,function(p){
    regElement=gsub("active enhancer","active",regElement)
    listFiles = list.files(paste0(encodeFolder,p,"Peaks/bgz_tissueConsensus"),
                           pattern = paste(paste0("(",paste(tissue,collapse = "|"),")"),
                                           paste0("(",paste(regElement,collapse = "|"),")"),
                                           "union",
                                           sep = "."))
    listFiles = listFiles[-grep(".tbi",listFiles)]
    if(length(listFiles)>0){
      data.frame(p,listFiles)
    }
  })
  return(as.data.frame(do.call(rbind,r)))
}
# The function returns the name of the cell-line selected in the tree
extractCellline <- function(cellList){
  unlist(lapply(cellList,function(x){
    if(list.depth(x)==3) {
      names(x[[1]][[1]])
    }
  }))
}
# Check the availability of tissue consensus for narrow and broad peak
celllineAvailability <- function(cellList,peaks,regElement) {
  cellline = extractCellline(cellList)
  r = lapply(peaks,function(p){
    regElement=gsub("active enhancer","active",regElement)
    listFiles = list.files(paste0(encodeFolder,p,"Peaks/bgz_consensus"),
                           pattern = paste(paste0("(",paste(cellline,collapse = "|"),")"),
                                           "intersect",
                                           paste0("(",paste(regElement,collapse = "|"),")"),
                                           sep = "."))
    listFiles = listFiles[-grep(".tbi",listFiles)]
    if(length(listFiles)>0){
      data.frame(p,listFiles)
    }
  })
  return(as.data.frame(do.call(rbind,r)))
}
# The function returns the files (combination of, in a DF) for global consensus selected by user
mappingInput <- function(consensus,peaks) {
  consensus = gsub("active enhancer","active",consensus)
  combFolder = expand.grid(peaks,consensus)
  return(combFolder)
}
# Empty track for TnT library
generateEmptyTrack <- function(chr){
  return(TnT::BlockTrack(GenomicRanges::GRanges(chr,IRanges::IRanges(0,1)),
                         color = "#EEEEEE",background = "#EEEEEE",
                         height = 15,label=NULL))
}
# Checks the format validity of the regions searched
checkValidRegion <- function(region) {
  elements <- convertPosition(region,TRUE)
  window_load <- window2Load(elements)
  start_geneView <- window_load[1]
  end_geneView <- window_load[2]
  chr_geneView <- elements[1]
  gene <- ensembldb::genes(EnsDb,
                            filter=AnnotationFilter::AnnotationFilterList(AnnotationFilter::GeneStartFilter(start_geneView,condition = '>'),
                                                                          AnnotationFilter::GeneEndFilter(end_geneView,condition = '<'),
                                                                          AnnotationFilter::SeqNameFilter(chr_geneView),
                                                                          logicOp = c('&','&')))
  return(length(gene)==0)
}
# Determines the colors of the regions
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
# returns the depth of a nested list
list.depth <- function(this, thisdepth = 0) {
  # http://stackoverflow.com/a/13433689/1270695
  if(!is.list(this)) {
    return(thisdepth)
  } else {
    return(max(unlist(lapply(this, list.depth, thisdepth = thisdepth+1))))
  }
}
# Load the SNPs data
loadSNPs <- function(elements,window_load){
  snps.filtered = system(paste0("tabix ",inputFolder,"dbSNP_v151_hg19/dbSNP_151_hg19_complete_TOPMED_ALL.bed.gz ",gsub("chr","",elements[1]),":",window_load[1],"-",window_load[2]),intern = T)
  snps.filtered = data.frame(do.call(rbind, strsplit(snps.filtered, "\t", fixed=TRUE)))
  
  ir <- IRanges::IRanges(start = as.numeric(as.character(snps.filtered$X2)), width = 1)
  gpos <- GenomicRanges::GRanges(elements[1], ir)
  gpos$score <- rep(1,nrow(snps.filtered))
  snps.filtered$alleles = paste(snps.filtered$X4,snps.filtered$X5,sep=">")
  gtooltip <- snps.filtered[,c(3,6,7,8)]
  colnames(gtooltip) = c("rs ID", "1000G MAF", "TOPMED MAF", "Alleles")
  pt = TnT::PinTrack(gpos,value = gpos$score,color = "black",label = "dbSNP v151", domain = c(0,1), tooltip = gtooltip,background = "#eeeeee")
  emptyTrack = TnT::BlockTrack(GenomicRanges::GRanges(elements[1],IRanges::IRanges(0,1)),
                               color = "#EEEEEE",background = "#EEEEEE",
                               height = 15,label=NULL)
  return(append(pt,emptyTrack))
}
# Laod the TSS data
loadTSS <- function(elements,window_load){
  tss.filtered = system(paste0("tabix ",inputFolder,"data/UCSC_TSS.bed.gz ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
  tss.filtered = data.frame(do.call(rbind, strsplit(tss.filtered, "\t", fixed=TRUE)))
  
  ir <- IRanges::IRanges(start = as.numeric(as.character(tss.filtered$X2)), end = as.numeric(as.character(tss.filtered$X3)))
  gpos <- GenomicRanges::GRanges(elements[1], ir)
  # gpos$score <- rep(1,nrow(tss.filtered))
  # snps.filtered$alleles = paste(snps.filtered$X4,snps.filtered$X5,sep=">")
  gtooltip <- tss.filtered[,c(4,5)]
  colnames(gtooltip) = c("strand","confScore")
  pt = TnT::BlockTrack(gpos,color = "black",label = "SwitchGear TSS", tooltip = gtooltip,background = "#eeeeee",
                       height = 15)
  emptyTrack = TnT::BlockTrack(GenomicRanges::GRanges(elements[1],IRanges::IRanges(0,1)),
                               color = "#EEEEEE",background = "#EEEEEE",
                               height = 15,label=NULL)
  return(append(pt,emptyTrack))
}

# With assing function in .GlobalEnv
loadAssemblyData <- function(assembly,session,mouse=F) {
  if(simpleDebug){cat(file=stderr(), paste0("Load...",assembly,'\n'))}

  removeUI(
    ## pass in appropriate div id
    selector = '#logo-lg'
  )
  removeUI(
    ## pass in appropriate div id
    selector = '#logo-sq'
  )
  
  if(simpleDebug){cat(file=stderr(), paste0('Preparation parameters and functions\n'))}
  
  assign('base_path',paste0("/CONREL_Data/",assembly,"/"),.GlobalEnv)
  assign('encodeFolder',paste0("/CONREL_Data/",assembly,"/consensus/"),.GlobalEnv)
  assign('inputFolder',paste0("/CONREL_Data/",assembly,"/"),.GlobalEnv)
  assign('TBA_folder',paste0("/CONREL_Data/",assembly,"/"),.GlobalEnv)
  assign('hg19',fread(paste0(inputFolder,"data/",assembly,".chrom.bed"),data.table = F),.GlobalEnv)
  assign('hg_assembly',assembly,.GlobalEnv)
  assign('encodeFolder',paste0(base_path,"consensus/"),.GlobalEnv)
  assign('TBA_folder',base_path,.GlobalEnv)
  assign('defaultGene',"BRCA2",.GlobalEnv)
  assign('elements',c("chr7","139,424,940","141,784,100"),.GlobalEnv)
  
  assign('annot',new.env(),.GlobalEnv)
  load(paste0(inputFolder,"data/MasterMotifs_20181017.RData"),envir = annot)
  
  
  if(assembly=='hg19') {
    EnsDb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    assign('regionEx1','chrX:66,764,465-66,950,461',.GlobalEnv)
    assign('regionEx2','chr19:51,353,417-51,367,543',.GlobalEnv)
  } else if(assembly=='hg38') {
    EnsDb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    assign('regionEx1','chrX:67,545,147-67,694,722',.GlobalEnv)
    assign('regionEx2','chr19:50,849,918-50,863,080',.GlobalEnv)
  } else if(assembly=='mm10') {
    EnsDb <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
    assign('regionEx1','chrX:98,149,782-98,316,880',.GlobalEnv)
    assign('regionEx2','chr6:125,157,240-125,171,081',.GlobalEnv)
  }

  assign('increaseWindow',1000000,.GlobalEnv)
  
  shinyjs::show("ok1") 
  if(simpleDebug){cat(file=stderr(), paste0('Load files mapping cell-lines, tissues, etc...\n'))}

  # TISSUES consensus for all ENCODE (the list in the bcglab DB)
  mapping <- fread(paste0(inputFolder,"consensus/mapping_cellLines.csv"),data.table = F)
  mapping$name <- gsub("'","",mapping$name)
  mapping$name <- gsub(" ","_",mapping$name)
  mapping$name <- gsub("/","_",mapping$name)
  mapping$tissue <- gsub(" ","_",mapping$tissue)
  
  assign('mapping',mapping,.GlobalEnv)
  assign('tissue',sort(unique(mapping$tissue)),.GlobalEnv)
  assign('tissueOG',tissue[-grep("/",tissue)],.GlobalEnv)
  
  
  # CELL-LINES consensus (list of cell-lines available) for narrow and broad peaks
  narrow_cell <- fread(paste0(inputFolder,"consensus/narrowPeaks/dataTableNoReplicates.tsv"),data.table = F)
  narrow_cell$Biosample.term.name <- gsub(" ","_",narrow_cell$Biosample.term.name)
  narrow_cell$tissueMapped <- mapping[match(narrow_cell$Biosample.term.name,mapping$name),"tissue"]
  narrow_cell <- unique(narrow_cell[,c("Biosample.term.name","tissueMapped")])
  idx <- grep("/",narrow_cell$tissueMapped)
  split_multiTissue <- strsplit(narrow_cell[idx,2],"/")
  res <- lapply(1:length(split_multiTissue),function(x){
    data.frame(Biosample.term.name=narrow_cell[idx[x],1],tissueMapped=split_multiTissue[[x]])
  })
  narrow_cell <- rbind(narrow_cell[-idx,],do.call(rbind,res))
  assign('narrow_cell',narrow_cell,.GlobalEnv)
  
  broad_cell <- fread(paste0(inputFolder,"consensus/broadPeaks/dataTableNoReplicates.tsv"),data.table = F)
  broad_cell$Biosample.term.name <- gsub(" ","_",broad_cell$Biosample.term.name)
  broad_cell$tissueMapped <- mapping[match(broad_cell$Biosample.term.name,mapping$name),"tissue"]
  broad_cell <- unique(broad_cell[,c("Biosample.term.name","tissueMapped")])
  idx <- grep("/",broad_cell$tissueMapped)
  split_multiTissue <- strsplit(broad_cell[idx,2],"/")
  res <- lapply(1:length(split_multiTissue),function(x){
    data.frame(Biosample.term.name=broad_cell[idx[x],1],tissueMapped=split_multiTissue[[x]])
  })
  broad_cell <- rbind(broad_cell[-idx,],do.call(rbind,res))
  assign('broad_cell',broad_cell,.GlobalEnv)
  
  shinyjs::show("ok2")
  if(simpleDebug){cat(file=stderr(), paste0('Load SNPs, TSSs, TBAs, etc..\n'))}
  
  ### COUNT PFM ###
  load(paste0(inputFolder,"data/countPFMS_v2.RData"),envir = .GlobalEnv)
  assign('minCount',50,.GlobalEnv)
  
  
  map_Ncells <- fread("/CONREL_Data/hg19/data/mapping_Ncells.csv")
  map_Ncells_Global <- data.table(
    rbind(data.table("peaks"="narrow",
                     "cre"="promoter",
                     "count"=sum(map_Ncells$`Narrow peaks`>=1),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=1,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=1,3])),
          data.table("peaks"="narrow",
                     "cre"="enhancer",
                     "count"=sum(map_Ncells$`Narrow peaks`>=2),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=2,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=2,3])),
          data.table("peaks"="narrow",
                     "cre"="active",
                     "count"=sum(map_Ncells$`Narrow peaks`>=3),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=3,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=3,3])),
          data.table("peaks"="broad",
                     "cre"="promoter",
                     "count"=sum(map_Ncells$`Broad peaks`>=1),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=1,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=1,5])),
          data.table("peaks"="broad",
                     "cre"="enhancer",
                     "count"=sum(map_Ncells$`Broad peaks`>=2),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=2,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=2,5])),
          data.table("peaks"="broad",
                     "cre"="active",
                     "count"=sum(map_Ncells$`Broad peaks`>=3),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=3,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=3,5]))
    ))
  
  map_Ncells_Global <- map_Ncells_Global[, list(cells=list(.SD)), by = list(peaks,cre,count)]
  for(i in 1:length(map_Ncells_Global$cells)) {
    colnames(map_Ncells_Global$cells[[i]]) <- c("cell line","Number of experiments")
  }
  assign('map_Ncells_Global',map_Ncells_Global,.GlobalEnv)
  
  map_Ncells_Tissue <- lapply(unique(map_Ncells$Tissue),function(x){
    rbind(data.table("peaks"="narrow",
                     "cre"="promoter",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=1&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=1&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=1&map_Ncells$Tissue==x,3])),
          data.table("peaks"="narrow",
                     "cre"="enhancer",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=2&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=2&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=2&map_Ncells$Tissue==x,3])),
          data.table("peaks"="narrow",
                     "cre"="active",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=3&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=3&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=3&map_Ncells$Tissue==x,3])),
          data.table("peaks"="broad",
                     "cre"="promoter",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=1&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=1&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=1&map_Ncells$Tissue==x,5])),
          data.table("peaks"="broad",
                     "cre"="enhancer",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=2&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=2&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=2&map_Ncells$Tissue==x,5])),
          data.table("peaks"="broad",
                     "cre"="active",
                     "tissue" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=3&map_Ncells$Tissue==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=3&map_Ncells$Tissue==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=3&map_Ncells$Tissue==x,5]))
    )})
  map_Ncells_Tissue_bind <- do.call(rbind,map_Ncells_Tissue)
  map_Ncells_Tissue_bind <- map_Ncells_Tissue_bind[map_Ncells_Tissue_bind$count!=0,]
  map_Ncells_Tissue_table <- map_Ncells_Tissue_bind[, list(cells=list(.SD)), by = list(peaks,cre,count,tissue)]
  for(i in 1:length(map_Ncells_Tissue_table$cells)) {
    colnames(map_Ncells_Tissue_table$cells[[i]]) <- c("cell line","Number of experiments")
  }
  assign('map_Ncells_Tissue_table',map_Ncells_Tissue_table,.GlobalEnv)
  
  map_Ncells_Cell <- lapply(unique(map_Ncells$`Cell line`),function(x){
    rbind(data.table("peaks"="narrow",
                     "cre"="promoter",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=1&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=1&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=1&map_Ncells$`Cell line`==x,3])),
          data.table("peaks"="narrow",
                     "cre"="enhancer",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=2&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=2&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=2&map_Ncells$`Cell line`==x,3])),
          data.table("peaks"="narrow",
                     "cre"="active",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Narrow peaks`>=3&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Narrow peaks`>=3&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Narrow peaks`>=3&map_Ncells$`Cell line`==x,3])),
          data.table("peaks"="broad",
                     "cre"="promoter",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=1&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=1&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=1&map_Ncells$`Cell line`==x,5])),
          data.table("peaks"="broad",
                     "cre"="enhancer",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=2&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=2&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=2&map_Ncells$`Cell line`==x,5])),
          data.table("peaks"="broad",
                     "cre"="active",
                     "cell" = x,
                     "count"=sum(map_Ncells$`Broad peaks`>=3&map_Ncells$`Cell line`==x),
                     "cells"=data.table(map_Ncells[map_Ncells$`Broad peaks`>=3&map_Ncells$`Cell line`==x,1],
                                        map_Ncells[map_Ncells$`Broad peaks`>=3&map_Ncells$`Cell line`==x,5]))
    )})
  map_Ncells_Cell_bind <- do.call(rbind,map_Ncells_Cell)
  map_Ncells_Cell_bind <- map_Ncells_Cell_bind[map_Ncells_Cell_bind$count!=0,]
  map_Ncells_Cell_table <- map_Ncells_Cell_bind[, list(cells=list(.SD)), by = list(peaks,cre,count,cell)]
  for(i in 1:length(map_Ncells_Cell_table$cells)) {
    colnames(map_Ncells_Cell_table$cells[[i]]) <- c("cell line","Number of experiments")
  }
  assign('map_Ncells_Cell_table',map_Ncells_Cell_table,.GlobalEnv)
  
  ###############
  shinyjs::show("ok3")
  if(simpleDebug){cat(file=stderr(), paste0('Loading Ensembl data!\n'))}

  ensembldb::seqlevelsStyle(EnsDb) <- "UCSC"
  assign('EnsDb',EnsDb,.GlobalEnv)
  assign('chrFilter',AnnotationFilter::AnnotationFilterList(AnnotationFilter::SeqNameFilter("chr1"),AnnotationFilter::SeqNameFilter("chr2"),AnnotationFilter::SeqNameFilter("chr3"),AnnotationFilter::SeqNameFilter("chr4"),
                                                       AnnotationFilter::SeqNameFilter("chr5"),AnnotationFilter::SeqNameFilter("chr6"),AnnotationFilter::SeqNameFilter("chr7"),AnnotationFilter::SeqNameFilter("chr8"),
                                                       AnnotationFilter::SeqNameFilter("chr9"),AnnotationFilter::SeqNameFilter("chr10"),AnnotationFilter::SeqNameFilter("chr11"),AnnotationFilter::SeqNameFilter("chr12"),
                                                       AnnotationFilter::SeqNameFilter("chr13"),AnnotationFilter::SeqNameFilter("chr14"),AnnotationFilter::SeqNameFilter("chr15"),AnnotationFilter::SeqNameFilter("chr16"),
                                                       AnnotationFilter::SeqNameFilter("chr17"),AnnotationFilter::SeqNameFilter("chr18"),AnnotationFilter::SeqNameFilter("chr19"),AnnotationFilter::SeqNameFilter("chr20"),
                                                       AnnotationFilter::SeqNameFilter("chr21"),AnnotationFilter::SeqNameFilter("chr22"),AnnotationFilter::SeqNameFilter("chrX"),AnnotationFilter::SeqNameFilter("chrY"),
                                                       AnnotationFilter::SeqNameFilter("chrM"),logicOp = rep(c("|"),24)),.GlobalEnv)
  
  genesEns <- suppressWarnings(as.data.table(ensembldb::genes(EnsDb,return.type='DataFrame')))
  genesEns <- genesEns[genesEns$seq_name%in%paste0('chr',c(1:22,'X','Y','M')),'symbol']
  colnames(genesEns) <- c("gene_symbol")
  
  updateSelectizeInput(session, "genes",
                       choices = as.character(genesEns$gene_symbol),
                       server = TRUE,selected = NULL)
  assign('genesEns',genesEns,.GlobalEnv)
  
  assign('allEns',ensembldb::genes(EnsDb,columns="gene_biotype",return.type="DataFrame")$gene_biotype,.GlobalEnv)
  assign('allEns_tx',ensembldb::genes(EnsDb,columns="tx_biotype",return.type="DataFrame")$tx_biotype,.GlobalEnv)
  assign('levelsColor',unique(allEns),.GlobalEnv)
  assign('levelsColor_tx',unique(allEns_tx),.GlobalEnv)
  #mapColor <- c(grDevices::rainbow(10),grDevices::rainbow(10,start = 0.2),grDevices::rainbow(length(levelsColor)-20,start = 0.4))
  mapColor <- grDevices::rainbow(length(levelsColor))
  c25 <- c(
    "#1c86ee", "#E31A1C", "#008b00", "#6A3D9A","#FF7F00",
    "#d4d4d4", "#ffd700", "#7ec0ee", "#FB9A99", "#90ee90",
    "#CAB2D6","#FDBF6F", "#b3b3b3", "#eee685", "#b03060",
    "#ff83fa", "#ff1493", "#0000ff", "#36648b", "#00ced1",
    "#00ff00", "#8b8b00", "#cdcd00", "#8b4500", "#a52a2a"
  )
  
  assign('mapColor',c(c25,c25)[1:length(levelsColor)],.GlobalEnv)
  assign('mapColor_tx',c(c25,c25)[1:length(levelsColor_tx)],.GlobalEnv)
  assign('tabEns',table(allEns)[match(levelsColor,names(table(allEns)))],.GlobalEnv)
  assign('tabEns_tx',table(allEns_tx)[match(levelsColor_tx,names(table(allEns_tx)))],.GlobalEnv)
  assign('legend_idx',order(tabEns,decreasing = T),.GlobalEnv)
  assign('legend_idx_tx',order(tabEns_tx,decreasing = T),.GlobalEnv)
  


  if(mouse){
    square_imm = paste0("mus.svg")
    corner_imm = paste0("mus",assembly,"_corner.png")
  }else{
    square_imm = paste0("david.svg")
    corner_imm = paste0("david",assembly,"_corner.png")
  }
  insertUI(
    selector = ".logo-lg",
    where = "afterEnd",
    ui = tagList(div(img(class = "logo-lg", 
                         id='logo-lg', 
                         src=corner_imm,
                         width = "100%", 
                         height = "auto",
                         style="float: left;display:inline-block;margin-left:-25px;")),
                 img(class = "logo-sq", 
                     id='logo-sq', 
                     src = square_imm)
    )
  )
  # removeModal()
  updateTabItems(session, "sideBar", "home")
  
  shinyjs::show("ok4")
  if(simpleDebug){cat(file=stderr(), paste0('End LoadAssembly\n'))}
}

if(simpleDebug){cat(file=stderr(), paste0('step4 globalFile\n'))}
