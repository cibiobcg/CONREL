# a function to check the presence of the folder with data, otherwise the buttons are disabled
checkFolderAssembly <- function() {
  if(!dir.exists('/CONREL_Data/hg19')){disable('hg19')}
  if(!dir.exists('/CONREL_Data/hg38')){disable('hg38')}
  if(!dir.exists('/CONREL_Data/mm10')){disable('mm10')}
}

# Return the UI for a modal dialog with data selection input. If 'failed' is
# TRUE, then display a message that the previous value was invalid.
dataModal <- function(failed = FALSE) {
  modalDialog(title = "",
              div(
                id = "main_content",
                fluidRow(
                  column(2),
                  column(10,align = "center",
                         p('Species',style = 'font-size: 24px;color: black;font-family: "Georgia", Times, "Times New Roman", serif;font-weight: bold;'))
                ),
                fluidRow(
                  column(2),
                  column(5,align = "center",p('Human',style = 'font-size: 21px;color: dimgrey;font-family: "Georgia", Times, "Times New Roman", serif;font-weight: bold;')),
                  column(5,align = "center",p('Mouse',style = 'font-size: 21px;color: dimgrey;font-family: "Georgia", Times, "Times New Roman", serif;font-weight: bold;'))
                ),
                fluidRow(
                  column(2,align = "center",style = "height:300px;",
                         p('Genome Assembly',style = "display:flex;height:300px;
                                                    font-size: 24px;color: black;font-family: 'Georgia', Times, 'Times New Roman', serif;font-weight: bold;
                                                    justify-content:center;
                                                    align-items:center;
                                                    writing-mode: vertical-lr;
                                                    transform: rotate(-180deg)")),
                  column(10,style = "height:300px;",
                         fluidRow(
                           column(6,align = "center",
                                  actionButton('hg38','',width = '100%',
                                               style = "color:red;font-weight: bold;width: 120px; height: 120px;background: url('human38.png');  background-size: cover; background-position: center;")),
                           column(6,align = "center",
                                  actionButton('mm10','',width='100%',
                                               style = "color:red;font-weight: bold;width: 120px; height: 120px;background: url('mouse.png');  background-size: cover; background-position: center;"))),
                         fluidRow(
                           column(6,align = "center",
                                  actionButton('hg19','',width='100%',
                                               style = "color:red;font-weight: bold;width: 120px; height: 120px;background: url('human19.png');  background-size: cover; background-position: center;")),
                           column(6)))
                )),
                div(
                  id = "loading_page",
                  fluidRow(align='center',
                           span('Please wait, loading data....'),hr(),
                           spin_2(),
                  )
              ),
              if (failed)
                div(tags$b("An error has occured", style = "color: red;")),
              
              footer = tagList(
              )
  )
}
# Show modal when button is clicked.
observeEvent(input$sideBar, ignoreInit = T, {
  if(input$sideBar=="changeORG") {
    showModal(dataModal())
    checkFolderAssembly()
    shinyjs::hide("loading_page")
  }
})

# a function in common with each of the observer
loadAssemblyData <- function(assembly) {
  if(simpleDebug){print(paste0("Load...",assembly))}

  removeUI(
    ## pass in appropriate div id
    selector = '#logo-lg'
  )
  removeUI(
    ## pass in appropriate div id
    selector = '#logo-sq'
  )
  
  base_path <<- paste0("/CONREL_Data/",assembly,"/")
  encodeFolder <<- paste0("/CONREL_Data/",assembly,"/consensus/")
  inputFolder <<- paste0("/CONREL_Data/",assembly,"/")
  TBA_folder <<- paste0("/CONREL_Data/",assembly,"/")
  hg19 <<- fread(paste0(inputFolder,"data/",assembly,".chrom.bed"),data.table = F)
  hg_assembly <<- assembly
  
  if(assembly=='hg19') {
    EnsDb <<- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
  } else if(assembly=='hg38') {
    EnsDb <<- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  } else if(assembly=='mm10') {
    EnsDb <<- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
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
    return(TnT::BlockTrack(GRanges(chr,IRanges(0,1)),
                           color = "#EEEEEE",background = "#EEEEEE",
                           height = 15,label=NULL))
  }
  
  
  
  
  
  checkValidRegion <<- function(region) {
    elements <<- convertPosition(region,TRUE)
    window_load <<- window2Load(elements)
    start_geneView <<- window_load[1]
    end_geneView <<- window_load[2]
    chr_geneView <<- elements[1]
    gene <<- genes(EnsDb,
                   filter=AnnotationFilterList(GeneStartFilter(start_geneView,condition = '>'),
                                               GeneEndFilter(end_geneView,condition = '<'),
                                               SeqNameFilter(chr_geneView),
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
  
  
  ### COUNT PFM ###
  load(paste0(inputFolder,"data/countPFMS.RData"))
  minCount <<- 50
  
  
  loadSNPs <<- function(elements,window_load){
    snps.filtered = system(paste0("tabix ",inputFolder,"dbSNP_v151_hg19/dbSNP_151_hg19_complete_TOPMED_ALL.bed.gz ",gsub("chr","",elements[1]),":",window_load[1],"-",window_load[2]),intern = T)
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
  
  
  
  loadTSS <<- function(elements,window_load){
    tss.filtered = system(paste0("tabix ",inputFolder,"data/UCSC_TSS.bed.gz ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
    tss.filtered = data.frame(do.call(rbind, strsplit(tss.filtered, "\t", fixed=TRUE)))
    
    ir <- IRanges(start = as.numeric(as.character(tss.filtered$X2)), end = as.numeric(as.character(tss.filtered$X3)))
    gpos <- GRanges(elements[1], ir)
    # gpos$score <- rep(1,nrow(tss.filtered))
    # snps.filtered$alleles = paste(snps.filtered$X4,snps.filtered$X5,sep=">")
    gtooltip <- tss.filtered[,c(4,5)]
    colnames(gtooltip) = c("strand","confScore")
    pt = BlockTrack(gpos,color = "black",label = "SwitchGear TSS", tooltip = gtooltip,background = "#eeeeee",
                    height = 15)
    emptyTrack = TnT::BlockTrack(GRanges(elements[1],IRanges(0,1)),
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
  
  seqlevelsStyle(EnsDb) <<- "UCSC"
  genesEns <<- suppressWarnings(as.data.table(genes(EnsDb)$symbol))
  allEns <<- genes(EnsDb,columns="gene_biotype",return.type="DataFrame")$gene_biotype
  allEns_tx <<- genes(EnsDb,columns="tx_biotype",return.type="DataFrame")$tx_biotype
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
  chrFilter <<- AnnotationFilterList(SeqNameFilter("chr1"),SeqNameFilter("chr2"),SeqNameFilter("chr3"),SeqNameFilter("chr4"),
                                     SeqNameFilter("chr5"),SeqNameFilter("chr6"),SeqNameFilter("chr7"),SeqNameFilter("chr8"),
                                     SeqNameFilter("chr9"),SeqNameFilter("chr10"),SeqNameFilter("chr11"),SeqNameFilter("chr12"),
                                     SeqNameFilter("chr13"),SeqNameFilter("chr14"),SeqNameFilter("chr15"),SeqNameFilter("chr16"),
                                     SeqNameFilter("chr17"),SeqNameFilter("chr18"),SeqNameFilter("chr19"),SeqNameFilter("chr20"),
                                     SeqNameFilter("chr21"),SeqNameFilter("chr22"),SeqNameFilter("chrX"),SeqNameFilter("chrY"),
                                     SeqNameFilter("chrM"),logicOp = rep(c("|"),24))
  
  updateSelectizeInput(session, "genes",
                       choices = as.character(genesEns$gene_symbol),
                       server = TRUE,selected = NULL)
  insertUI(
    selector = ".logo-lg",
    where = "afterEnd",
    ui = tagList(div(img(class = "logo-lg", 
                         id='logo-lg', 
                         src=paste0("human",assembly,"_corner.png"),
                         width = "100%", 
                         height = "auto",
                         style="float: left;display:inline-block;margin-left:-25px;")),
                 img(class = "logo-sq", 
                     id='logo-sq', 
                     src = paste0("human",assembly,"_square.svg"))
    )
  )
  # removeModal()
  updateTabItems(session, "sideBar", "home")
}


observeEvent(input$hg19,ignoreInit = T, {
  # Check that data object exists and is data frame.
  shinyjs::hide("main_content")
  shinyjs::show("loading_page")
  
  tagWorkInProgress = tagList(NULL)
  
  loadAssemblyData('hg19')
  
  removeModal()
  
})
observeEvent(input$hg38,ignoreInit = T, {
  # Check that data object exists and is data frame.
  shinyjs::hide("main_content")
  shinyjs::show("loading_page")
  
  tagWorkInProgress <<- tagList(div(style="border:2px red dotted;",
                                    p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
                                      "ATTENTION! WORK IN PROGESS"),
                                    p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
                                      "The GRCh38 version of CONREL is online but not fully functional. Some TBA annotations are still on processing, then the TBA informations will not show for those CREs.")))
  tagWorkInProgress <<- tagList(NULL)
  loadAssemblyData('hg38')
  removeModal()
  
})
observeEvent(input$mm10,ignoreInit = T, {
  # Check that data object exists and is data frame.
  shinyjs::hide("main_content")
  shinyjs::show("loading_page")
  
  tagWorkInProgress <<- tagList(div(style="border:2px red dotted;",
                                    p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
                                      "ATTENTION! WORK IN PROGESS"),
                                    p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
                                      "The GRCh38 version of CONREL is online but not fully functional. Some TBA annotations are still on processing, then the TBA informations will not show for those CREs.")))
  tagWorkInProgress <<- tagList(NULL)
  loadAssemblyData('mm10')
  removeModal()
  
})
