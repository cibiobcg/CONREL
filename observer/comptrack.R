# Observer that generates the Composite track by position
observeEvent(c(input$searchSetting,input$searchSetting_bottom), ignoreInit = T,{

  output$peaksError <- NULL
  output$regElementError <- NULL
  output$consensusError <- NULL
  output$tissueError <- NULL
  output$tbaError <- NULL
  output$cellLineError <- NULL
  
  if(simpleDebug){print(paste0("Search by position: ",clickPos))}
  if(simpleDebug){print(paste0("Search by gene: ",clickGene))}
  # source(file.path("scripts", "tab.R"),  local = TRUE)$value
  # Check the input on the search2 tab and print the error
  if(simpleDebug){print("Check Conditions")}
  peak <- length(input$peaks)==0
  regElement <- length(input$regElement)==0
  consensus <- length(input$consensus)==0
  tissue <- "tissue" %in% input$consensus & length(input$tissue)==0
  cellLine <- "cell-line" %in% input$consensus & (is.null(input$tree) | length(get_selected(input$tree, format = "slices"))==0)
  tba <- length(input$tba)>0 & length(input$pvalues)==0
  
  if(peak){
    output$peaksError <- renderUI({
      p("Please select at least one type of peaks")
    })
  }
  if(regElement) {
    output$regElementError <- renderUI({
      p("Please select at least one regulatory element")
    })
  }
  if(consensus){
    output$consensusError <- renderUI({
      p("Please select at least one type of consensus regions")
    })
  }
  if(tissue){
    output$tissueError <- renderUI({
      p("Please select at least one tissue, or deselect 'Tissue consensus regions'")
    })
  }
  if(cellLine){
    output$cellLineError <- renderUI({
      p("Please select at least one cell line, or deselect 'Cell-line consensus regions'")
    })
  }
  if(tba){
    output$tbaError <- renderUI({
      p("Please select at least one pvalue threshold")
    })
  }
  
  if(!any(peak,regElement,consensus,tissue,tba,cellLine)) {
    if(simpleDebug){print("View genome tab")}
    source(file.path("scripts", "genomeMenu.R"),  local = TRUE)$value
    updateTabItems(session, "sideBar", "genome")
    withProgress(message = 'Generating plot', value = 0, {
      if(simpleDebug){print("Composite track start")}
      incProgress(0.1,detail = "Loading genes/trascripts track")
      
      ##############################################################
      if(clickPos){
        region = input$region
        
        elements <<- convertPosition(region,TRUE)
      }
      if(clickGene){
        f = AnnotationFilter(~ (symbol == input$genes & chrFilter))
        gr_tmp <- biovizBase::crunch(EnsDb.Hsapiens.v75, f)
        
        elements <<- c(seqlevels(gr_tmp),min(start(gr_tmp)),max(end(gr_tmp)))
      }
      
      window_load = window2Load(elements)
      
      if(input$choice_track=='gene') {
        if(simpleDebug){print("Load gene track view")}
        # GENE track
        if(clickPos){
          start_geneView = window_load[1]
          end_geneView = window_load[2]
          chr_geneView = elements[1]
        }
        if(clickGene){
          start_geneView = min(start(gr_tmp))-increaseWindow
          end_geneView = max(end(gr_tmp))+increaseWindow
          chr_geneView = seqlevels(gr_tmp)
        }
        gene <- genes(EnsDb.Hsapiens.v75,
                      filter=AnnotationFilterList(GeneStartFilter(start_geneView,condition = '>'),
                                                  GeneEndFilter(end_geneView,condition = '<'),
                                                  SeqNameFilter(chr_geneView),
                                                  logicOp = c('&','&')))
        tx <- TnT::FeatureTrack(gene, tooltip = as.data.frame(gene),
                                names = paste(gene$symbol, " (", gene$gene_biotype, ")", sep = ""),
                                color = mapColor[match(gene$gene_biotype,levelsColor)],
                                #color = TnT::mapcol(gene$gene_biotype, palette.fun = grDevices::rainbow),
                                background = "#eeeeee",label = NULL)
      } else {
        if(simpleDebug){print("Load transcript track view")}
        # TRANSCRIPT track
        if(clickPos){
          start_transcriptView = window_load[1]
          end_transcriptView = window_load[2]
          chr_transcriptView = elements[1]
        }
        if(clickGene){
          start_transcriptView = min(start(gr_tmp))-increaseWindow
          end_transcriptView = max(end(gr_tmp))+increaseWindow
          chr_transcriptView = seqlevels(gr_tmp)
        }
        
        gr <- biovizBase::crunch(EnsDb.Hsapiens.v75,
                                 GRanges(chr_transcriptView,IRanges(start_transcriptView,end_transcriptView)))
        tx <- TnT::TxTrackFromGRanges(gr, color = "grey2",label = NULL,background = "#eeeeee")
        trackData(tx)$tooltip <- select(EnsDb.Hsapiens.v75,
                                        keys = tx$tooltip$tx_id,
                                        keytype = "TXID",
                                        columns = c("GENEID", "SYMBOL", "TXBIOTYPE"))
        trackData(tx)$color <- TnT::mapcol(tx$tooltip$TXBIOTYPE)
        trackData(tx)$display_label <- TnT::strandlabel(
          paste(tx$tooltip$SYMBOL, tx$tooltip$TXBIOTYPE),
          strand(TnT::trackData(tx)))
      }
      
      if(clickGene){
        elements <<- c(paste0(seqlevelsInUse(gr_tmp)),min(start(gr_tmp)),max(end(gr_tmp)))  
      }
      
      # SNPs track
      if(length(input$snps)>0){
        incProgress(0.1,detail = "Loading SNPs track")
        track_snps = loadSNPs(elements,window_load)
        tx = append(tx,track_snps)
      }
      
      ##############################################
      emptyTrack = TnT::BlockTrack(GRanges(elements[1],IRanges(0,1)),
                                   color = "#EEEEEE",background = "#EEEEEE",
                                   height = 12,label=NULL)
      tx = append(emptyTrack,tx)
      incStep = 0.8/length(input$consensus)
      
      if('global' %in% input$consensus) {
        if(simpleDebug){print("Load global consensus")}
        source(file.path("scripts", "globalConsensus.R"),  local = TRUE)$value
        
        tx = append(tx,tx_gCons)
      }
      if('tissue' %in% input$consensus) {
        if(simpleDebug){print("Load tissue consensus")}
        source(file.path("scripts", "tissueConsensus.R"),  local = TRUE)$value
        
        tx = append(tx,tx_tCons)
      }
      if('cell-line' %in% input$consensus) {
        if(simpleDebug){print("Load cell-line consensus")}
        source(file.path("scripts", "clineConsensus.R"),  local = TRUE)$value
        
        tx = append(tx,tx_cCons)
      }
      
      
      trackRendered = TnTBoard(tx,
                               view.range = GRanges(elements[1],IRanges(as.numeric(elements[2]),as.numeric(elements[3]))),
                               use.tnt.genome = T,
                               coord.range = IRanges(window_load[1],window_load[2]))
      track(trackRendered)
    })
    clickPos <<- FALSE
    clickGene <<- FALSE
    removeUI(
      ## pass in appropriate div id
      selector = '#tbaTab'
    )
    removeUI(
      ## pass in appropriate div id
      selector = '#regionTab'
    )
  }
})



observeEvent(input$example1, ignoreInit = T,{
  
  output$peaksError <- NULL
  output$regElementError <- NULL
  output$consensusError <- NULL
  output$tissueError <- NULL
  output$tbaError <- NULL
  output$cellLineError <- NULL
  
  clickPos <<- TRUE
  region = c("chrX:66,764,465-66,950,461")
  
  if(simpleDebug){print(paste0("Search by position: ",clickPos))}
  if(simpleDebug){print(paste0("Search by gene: ",clickGene))}
  
  
  if(simpleDebug){print("View genome tab")}
  source(file.path("scripts", "genomeMenu.R"),  local = TRUE)$value
  updateTabItems(session, "sideBar", "genome")
  withProgress(message = 'Generating plot', value = 0, {
    if(simpleDebug){print("Composite track start")}
    incProgress(0.1,detail = "Loading genes/trascripts track")
    
    ##############################################################
    
    elements <<- convertPosition(region,TRUE)
    
    
    window_load = window2Load(elements)
    
    
    if(simpleDebug){print("Load gene track view")}
    # GENE track
    if(clickPos){
      start_geneView = window_load[1]
      end_geneView = window_load[2]
      chr_geneView = elements[1]
    }
    
    gene <- genes(EnsDb.Hsapiens.v75,
                  filter=AnnotationFilterList(GeneStartFilter(start_geneView,condition = '>'),
                                              GeneEndFilter(end_geneView,condition = '<'),
                                              SeqNameFilter(chr_geneView),
                                              logicOp = c('&','&')))
    tx <- TnT::FeatureTrack(gene, tooltip = as.data.frame(gene),
                            names = paste(gene$symbol, " (", gene$gene_biotype, ")", sep = ""),
                            color = mapColor[match(gene$gene_biotype,levelsColor)],
                            #color = TnT::mapcol(gene$gene_biotype, palette.fun = grDevices::rainbow),
                            background = "#eeeeee",label = NULL)
    
    
    
    
    # SNPs track
    incProgress(0.1,detail = "Loading SNPs track")
    track_snps = loadSNPs(elements,window_load)
    tx = append(tx,track_snps)
    
    ##############################################
    emptyTrack = TnT::BlockTrack(GRanges(elements[1],IRanges(0,1)),
                                 color = "#EEEEEE",background = "#EEEEEE",
                                 height = 12,label=NULL)
    tx = append(emptyTrack,tx)
    incStep = 0.8/1
    
    
    if(simpleDebug){print("Load global consensus")}
    
    loadDF = mappingInput("promoter","narrow")
    tx_gCons = lapply(1:nrow(loadDF),function(i){
      incProgress(incStep/nrow(loadDF),detail = "Loading global consensus regions")
      # gCons = fread(paste0(encodeFolder,loadDF[i,1],"Peaks/globalConsensus/",loadDF[i,2],".union.bed"),data.table = F)
      # gCons = gCons[gCons$V1==elements[1],]
      gCons = system(paste0("tabix ",inputFolder,"consensus/",loadDF[i,1],"Peaks/bgz_globalConsensus/",loadDF[i,2],".union.bed.gz ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
      gCons = data.frame(do.call(rbind, strsplit(gCons, "\t", fixed=TRUE)))
      if(nrow(gCons)==0){
        generateEmptyTrack(elements[1])
      } else {
        GRanges_gcons = makeGRangesFromDataFrame(gCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = T)
        FeatureTrack(
          GRanges_gcons,
          tooltip = as.data.frame(GRanges_gcons),
          label = paste0("global consensus tracks - ",gsub("active","active enhancer",loadDF[i,2])," , ",loadDF[i,1]," peaks"),
          names = "",height = 25,color = consensusColor(loadDF[i,2],loadDF[i,1]),
          background = "#eeeeee")
      }
    })
    tx = append(tx,tx_gCons)
    
    
    trackRendered = TnTBoard(tx,
                             view.range = GRanges(elements[1],IRanges(as.numeric(elements[2]),as.numeric(elements[3]))),
                             use.tnt.genome = T,
                             coord.range = IRanges(window_load[1],window_load[2]))
    track(trackRendered)
  })
  clickPos <<- FALSE
  clickGene <<- FALSE
  removeUI(
    ## pass in appropriate div id
    selector = '#tbaTab'
  )
  removeUI(
    ## pass in appropriate div id
    selector = '#regionTab'
  )
})


observeEvent(input$example2, ignoreInit = T,{
  
  output$peaksError <- NULL
  output$regElementError <- NULL
  output$consensusError <- NULL
  output$tissueError <- NULL
  output$tbaError <- NULL
  output$cellLineError <- NULL
  
  clickPos <<- TRUE
  region = c("chr19:51,353,417-51,367,543")
  
  if(simpleDebug){print(paste0("Search by position: ",clickPos))}
  if(simpleDebug){print(paste0("Search by gene: ",clickGene))}
  
  
  if(simpleDebug){print("View genome tab")}
  source(file.path("scripts", "genomeMenu.R"),  local = TRUE)$value
  updateTabItems(session, "sideBar", "genome")
  withProgress(message = 'Generating plot', value = 0, {
    if(simpleDebug){print("Composite track start")}
    incProgress(0.1,detail = "Loading genes/trascripts track")
    
    ##############################################################
    
    elements <<- convertPosition(region,TRUE)
    
    
    window_load = window2Load(elements)
    
    
    if(simpleDebug){print("Load gene track view")}
    # GENE track
    if(clickPos){
      start_geneView = window_load[1]
      end_geneView = window_load[2]
      chr_geneView = elements[1]
    }
    
    gene <- genes(EnsDb.Hsapiens.v75,
                  filter=AnnotationFilterList(GeneStartFilter(start_geneView,condition = '>'),
                                              GeneEndFilter(end_geneView,condition = '<'),
                                              SeqNameFilter(chr_geneView),
                                              logicOp = c('&','&')))
    tx <- TnT::FeatureTrack(gene, tooltip = as.data.frame(gene),
                            names = paste(gene$symbol, " (", gene$gene_biotype, ")", sep = ""),
                            color = mapColor[match(gene$gene_biotype,levelsColor)],
                            #color = TnT::mapcol(gene$gene_biotype, palette.fun = grDevices::rainbow),
                            background = "#eeeeee",label = NULL)
    
    ##############################################
    emptyTrack = TnT::BlockTrack(GRanges(elements[1],IRanges(0,1)),
                                 color = "#EEEEEE",background = "#EEEEEE",
                                 height = 12,label=NULL)
    tx = append(emptyTrack,tx)
    incStep = 0.8/1
    
    
    if(simpleDebug){print("Load global consensus")}
    
    loadDF = mappingInput(c("promoter","enhancer","active"),c("narrow","broad"))
    tx_gCons = lapply(1:nrow(loadDF),function(i){
      incProgress(incStep/nrow(loadDF),detail = "Loading global consensus regions")
      # gCons = fread(paste0(encodeFolder,loadDF[i,1],"Peaks/globalConsensus/",loadDF[i,2],".union.bed"),data.table = F)
      # gCons = gCons[gCons$V1==elements[1],]
      gCons = system(paste0("tabix ",inputFolder,"consensus/",loadDF[i,1],"Peaks/bgz_globalConsensus/",loadDF[i,2],".union.bed.gz ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
      gCons = data.frame(do.call(rbind, strsplit(gCons, "\t", fixed=TRUE)))
      if(nrow(gCons)==0){
        generateEmptyTrack(elements[1])
      } else {
        GRanges_gcons = makeGRangesFromDataFrame(gCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = T)
        FeatureTrack(
          GRanges_gcons,
          tooltip = as.data.frame(GRanges_gcons),
          label = paste0("global consensus tracks - ",gsub("active","active enhancer",loadDF[i,2])," , ",loadDF[i,1]," peaks"),
          names = "",height = 25,color = consensusColor(loadDF[i,2],loadDF[i,1]),
          background = "#eeeeee")
      }
    })
    tx = append(tx,tx_gCons)
    
    loadDF = tissueAvailability("prostate",c("narrow","broad"),c("promoter","enhancer","active"))
    tx_tCons = lapply(1:nrow(loadDF),function(i){
      splitFileName = strsplit(as.character(loadDF[i,2]),"\\.")[[1]]
      incProgress(incStep/nrow(loadDF),detail = "Loading tissue consensus regions")
      # tCons = fread(paste0(encodeFolder,loadDF[i,1],"Peaks/tissueConsensus/",as.character(loadDF[i,2])),data.table = F)
      # tCons = tCons[tCons$V1==elements[1],]
      tCons = system(paste0("tabix ",inputFolder,"consensus/",loadDF[i,1],"Peaks/bgz_tissueConsensus/",as.character(loadDF[i,2])," ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
      tCons = data.frame(do.call(rbind, strsplit(tCons, "\t", fixed=TRUE)))
      if(nrow(tCons)==0){
        generateEmptyTrack(elements[1])
      } else {
        GRanges_tcons = makeGRangesFromDataFrame(tCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = TRUE)
        FeatureTrack(
          GRanges_tcons,
          tooltip = as.data.frame(GRanges_tcons),
          label = paste0("tissue consensus tracks - ",splitFileName[1]," - ",splitFileName[2]," , ",loadDF[i,1]," peaks"),
          names = "",height = 25,color = consensusColor(loadDF[i,2],loadDF[i,1]),
          background = "#eeeeee")
      }
    })
    tx = append(tx,tx_tCons)

    
    trackRendered = TnTBoard(tx,
                             view.range = GRanges(elements[1],IRanges(as.numeric(elements[2]),as.numeric(elements[3]))),
                             use.tnt.genome = T,
                             coord.range = IRanges(window_load[1],window_load[2]))
    track(trackRendered)
    
    
    
    
  })
  clickPos <<- FALSE
  clickGene <<- FALSE
  removeUI(
    ## pass in appropriate div id
    selector = '#tbaTab'
  )
  removeUI(
    ## pass in appropriate div id
    selector = '#regionTab'
  )
})