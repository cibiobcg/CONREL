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
        gr_tmp <- biovizBase::crunch(EnsDb, f)
        
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
        gene <- genes(EnsDb,
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
        
        gr <- biovizBase::crunch(EnsDb,
                                 GRanges(chr_transcriptView,IRanges(start_transcriptView,end_transcriptView)))
        tx <- TnT::TxTrackFromGRanges(gr, color = "grey2",label = NULL,background = "#eeeeee")
        trackData(tx)$tooltip <- select(EnsDb,
                                        keys = tx$tooltip$tx_id,
                                        keytype = "TXID",
                                        columns = c("GENEID", "SYMBOL", "TXBIOTYPE"))
        trackData(tx)$color <- mapColor_tx[match(tx$tooltip$TXBIOTYPE,levelsColor_tx)]
        # trackData(tx)$color <- TnT::mapcol(tx$tooltip$TXBIOTYPE)
        trackData(tx)$display_label <- TnT::strandlabel(
          paste(tx$tooltip$SYMBOL, tx$tooltip$TXBIOTYPE),
          strand(TnT::trackData(tx)))
      }
      
      if(clickGene){
        elements <<- c(paste0(seqlevelsInUse(gr_tmp)),min(start(gr_tmp)),max(end(gr_tmp)))  
      }
      
      # TSS track
      if(length(input$tss)>0){
        incProgress(0.1,detail = "Loading TSS track")
        track_tss = loadTSS(elements,window_load)
        tx = append(tx,track_tss)
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