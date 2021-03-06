# Observer that generates the Composite track by position
observeEvent(input$Run, ignoreInit = T,{
  source(file.path("scripts", "tab.R"),  local = TRUE)$value
  clickPos <<- TRUE
  if(simpleDebug){print("View genome tab")}
  source(file.path("scripts", "genomeMenu.R"),  local = TRUE)$value
  updateTabItems(session, "sideBar", "genome")
  withProgress(message = 'Generating plot', value = 0, {
    if(simpleDebug){print("Composite track start")}
    validate(
      validateRegion(input),
      validateInput(input)
    )
    incProgress(0.1,detail = "Loading genes/trascripts track")
    
    ##############################################################
    if(simpleDebug){print("Search by position")}

    region = input$region
    
    elements <<- convertPosition(region,TRUE)
    window_load = window2Load(elements)
    
    if(input$choice_track=='gene') {
      if(simpleDebug){print("Search by position - gene view")}
      # GENE track
      gene <- genes(EnsDb.Hsapiens.v75,
                    filter=AnnotationFilterList(GeneStartFilter(window_load[1],condition = '>'),
                                                GeneEndFilter(window_load[2],condition = '<'),
                                                SeqNameFilter(elements[1]),
                                                logicOp = c('&','&')))
      tx <- TnT::FeatureTrack(gene, tooltip = as.data.frame(gene),
                              names = paste(gene$symbol, " (", gene$gene_biotype, ")", sep = ""),
                              color = TnT::mapcol(gene$gene_biotype, palette.fun = grDevices::rainbow),
                              background = "#eeeeee",label = NULL)
    } else {
      if(simpleDebug){print("Search by position - transcript view")}
      gr <- biovizBase::crunch(EnsDb.Hsapiens.v75,
                               GRanges(elements[1],IRanges(window_load[1],window_load[2])))
      
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
    
    # SNPs track
    if(length(input$snps)>0){
      incProgress(0.1,detail = "Loading SNPs track")
      track_snps = loadSNPs(elements,window_load)
      tx = append(tx,track_snps)
    }
    ##############################################
    emptyTrack = TnT::BlockTrack(GRanges(elements[1],IRanges(0,1)),
                                 color = "#EEEEEE",background = "#EEEEEE",
                                 height = 15,label=NULL)
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
  
})