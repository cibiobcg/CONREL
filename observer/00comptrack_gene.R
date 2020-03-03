# Observer that generates the Composite track by gene selection
observeEvent(input$searchGenes, ignoreInit = T,{
  source(file.path("scripts", "tab.R"),  local = TRUE)$value
  clickGene <<- TRUE
  if(simpleDebug){print("View genome tab")}
  source(file.path("scripts", "genomeMenu.R"),  local = TRUE)$value
  updateTabItems(session, "sideBar", "genome")
  withProgress(message = 'Generating plot', value = 0, {
    if(simpleDebug){print("Composite track start")}
    validate(
      validateGenes(input),
      validateInput(input)
    )
    incProgress(0.1,detail = "Loading genes/trascripts track")
    
    ###########################
    if(simpleDebug){print("Search by gene")}

    f = AnnotationFilter(~ (symbol == input$genes & chrFilter))
    gr_tmp <- biovizBase::crunch(EnsDb.Hsapiens.v75, f)
    
    elements <<- c(seqlevels(gr_tmp),min(start(gr_tmp)),max(end(gr_tmp)))
    window_load = window2Load(elements)
    
    if(input$choice_track=='gene') {
      if(simpleDebug){print("Search by gene - gene view")}
      # GENE track
      gene <- genes(EnsDb.Hsapiens.v75,
                    filter=AnnotationFilterList(GeneStartFilter(min(start(gr_tmp))-1000000,condition = '>'),
                                                GeneEndFilter(max(end(gr_tmp))+1000000,condition = '<'),
                                                SeqNameFilter(seqlevels(gr_tmp)),
                                                logicOp = c('&','&')))
      tx <- TnT::FeatureTrack(gene, tooltip = as.data.frame(gene),
                              names = paste(gene$symbol, " (", gene$gene_biotype, ")", sep = ""),
                              color = TnT::mapcol(gene$gene_biotype, palette.fun = grDevices::rainbow),
                              background = "#eeeeee",label = NULL)
    } else {
      if(simpleDebug){print("Search by gene - transcript view")}
      # TRANSCRIPT track
      gr <- biovizBase::crunch(EnsDb.Hsapiens.v75,
                               GRanges(seqlevels(gr_tmp),IRanges(min(start(gr_tmp))-1000000,
                                                                 max(end(gr_tmp))+1000000)))
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
    
    elements <<- c(paste0(seqlevelsInUse(gr_tmp)),min(start(gr_tmp)),max(end(gr_tmp)))
    
    # SNPs track
    if(length(input$snps)>0){
      incProgress(0.1,detail = "Loading SNPs track")
      track_snps = loadSNPs(elements,window_load)
      tx = append(tx,track_snps)
    }
    #########################################
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