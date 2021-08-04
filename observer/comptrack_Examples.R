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
    
    elements <<- convertPosition(regionEx1,TRUE)
    
    
    window_load = window2Load(elements)
    
    
    if(simpleDebug){print("Load gene track view")}
    # GENE track
    if(clickPos){
      start_geneView = window_load[1]
      end_geneView = window_load[2]
      chr_geneView = elements[1]
    }
    
    gene <- ensembldb::genes(EnsDb,
                             filter=AnnotationFilter::AnnotationFilterList(AnnotationFilter::GeneStartFilter(start_geneView,condition = '>'),
                                                                           AnnotationFilter::GeneEndFilter(end_geneView,condition = '<'),
                                                                           AnnotationFilter::SeqNameFilter(chr_geneView),
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
    emptyTrack = TnT::BlockTrack(GenomicRanges::GRanges(elements[1],IRanges::IRanges(0,1)),
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
        GRanges_gcons = GenomicRanges::makeGRangesFromDataFrame(gCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = T)
        TnT::FeatureTrack(
          GRanges_gcons,
          tooltip = as.data.frame(GRanges_gcons),
          label = paste0("global consensus tracks - ",gsub("active","active enhancer",loadDF[i,2])," , ",loadDF[i,1]," peaks"),
          names = "",height = 25,color = consensusColor(loadDF[i,2],loadDF[i,1]),
          background = "#eeeeee")
      }
    })
    tx = append(tx,tx_gCons)
    
    
    trackRendered = TnT::TnTBoard(tx,
                             view.range = GenomicRanges::GRanges(elements[1],IRanges::IRanges(as.numeric(elements[2]),as.numeric(elements[3]))),
                             use.tnt.genome = T,
                             coord.range = IRanges::IRanges(window_load[1],window_load[2]))
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
    
    elements <<- convertPosition(regionEx2,TRUE)
    
    
    window_load = window2Load(elements)
    
    
    if(simpleDebug){print("Load gene track view")}
    # GENE track
    if(clickPos){
      start_geneView = window_load[1]
      end_geneView = window_load[2]
      chr_geneView = elements[1]
    }
    
    gene <- ensembldb::genes(EnsDb,
                             filter=AnnotationFilter::AnnotationFilterList(AnnotationFilter::GeneStartFilter(start_geneView,condition = '>'),
                                                                           AnnotationFilter::GeneEndFilter(end_geneView,condition = '<'),
                                                                           AnnotationFilter::SeqNameFilter(chr_geneView),
                                                                           logicOp = c('&','&')))
    tx <- TnT::FeatureTrack(gene, tooltip = as.data.frame(gene),
                            names = paste(gene$symbol, " (", gene$gene_biotype, ")", sep = ""),
                            color = mapColor[match(gene$gene_biotype,levelsColor)],
                            #color = TnT::mapcol(gene$gene_biotype, palette.fun = grDevices::rainbow),
                            background = "#eeeeee",label = NULL)
    
    ##############################################
    emptyTrack = TnT::BlockTrack(GenomicRanges::GRanges(elements[1],IRanges::IRanges(0,1)),
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
        GRanges_gcons = GenomicRanges::makeGRangesFromDataFrame(gCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = T)
        TnT::FeatureTrack(
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
        GRanges_tcons = GenomicRanges::makeGRangesFromDataFrame(tCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = TRUE)
        TnT::FeatureTrack(
          GRanges_tcons,
          tooltip = as.data.frame(GRanges_tcons),
          label = paste0("tissue consensus tracks - ",splitFileName[1]," - ",splitFileName[2]," , ",loadDF[i,1]," peaks"),
          names = "",height = 25,color = "#333333",
          background = "#eeeeee")
      }
    })
    tx = append(tx,tx_tCons)
    
    
    trackRendered = TnT::TnTBoard(tx,
                             view.range = GenomicRanges::GRanges(elements[1],IRanges::IRanges(as.numeric(elements[2]),as.numeric(elements[3]))),
                             use.tnt.genome = T,
                             coord.range = IRanges::IRanges(window_load[1],window_load[2]))
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