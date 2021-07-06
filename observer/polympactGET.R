#Observer for GET
observe({
  query <- parseQueryString(session$clientData$url_search)
  if ((!is.null(query[['rsID']]))&(!is.null(query[['tissue']]))&(!is.null(query[['peak']]))) {
    Sys.sleep(2)
    withProgress(message = 'Generating plot', value = 0, {
      incProgress(0.2,detail = "Searching variant position")
      
      polympact_rsID = query[['rsID']]
      polympact_tissue = query[['tissue']]
      polympact_tissue = c(strsplit(polympact_tissue,',')[[1]])
      polympact_peak = query[['peak']]
      polympact_peak = c(strsplit(polympact_peak,',')[[1]])
      polympact_cres = c("promoter","enhancer","active")
      
      res_rsPos = system(paste0('zgrep "',polympact_rsID,'\t" ',inputFolder,'dbSNP_v151_hg19/dbSNP_151_hg19_complete_TOPMED_ALL.bed.gz'),intern=T)
      if(length(res_rsPos)==1){
        snpWindow = 50
        res_rsPos = strsplit(res_rsPos,'\t')[[1]]
        chr_poly = as.integer(res_rsPos[1])
        pos_poly = as.integer(res_rsPos[2])
        
        
        output$peaksError <- NULL
        output$regElementError <- NULL
        output$consensusError <- NULL
        output$tissueError <- NULL
        output$tbaError <- NULL
        output$cellLineError <- NULL
        
        clickPos <<- TRUE
        region = c(paste0('chr',chr_poly,':',pos_poly-snpWindow,'-',pos_poly+snpWindow))
        # region = c("chrX:66,764,465-66,950,461")
        
        
        if(simpleDebug){print(paste0("Search by position: ",clickPos))}
        if(simpleDebug){print(paste0("Search by gene: ",clickGene))}
        
        
        if(simpleDebug){print("View genome tab")}
        source(file.path("scripts", "genomeMenu.R"),  local = TRUE)$value
        updateTabItems(session, "sideBar", "genome")
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
        
        
        
        
        # SNPs track
        incProgress(0.1,detail = "Loading SNPs track")
### Instead of loadSNPs function
        snps.filtered = system(paste0("tabix ",inputFolder,"dbSNP_v151_hg19/dbSNP_151_hg19_complete_TOPMED_ALL.bed.gz ",gsub("chr","",elements[1]),":",window_load[1],"-",window_load[2]),intern = T)
        snps.filtered = data.frame(do.call(rbind, strsplit(snps.filtered, "\t", fixed=TRUE)))
        ir <- IRanges(start = as.numeric(as.character(snps.filtered$X2)), width = 1)
        gpos <- GRanges(elements[1], ir)
        gpos$score <- rep(0.4,nrow(snps.filtered))
        idx_poly_snp = which(as.data.frame(gpos)$start == pos_poly)
        gpos$score[idx_poly_snp] <- 0.75
        snps_colors = rep("black",length(gpos))
        snps_colors[idx_poly_snp] <- 'red'
        snps.filtered$alleles = paste(snps.filtered$X4,snps.filtered$X5,sep=">")
        gtooltip <- snps.filtered[,c(3,6,7,8)]
        colnames(gtooltip) = c("rs ID", "1000G MAF", "TOPMED MAF", "Alleles")
        pt = PinTrack(gpos,value = gpos$score,color = snps_colors,label = "dbSNP v151", domain = c(0,1), tooltip = gtooltip,background = "#eeeeee")
        emptyTrack = TnT::BlockTrack(GRanges(elements[1],IRanges(0,1)),
                                     color = "#EEEEEE",background = "#EEEEEE",
                                     height = 15,label=NULL)
        track_snps = append(pt,emptyTrack)
###
        tx = append(tx,track_snps)
        
        ##############################################
        emptyTrack = TnT::BlockTrack(GRanges(elements[1],IRanges(0,1)),
                                     color = "#EEEEEE",background = "#EEEEEE",
                                     height = 12,label=NULL)
        tx = append(emptyTrack,tx)
        incStep = 0.6/2
        
        
        if(simpleDebug){print("Load global consensus")}
        
        loadDF = mappingInput(polympact_cres,polympact_peak)
        tx_gCons = lapply(1:nrow(loadDF),function(i){
          incProgress(incStep/nrow(loadDF),detail = "Loading global consensus regions")
          # gCons = fread(paste0(encodeFolder,loadDF[i,1],"Peaks/globalConsensus/",loadDF[i,2],".union.bed"),data.table = F)
          # gCons = gCons[gCons$V1==elements[1],]
          gCons = system(paste0("tabix ",inputFolder,"consensus/",loadDF[i,1],"Peaks/bgz_globalConsensus/",loadDF[i,2],".union.bed.gz ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
          gCons = data.frame(do.call(rbind, strsplit(gCons, "\t", fixed=TRUE)))
          if(nrow(gCons)==0){
            generateEmptyTrack(elements[1])
          } else {
            refDB = system(paste0("tabix ",inputFolder,"consensus/",loadDF[i,1],"Peaks/bgz_globalConsensus/",loadDF[i,2],".union.bed.gz ",elements[1],":",window_load[1],"-",window_load[2],
                                  " | bedtools intersect -wao -a stdin -b ",inputFolder,"data/referenceDB/reference_",loadDF[i,2],".gz"),intern = T)
            refDB = data.frame(do.call(rbind, strsplit(refDB, "\t", fixed=TRUE)))
            overlap = refDB %>% group_by(X1,X2,X3) %>% summarize(referenceDB = paste(sort(unique(X7)), collapse = ', '))
            
            gCons = cbind(gCons,overlap$referenceDB,paste0(c(as.character(loadDF[i,1]),as.character(loadDF[i,2]),"global"),collapse="%"))
            colnames(gCons) = c(colnames(gCons)[1:(length(colnames(gCons))-1)],"typeCRE")
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
        
        
        loadDF = tissueAvailability(polympact_tissue,polympact_peak,polympact_cres)
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
            refDB = system(paste0("tabix ",inputFolder,"consensus/",loadDF[i,1],"Peaks/bgz_tissueConsensus/",as.character(loadDF[i,2])," ",elements[1],":",window_load[1],"-",window_load[2],
                                  " | bedtools intersect -wao -a stdin -b ",inputFolder,"data/referenceDB/reference_",splitFileName[2],".gz"),intern = T)
            refDB = data.frame(do.call(rbind, strsplit(refDB, "\t", fixed=TRUE)))
            overlap = refDB %>% group_by(X1,X2,X3) %>% summarize(referenceDB = paste(sort(unique(X7)), collapse = ', '))
            
            tCons = cbind(tCons,overlap$referenceDB,paste0(c(as.character(loadDF[i,1]),as.character(loadDF[i,2]),"tissue"),collapse="%"))
            colnames(tCons) = c(colnames(tCons)[1:(length(colnames(tCons))-1)],"typeCRE")
            GRanges_tcons = makeGRangesFromDataFrame(tCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = TRUE)
            FeatureTrack(
              GRanges_tcons,
              tooltip = as.data.frame(GRanges_tcons),
              label = paste0("tissue consensus tracks - ",splitFileName[1]," - ",splitFileName[2]," , ",loadDF[i,1]," peaks"),
              names = "",height = 25,color = "#333333",
              background = "#eeeeee")
          }
        })
        tx = append(tx,tx_tCons)
        
        
        trackRendered = TnTBoard(tx,
                                 view.range = GRanges(elements[1],IRanges(as.numeric(elements[2]),as.numeric(elements[3]))),
                                 use.tnt.genome = T,
                                 coord.range = IRanges(window_load[1],window_load[2]))
        track(trackRendered)
      } else {
        updateTabItems(session, "sideBar", "errorTab")
        
      }
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