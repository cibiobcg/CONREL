loadDF = tissueAvailability(input$tissue,input$peaks,input$regElement)
tx_tCons = lapply(1:nrow(loadDF),function(i){
  splitFileName = strsplit(as.character(loadDF[i,2]),"\\.")[[1]]
  incProgress(incStep/nrow(loadDF),detail = "Loading tissue consensus regions")
  # tCons = fread(paste0(encodeFolder,loadDF[i,1],"Peaks/tissueConsensus/",as.character(loadDF[i,2])),data.table = F)
  # tCons = tCons[tCons$V1==elements[1],]
  tCons = system(paste0("tabix ",encodeFolder,loadDF[i,1],"Peaks/bgz_tissueConsensus/",as.character(loadDF[i,2])," ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
  tCons = data.frame(do.call(rbind, strsplit(tCons, "\t", fixed=TRUE)))
  if(nrow(tCons)==0){
    generateEmptyTrack(elements[1])
  } else {
    if(length(input$tba)>0){
      tba = system(paste0("tabix ",TBA_folder,"TBA_consensus/tissue_TBA/",gsub("union",loadDF[i,1],as.character(loadDF[i,2]))," ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
      tba_split = lapply(strsplit(tba, "\t", fixed=TRUE),function(x)if(length(x)<9){c(x,"")}else{x})
      tba = data.frame(do.call(rbind, tba_split),stringsAsFactors = F)
        for(pvalue in input$pvalues) {
        pvalue_col = match(pvalue,c("0.1","0.05","0.01","0.001","0.0001","0.00001"))+3
        tCons = cbind(tCons,unlist(lapply(tba[,pvalue_col],function(row){
          pfms = as.numeric(strsplit(row,"/")[[1]])
          paste0(annot$annot[pfms[!countPFMS[pfms]<minCount],"GeneSymbol"],collapse = " - ")
        })))
      }
      colnames(tCons) = c(colnames(tCons)[1:3],paste0("TBA_pvalue<",input$pvalues))
      # gCons$'TBA_pvalue<0.00001' = unlist(lapply(tba$X9,function(row){
      #   pfms = as.numeric(strsplit(row,"/")[[1]])
      #   paste0(annot$annot[pfms[!countPFMS[pfms]<minCount],"GeneSymbol"],collapse = " - ")
      # }))
    }
    GRanges_tcons = makeGRangesFromDataFrame(tCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = TRUE)
    FeatureTrack(
      GRanges_tcons,
      tooltip = as.data.frame(GRanges_tcons),
      label = paste0("tissue consensus tracks - ",splitFileName[1]," - ",splitFileName[2]," , ",loadDF[i,1]," peaks"),
      names = "",height = 50,color = consensusColor(loadDF[i,2],loadDF[i,1]),
      background = "#eeeeee")
  }
})