cellList = get_selected(input$tree, format = "slices")
loadDF = celllineAvailability(cellList,input$peaks,input$regElement)
tx_cCons = lapply(1:nrow(loadDF),function(i){
  splitFileName = strsplit(as.character(loadDF[i,2]),"\\.")[[1]]
  incProgress(incStep/nrow(loadDF),detail = "Loading cell-lines consensus regions")
  # cCons = fread(paste0(encodeFolder,loadDF[i,1],"Peaks/consensus/",as.character(loadDF[i,2])),data.table = F)
  # cCons = cCons[gCons$V1==elements[1],]
  cCons = system(paste0("tabix ",inputFolder,"consensus/",loadDF[i,1],"Peaks/bgz_consensus/",as.character(loadDF[i,2])," ",elements[1],":",window_load[1],"-",window_load[2]," | bedtools merge -i stdin"),intern = T)
  cCons = data.frame(do.call(rbind, strsplit(cCons, "\t", fixed=TRUE)))
  if(nrow(cCons)==0){
    generateEmptyTrack(elements[1])
  } else {
    refDB = system(paste0("tabix ",inputFolder,"consensus/",loadDF[i,1],"Peaks/bgz_consensus/",as.character(loadDF[i,2])," ",elements[1],":",window_load[1],"-",window_load[2]," | bedtools merge -i stdin",
                          " | bedtools intersect -wao -a stdin -b ",inputFolder,"data/referenceDB/reference_",splitFileName[3],".gz"),intern = T)
    refDB = data.frame(do.call(rbind, strsplit(refDB, "\t", fixed=TRUE)))
    overlap = refDB %>% dplyr::group_by(X1,X2,X3) %>% dplyr::summarize(referenceDB = paste(sort(unique(X7)), collapse = ', '))
    cCons = cbind(cCons,overlap$referenceDB,paste0(c(as.character(loadDF[i,1]),as.character(loadDF[i,2]),"cell"),collapse="%"))
    colnames(cCons) = c(colnames(cCons)[1:(length(colnames(cCons))-1)],"typeCRE")
    GRanges_ccons = GenomicRanges::makeGRangesFromDataFrame(cCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = TRUE)
    TnT::FeatureTrack(
      GRanges_ccons,
      tooltip = as.data.frame(GRanges_ccons),
      label = paste0("cell line consensus tracks - ",splitFileName[1]," - ",splitFileName[3]," , ",loadDF[i,1]," peaks"),
      names = "",height = 25,color = consensusColor(loadDF[i,2],loadDF[i,1]),
      background = "#eeeeee")
  }
})