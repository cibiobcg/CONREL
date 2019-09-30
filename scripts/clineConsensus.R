cellList = get_selected(input$tree, format = "slices")
loadDF = celllineAvailability(cellList,input$peaks,input$regElement)
tx_cCons = lapply(1:nrow(loadDF),function(i){
  splitFileName = strsplit(as.character(loadDF[i,2]),"\\.")[[1]]
  incProgress(incStep/nrow(loadDF),detail = "Loading cell-lines consensus regions")
  # cCons = fread(paste0(encodeFolder,loadDF[i,1],"Peaks/consensus/",as.character(loadDF[i,2])),data.table = F)
  # cCons = cCons[gCons$V1==elements[1],]
  cCons = system(paste0("tabix ",inputFolder,"consensus/",loadDF[i,1],"Peaks/bgz_consensus/",as.character(loadDF[i,2])," ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
  cCons = data.frame(do.call(rbind, strsplit(cCons, "\t", fixed=TRUE)))
  if(nrow(cCons)==0){
    generateEmptyTrack(elements[1])
  } else {
    GRanges_ccons = makeGRangesFromDataFrame(cCons,seqnames.field="X1",start.field = "X2",end.field = "X3")
    FeatureTrack(
      GRanges_ccons,
      tooltip = as.data.frame(GRanges_ccons),
      label = paste0("cell line consensus tracks - ",splitFileName[1]," - ",splitFileName[3]," , ",loadDF[i,1]," peaks"),
      names = "",height = 50,color = consensusColor(loadDF[i,2],loadDF[i,1]),
      background = "#eeeeee")
  }
})