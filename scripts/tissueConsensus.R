loadDF = tissueAvailability(input$tissue,input$peaks,input$regElement)
minCount = input$minCount
if(minCount<1){
  minCount = 1
}
fractionMotifsAssoc = input$fractionMotifsAssoc
if (fractionMotifsAssoc > 1) {
  fractionMotifsAssoc = 1
} else if(fractionMotifsAssoc<0){
  fractionMotifsAssoc = 0
}
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
    if(length(input$tba)>0){
      tba = system(paste0("tabix ",inputFolder,"TBA_consensus/tissue_TBA/",gsub("union",loadDF[i,1],as.character(loadDF[i,2]))," ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
      tba_split = lapply(strsplit(tba, "\t", fixed=TRUE),function(x)if(length(x)<9){c(x,"")}else{x})
      tba = data.frame(do.call(rbind, tba_split),stringsAsFactors = F)
      for(pvalue in input$pvalues) {
      # for(pvalue in c("0.1","0.05","0.01","0.001","0.0001","0.00001")) {
        pvalue_col = match(pvalue,c("0.1","0.05","0.01","0.001","0.0001","0.00001"))+3
        # tCons = cbind(tCons,unlist(lapply(tba[,pvalue_col],function(row){
        #   pfms = as.numeric(strsplit(row,"/")[[1]])
        #   paste0(annot$annot[pfms[!countPFMS[pfms]<minCount],"GeneSymbol"],collapse = " - ")
        # })))
        tCons = cbind(tCons,unlist(lapply(tba[,pvalue_col],function(row){
          pfms = as.numeric(unlist(lapply(strsplit(row,"/")[[1]],function(split){strsplit(split,"\\(")[[1]][1]})))
          fract_1000GP = gsub(")","",unlist(lapply(strsplit(row,"/")[[1]],function(split){strsplit(split,"\\(")[[1]][2]})))
          # paste0(paste0(annot$annot[pfms[!countPFMS[pfms]<minCount],"GeneSymbol"],"-",
          #               fract_1000GP[!countPFMS[pfms]<minCount]),collapse = "/")
          load(paste0(inputFolder,"data/motifs_CountsFreq_Tissue/",gsub(".bed.gz",".counts_freq.RData",gsub("union",loadDF[i,1],as.character(loadDF[i,2])))))
          validPFMS = countPFMS[pfms]>=minCount & table[,pvalue_col-1][pfms]<=fractionMotifsAssoc
          paste0(paste0(paste0(annot$annot[pfms[validPFMS],"GeneSymbol"]," (",annot$annot[pfms[validPFMS],"Code"],")"),collapse="/"),"%",
                 paste0(fract_1000GP[validPFMS],collapse = "/"))
        })))
      }
      colnames(tCons) = c(colnames(tCons)[1:3],paste0("TBA_pvalue<",input$pvalues))
      # colnames(tCons) = c(colnames(tCons)[1:3],paste0("TBA_pvalue<",c("0.1","0.05","0.01","0.001","0.0001","0.00001")))
    }
    tCons = cbind(tCons,paste0(c(as.character(loadDF[i,1]),as.character(loadDF[i,2]),"tissue"),collapse="%"))
    colnames(tCons) = c(colnames(tCons)[1:(length(colnames(tCons))-1)],"typeCRE")
    GRanges_tcons = makeGRangesFromDataFrame(tCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = TRUE)
    FeatureTrack(
      GRanges_tcons,
      tooltip = as.data.frame(GRanges_tcons),
      label = paste0("tissue consensus tracks - ",splitFileName[1]," - ",splitFileName[2]," , ",loadDF[i,1]," peaks"),
      names = "",height = 50,color = consensusColor(loadDF[i,2],loadDF[i,1]),
      background = "#eeeeee")
  }
})