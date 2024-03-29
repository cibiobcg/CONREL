loadDF = mappingInput(input$regElement,input$peaks)
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
    overlap = refDB %>% dplyr::group_by(X1,X2,X3) %>% dplyr::summarize(referenceDB = paste(sort(unique(X7)), collapse = ', '))
    if(length(input$tba)>0 & file.exists(paste0(inputFolder,"TBA_consensus/global_TBA/",loadDF[i,2],".",loadDF[i,1],".bed.gz"))){
      tba = system(paste0("tabix ",inputFolder,"TBA_consensus/global_TBA/",loadDF[i,2],".",loadDF[i,1],".bed.gz ",elements[1],":",window_load[1],"-",window_load[2]),intern = T)
      tba_split = lapply(strsplit(tba, "\t", fixed=TRUE),function(x)if(length(x)<9){c(x,"")}else{x})
      tba = data.frame(do.call(rbind, tba_split),stringsAsFactors = F)
      for(pvalue in input$pvalues) {
        # for(pvalue in c("0.1","0.05","0.01","0.001","0.0001","0.00001")) {
        pvalue_col = match(pvalue,c("0.1","0.05","0.01","0.001","0.0001","0.00001"))+3
        gCons = cbind(gCons,unlist(lapply(tba[,pvalue_col],function(row){
          pfms = as.numeric(unlist(lapply(strsplit(row,"/")[[1]],function(split){strsplit(split,"\\(")[[1]][1]})))
          fract_1000GP = gsub(")","",unlist(lapply(strsplit(row,"/")[[1]],function(split){strsplit(split,"\\(")[[1]][2]})))
          # paste0(paste0(annot$annot[pfms[!countPFMS[pfms]<minCount],"GeneSymbol"],"-",
          #               fract_1000GP[!countPFMS[pfms]<minCount]),collapse = "/")
          load(paste0(inputFolder,"data/motifs_CountsFreq_Global/",loadDF[i,2],".",loadDF[i,1],".counts_freq.RData"))
          validPFMS = countPFMS[pfms]>=minCount & table[,pvalue_col-1][pfms]<=fractionMotifsAssoc
          
          paste0(paste0(paste0(annot$annot[pfms[validPFMS],"GeneSymbol"]," (",annot$annot[pfms[validPFMS],"Code"],")"),collapse="/"),"%",
                 paste0(fract_1000GP[validPFMS],collapse = "/"))
          
        })))
      }
      colnames(gCons) = c(colnames(gCons)[1:3],paste0("TBA_pvalue<",input$pvalues))
      # colnames(gCons) = c(colnames(gCons)[1:3],paste0("TBA_pvalue<",c("0.1","0.05","0.01","0.001","0.0001","0.00001")))
    }
    gCons = cbind(gCons,overlap$referenceDB,paste0(c(as.character(loadDF[i,1]),as.character(loadDF[i,2]),"global"),collapse="%"))
    colnames(gCons) = c(colnames(gCons)[1:(length(colnames(gCons))-1)],"typeCRE")
    GRanges_gcons = GenomicRanges::makeGRangesFromDataFrame(gCons,seqnames.field="X1",start.field = "X2",end.field = "X3",keep.extra.columns = T)
    TnT::FeatureTrack(
      GRanges_gcons,
      tooltip = as.data.frame(GRanges_gcons),
      label = paste0("global consensus tracks - ",gsub("active","active enhancer",loadDF[i,2])," , ",loadDF[i,1]," peaks"),
      names = "",height = 25,color = consensusColor(loadDF[i,2],loadDF[i,1]),
      background = "#eeeeee")
  }
})
