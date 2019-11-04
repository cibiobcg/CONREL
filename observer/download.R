output$downloadTab <- downloadHandler(
  filename = function(){paste0("consensusInfo_",Sys.Date(),".csv")},
  content = function(fname){
    df = data.frame(label = input$tooltipTable[names(input$tooltipTable)=="label"],
                    value = input$tooltipTable[names(input$tooltipTable)%in%c("value","value1")])
    idx_tba = grep("^TBA",df$label)
    df = df[-idx_tba,]
    write.csv(df, fname,
              quote = F,sep = "\t",row.names = F)
  }
)

output$downloadAll <- downloadHandler(
  filename = function(){paste0("regionInfo_",Sys.Date(),".csv")},
  content = function(fname){
    df = data.frame(label = input$tooltipTable[names(input$tooltipTable)=="label"],
                    value = input$tooltipTable[names(input$tooltipTable)%in%c("value","value1")])
    idx_tba = grep("^TBA",df$label)
    df = df[-idx_tba,]
    write.csv(df, fname,
              quote = F,sep = "\t",row.names = F)
  }
)

output$downloadTBA <- downloadHandler(
  filename = function(){paste0("tbaInfo_",Sys.Date(),".csv")},
  content = function(fname){
    df = data.frame(label = input$tooltipTable[names(input$tooltipTable)=="label"],
                    value = input$tooltipTable[names(input$tooltipTable)%in%c("value","value1")])
    idx_tba = grep("^TBA",df$label)
    df_tba = df[idx_tba,]
    res = lapply(rev(1:nrow(df_tba)),function(i) {
      data.frame(p.value = gsub("TBA_pvalue.","< ",df_tba[i,1]),TF = strsplit(as.character(df_tba[i,2])," - ")[[1]])
    })
    df = do.call(rbind,res)
    write.csv(df, fname,
              quote = F,sep = "\t",row.names = F)
  }
)