output$downloadTab <- downloadHandler(
  filename = function(){paste0("consensusInfo_",Sys.Date(),".csv")},
  content = function(fname){
    write.csv(data.frame(label = input$tooltipTable[names(input$tooltipTable)=="label"],
                         value = input$tooltipTable[names(input$tooltipTable)%in%c("value","value1")]), fname,
              quote = F,sep = "\t",row.names = F)
  }
)

output$downloadAll <- downloadHandler(
  filename = function(){paste0("regionInfo_",Sys.Date(),".csv")},
  content = function(fname){
    write.csv(data.frame(label = input$tooltipTable[names(input$tooltipTable)=="label"],
                         value = input$tooltipTable[names(input$tooltipTable)%in%c("value","value1")]), fname,
              quote = F,sep = "\t",row.names = F)
  }
)