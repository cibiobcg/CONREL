# Observer that listens to changes in tooltip (click on a regions, Shiny.onInputChange() in javascript)
observeEvent(input$tooltipTable, {
  tooltip = input$tooltipTable
  df = data.frame(label = tooltip[names(tooltip)=="label"],
                  value = tooltip[names(tooltip)%in%c("value","value1")])
  if(length(input$tba)>0){
    idx_tba = grep("^TBA",df$label)
    df_tba = df[idx_tba,]
    if(nrow(df_tba)>0){
      res = lapply(rev(1:nrow(df_tba)),function(i) {
        if(df_tba[i,2]!="%"){
          data.frame(p.value = paste0("< ",formatC(as.numeric(gsub("TBA_pvalue.","",df_tba[i,1])), format = "e", digits = 2)),
                     TF_Symbol_and_Code = strsplit(strsplit(as.character(df_tba[i,2]),"%")[[1]][1],"/")[[1]],
                     ALL_1000GP = strsplit(strsplit(as.character(df_tba[i,2]),"%")[[1]][2],"/")[[1]])
        }
      })
    
    df_tba = do.call(rbind,res)
    }
    df_tba = datatable(df_tba,selection = 'none',filter = 'top',rownames = FALSE,
                       options = list(paging = TRUE,
                                      searching = TRUE))
    tbaDF(df_tba)
    if(length(idx_tba)>0){
      df = df[-idx_tba,]
    }
  }
  # df = as.data.frame(do.call(rbind,split(input$tooltipTable,sort(rep(1:(length(input$tooltipTable)/2),2)))))
  df = datatable(df,selection = 'none',rownames = FALSE,
                 options = list(paging = FALSE,
                                searching = FALSE,
                                dom="t"))
  tooltipDF(df) # set reactiveVal to new value
})

# Observer that listens the click on search by position or genes and reset the ata.frame with tooltip info
observeEvent(c(input$Run,input$searchGenes), {
  tooltipDF(NULL) # set reactiveVal to new value
  tbaDF(NULL)
})
