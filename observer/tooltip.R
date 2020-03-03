# Observer that listens to changes in tooltip (click on a regions, Shiny.onInputChange() in javascript)
observeEvent(input$tooltipTable, {
  removeUI(
    ## pass in appropriate div id
    selector = '#tbaTab'
  )
  removeUI(
    ## pass in appropriate div id
    selector = '#regionTab'
  )
  
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
    } else {
      insertUI(
        selector = "#placeholder",
        where = "beforeBegin",
        ui = tags$div(
          boxPlus(width = 12,collapsible = T,closable = F,title="Gene info:",
                  dataTableOutput("tTable")),
          id="regionTab"
        )
      )
    }
    df_tba = datatable(df_tba,selection = 'none',filter = 'top',rownames = FALSE,extensions = 'Buttons',
                       options = list(paging = TRUE,
                                      searching = TRUE,
                                      dom = 'tBfrtip',
                                      buttons = list('copy',
                                                     list(extend='csv',
                                                          filename = 'tbaInfo'),
                                                     list(extend='excel',
                                                          filename = 'tbaInfo'),
                                                     list(extend='pdf',
                                                          filename= 'tbaInfo'))),
                       )
    tbaDF(df_tba)
    if(length(idx_tba)>0){
      if(length(input$tba)>0) {
        insertUI(
          selector = "#placeholder",
          where = "afterEnd",
          ui = tags$div(
            boxPlus(width = 6,title = "TBA info:",collapsible = T,closable = F,
                    dataTableOutput("tbaTable")),
            id="tbaTab"
          )
        )
        insertUI(
          selector = "#placeholder",
          where = "beforeBegin",
          ui = tags$div(
            boxPlus(width = 6,collapsible = T,closable = F,title="Region info:",
                    dataTableOutput("tTable")),
            id="regionTab"
          )
        )
      }
      df = df[-idx_tba,]
    }
    df = datatable(df,selection = 'none',rownames = FALSE,colnames = c("",""),extensions = 'Buttons',
                   options = list(paging = FALSE,
                                  searching = FALSE,
                                  dom="tBfrtip",
                                  buttons = list('copy',
                                                 list(extend='csv',
                                                             filename = 'geneInfo'),
                                                 list(extend='excel',
                                                      filename = 'geneInfo'),
                                                 list(extend='pdf',
                                                      filename= 'geneInfo'))))
    tooltipDF(df) # set reactiveVal to new value
    # df <- df[, list(cars=list(.SD)), by = list(mpg,cyl)]
  } else {
    if("entrezid"%in%df$label) {
      title = "Gene info:"
    } else {
      title = "Region info:"
    }
    insertUI(
      selector = "#placeholder",
      where = "beforeBegin",
      ui = tags$div(
        boxPlus(width = 12,collapsible = T,closable = F,title=title,
                dataTableOutput("tTable")),
        id="regionTab"
      )
    )
    df = datatable(df,selection = 'none',rownames = FALSE,colnames = c("",""),extensions = 'Buttons',
                   options = list(paging = FALSE,
                                  searching = FALSE,
                                  dom="tBfrtip",
                                  buttons = list('copy',
                                                 list(extend='csv',
                                                      filename = 'geneInfo'),
                                                 list(extend='excel',
                                                      filename = 'geneInfo'),
                                                 list(extend='pdf',
                                                      filename= 'geneInfo'))))
    tooltipDF(df) # set reactiveVal to new value
  }
})

# # Observer that listens the click on search by position or genes and reset the ata.frame with tooltip info
# observeEvent(c(input$Run,input$searchGenes), {
#   tooltipDF(NULL) # set reactiveVal to new value
#   tbaDF(NULL)
# })
