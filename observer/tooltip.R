# Observer that listens to changes in tooltip (click on a regions, Shiny.onInputChange() in javascript)
observeEvent(input$tooltipTable, {
  tooltip = input$tooltipTable
  df = data.frame(label = tooltip[names(tooltip)=="label"],
                  value = tooltip[names(tooltip)%in%c("value","value1")])
  # df = as.data.frame(do.call(rbind,split(input$tooltipTable,sort(rep(1:(length(input$tooltipTable)/2),2)))))
  tooltipDF(df) # set reactiveVal to new value
})

# Observer that listens the click on search by position or genes and reset the ata.frame with tooltip info
observeEvent(c(input$Run,input$searchGenes), {
  tooltipDF(NULL) # set reactiveVal to new value
})