# Observer that generates the Composite track by position
observeEvent(input$searchButtonPos, ignoreInit = T,{
  output$searchRegionError <- NULL
  output$searchGeneError <- NULL
  if(!grepl(paste0("^chr([1-9]|1[0-9]|2[0-2]|[MXY]):[0-9][0-9,]*-[0-9][0-9,]*$"),input$region)){
    output$searchRegionError <- renderUI({
      p("Please provide a region in valid format. e.g. chr7:139424940-141784100")
    })
  } else {
    el = convertPosition(input$region,TRUE)
    chr = elements[1]
    start = as.numeric(el[2])
    end = as.numeric(el[3])
    if (end<=start) {
      output$searchRegionError <- renderUI({
        p("Please provide a region in valid format. END position need to be greater than START position")
      })
    } else if(!(data.table::between(start, 0, hg19[which(hg19$V1==chr),3]) | data.table::between(end, 0, hg19[which(hg19$V1==chr),2]))) {
      output$searchRegionError <- renderUI({
        p("Please provide a region in valid format. Positions provided are outside chromosome boundaries")
      })
    } else if(checkValidRegion(input$region)){
      output$searchRegionError <- renderUI({
        p("There is nothing to show in the region provided")
      })
    } else {
      if(simpleDebug){print("Change search tab")}
      clickPos <<- TRUE
      updateTabItems(session, "sideBar", "search2")
    }
  }
})

observeEvent(input$searchButtonGene, ignoreInit = T,{
  output$searchRegionError <- NULL
  output$searchGeneError <- NULL
  if(input$genes == ""){
    output$searchGeneError <- renderUI({
      p("Please provide a gene")
    })
  } else if (input$genes %in% ensembldb::genes(EnsDb)$symbol) {
    if(simpleDebug){print("Change search tab")}
    clickGene <<- TRUE
    updateTabItems(session, "sideBar", "search2")
  } else {
    output$searchGeneError <- renderUI({
      p("Please select a gene from the list")
    })
  }
})