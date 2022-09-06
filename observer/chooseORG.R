# Show modal when button is clicked.
observeEvent(input$sideBar, ignoreInit = T, {
  if(input$sideBar=="changeORG") {
    showModal(dataModal())
    checkFolderAssembly()
    shinyjs::hide("loading_page")
    shinyjs::hide('ok1')
    shinyjs::hide('ok2')
    shinyjs::hide('ok3')
    shinyjs::hide('ok4')
  }
})

# a function in common with each of the observer

observeEvent(input$hg19,ignoreInit = T, {
  # Check that data object exists and is data frame.
  shinyjs::hide("main_content")
  shinyjs::show("loading_page")
  
  tagWorkInProgress = tagList(NULL)
  
  loadAssemblyData('hg19',session)
  
  removeModal()
  
})

observeEvent(input$hg38,ignoreInit = T, {
  # Check that data object exists and is data frame.
  shinyjs::hide("main_content")
  shinyjs::show("loading_page")
  
  # tagWorkInProgress <<- tagList(div(style="border:2px red dotted;",
  #                                   p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
  #                                     "ATTENTION! WORK IN PROGESS"),
  #                                   p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
  #                                     "The GRCh38 version of CONREL is online but not fully functional. Some TBA annotations are still on processing, then the TBA informations will not show for those CREs.")))
  tagWorkInProgress <<- tagList(NULL)
  loadAssemblyData('hg38',session)
  removeModal()
  
})

observeEvent(input$mm10,ignoreInit = T, {
  # Check that data object exists and is data frame.
  shinyjs::hide("main_content")
  shinyjs::show("loading_page")
  
  tagWorkInProgress <<- tagList(div(style="border:2px red dotted;",
                                    p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
                                      "ATTENTION! WORK IN PROGESS"),
                                    p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
                                      "The GRCm38 version of CONREL is not fully functional and still in progress. We are generating and upgrading all the data")))
  # tagWorkInProgress <<- tagList(NULL)
  loadAssemblyData('mm10',session,mouse=T)
  removeModal()
  
})
