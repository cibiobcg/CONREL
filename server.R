server <- function(input, output, session) {
  if(simpleDebug){cat(file=stderr(), "Start server\n")}
  # updateTF <- eventReactive(input$actbtn_unicode,{
  #   input$unicode
  # })
  # update_tfbs_value <- reactiveVal(tfbs_value)
  # update_tf_value <- reactiveVal(tf_value)
  # update_cell_value <- reactiveVal(cell_value)
  # update_cre_value <- reactiveVal(cre_value)
  # 
  # output$tfbs_valueRender <- renderValueBox({
  #   valueBox(width=3,
  #     update_tfbs_value(tfbs_value),
  #     "TF DNA-binding sites motifs", color = "yellow", icon = icon("level-down-alt")
  #   )
  # })
  # output$tf_valueRender <- renderValueBox({
  #   valueBox(width=3,
  #     update_tf_value(tf_value),
  #     "Transcription factor", color = "purple", icon = icon("share")
  #   )
  # })
  # output$cell_valueRender <- renderValueBox({
  #   valueBox(width=3,
  #     update_cell_value(cell_value),
  #     "ChIP-seq data", color = "green", icon = icon("signal")
  #   )
  # })
  # output$cre_valueRender <- renderValueBox({
  #   valueBox(width=3,
  #     update_cre_value(cre_value),
  #     "genome regulatory elements", color = "blue", icon = icon("align-center")
  #   )
  # })
  ### modal for choosing 
  source(file.path("observer", "chooseORG.R"),  local = TRUE)$value
  
  
  ### Observed events for back button in search tab
  source(file.path("observer", "backButton.R"),  local = TRUE)$value
  ### Observed events for zoom/drag/update on track and update the input field
  source(file.path("observer", "regionUpdate.R"),  local = TRUE)$value
  
  ###############
  ### TOOLTIP ###
  ###############
  
  # The tooltip are disable. Instead the click on a region trigger an event that generates a
  # data.frame and it is rendered below the track
  tooltipDF = reactiveVal(df_tooltip)
  # Same for tba table
  tbaDF = reactiveVal(df_tba)
  # Same for cell line list table
  cellDF = reactiveVal(df_cell)
  # Observer that listens to changes in tooltip (click on a regions, Shiny.onInputChange() in javascript)
  source(file.path("observer", "tooltip.R"),  local = TRUE)$value
  
  ##############
  ### LEGEND ###
  ##############
  source(file.path("observer", "legendMenu.R"), local = TRUE)$value
  
  
  ################
  ### ASSEMBLY ###
  ################
  # source(file.path("observer","assembly.R"), local = TRUE)$value
  
  
  ###############################
  ### SELECT CONSENSUS & GENE ###
  ###############################
  updateSelectizeInput(session, 'regElement', choices = c("promoter","enhancer","active enhancer"), server = TRUE)
  # The input gene element is run on the server!
  # It increase the speed by at least 50x (especially during first rendering)
  source(file.path("observer", "selectGene.R"),  local = TRUE)$value
  
  
  ####################################
  ### SELECT TISSUE AND CELL LINES ###
  ####################################
  
  # Observer that changes the tissue consensus list based on what is selected
  source(file.path("observer", "selectTissues.R"),  local = TRUE)$value
  
  # Observer that changes the cell lines consensus list based on what is selected
  source(file.path("observer", "selectCellLines.R"),  local = TRUE)$value

  #############
  ### MODAL ###
  #############
  # source(file.path("observer", "modal.R"),  local = TRUE)$value
  
  
  
  ##############
  ### SEARCH ###
  ##############
  
  # Observer that check the search input
  source(file.path("observer", "search.R"),  local = TRUE)$value
  
  
  #######################
  ### COMPOSITE TRACK ###
  #######################         
  # The track response to generation of track by position/gene searches
  track = reactiveVal(trackOut)
  source(file.path("observer", "comptrack.R"),  local = TRUE)$value
  source(file.path("observer", "comptrack_Examples.R"),  local = TRUE)$value
  # Observer that generates the Composite track by position
  #source(file.path("observer", "comptrack_region.R"),  local = TRUE)$value
  # Observer that generates the Composite track by gene selection
  #source(file.path("observer", "comptrack_gene.R"),  local = TRUE)$value
  ### Observed for genome load position given the position from the GET
  source(file.path("observer", "polympactGET.R"),  local = TRUE)$value
  
  ################
  ### DOWNLOAD ###
  ################
  
  # # source(file.path("observer", "download.R"),  local = TRUE)$value
  # output$downloadData <- downloadHandler(
  #   filename = function() {
  #     paste('genomeBrowser_v1.tar', sep='')
  #   },
  #   content = function(con) {
  #     file.copy("/shares/CIBIO-Storage/CIBIO/sharedRL/Projects/genomeBrowser/source/container/genomeBrowser_v1.tar", con)
  #   },
  #   contentType = "application/zip"
  # )
  
  
  #################
  ### RENDERING ###
  #################
  # Rendering of data.frame tooltip
  output$tTable = renderDataTable(
    tooltipDF(),
    options = list(paging = FALSE,
                   searching = FALSE,
                   dom="t")
  )
  
  # Rendering of data.frame tba
  output$tbaTable = renderDataTable(
    tbaDF()
    # options = list(paging = TRUE,
    #                searching = TRUE)
    #                # dom="t")
  )
  
  output$cellTable = renderDataTable(
    cellDF()
    # options = list(paging = FALSE,
    #                searching = FALSE,
    #                dom="t")
  )
  
  # Rendering of composite track
  output$outcomp = TnT::renderTnT({
    if(simpleDebug){cat(file=stderr(), "Rendering...\n")}
    track()
  })
  
  showModal(dataModal())
  shinyjs::hide("loading_page")
  shinyjs::hide('ok1')
  shinyjs::hide('ok2')
  shinyjs::hide('ok3')
  shinyjs::hide('ok4')
  checkFolderAssembly()
  
}

