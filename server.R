server <- function(input, output, session) {
    if(simpleDebug){print("Start server")}
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
    # Observer that listens to changes in tooltip (click on a regions, Shiny.onInputChange() in javascript)
    source(file.path("observer", "tooltip.R"),  local = TRUE)$value
    
    
    ###################
    ### SELECT GENE ###
    ###################
    updateSelectizeInput(session, 'regElement', choices = c("promoter","enhancer","active enhancer"), server = TRUE)
    # The input gene element is run on the server!
    # It increase the speed by at least 50x (especially during first rendering)
    source(file.path("observer", "selectGene.R"),  local = TRUE)$value
    
    
    #####################
    ### SELECT TISSUE ###
    #####################
    
    # Observer that changes the tissue consensus list based on what is selected
    source(file.path("observer", "selectTissues.R"),  local = TRUE)$value
    
    
    #############
    ### MODAL ###
    #############
    source(file.path("observer", "modal.R"),  local = TRUE)$value
    
    
    #######################
    ### COMPOSITE TRACK ###
    #######################         
    # The track response to generation of track by position/gene searches
    track = reactiveVal(trackOut)
    # Observer that generates the Composite track by position
    source(file.path("observer", "comptrack_region.R"),  local = TRUE)$value
    # Observer that generates the Composite track by gene selection
    source(file.path("observer", "comptrack_gene.R"),  local = TRUE)$value
    
    
    ################
    ### DOWNLOAD ###
    ################
    
    source(file.path("observer", "download.R"),  local = TRUE)$value
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('genomeBrowser_v1.tar', sep='')
      },
      content = function(con) {
        file.copy("/shares/CIBIO-Storage/CIBIO/sharedRL/Projects/genomeBrowser/source/container/genomeBrowser.tar", con)
      },
      contentType = "application/zip"
    )

    
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
    
    # Rendering of the tree for cell-line consensus
    output$tree <- renderTree({
        if(length(input$peaks)==0){
            
        } else if(length(input$peaks)>1){
            union_cell = unique(rbind(narrow_cell,broad_cell))
            idx_tissue = tissue[tissue%in%union_cell$tissueMapped]
            list_tissue = lapply(idx_tissue,function(x){
                idx_name = union_cell[union_cell$tissueMapped==x,"Biosample.term.name"]
                list_cell = lapply(idx_name,function(y){y})
                names(list_cell) = idx_name
                # if('tissue' %in% input$consensus){
                #     if(x %in% input$tissue){
                #         structure(list_cell, stopened=TRUE, stselected=TRUE)
                #     } else{
                #         structure(list_cell, stopened=FALSE)
                #     }
                # } else {
                    structure(list_cell, stopened=FALSE)
                # }
            })
            names(list_tissue) = idx_tissue
            outList = list('Cell-line consensus' = list_tissue)
            attr(outList[[1]],"stopened")=TRUE
            outList
        } else if(input$peaks=='narrow'){
            idx_tissue = tissue[tissue%in%narrow_cell$tissueMapped]
            list_tissue = lapply(idx_tissue,function(x){
                idx_name = narrow_cell[narrow_cell$tissueMapped==x,"Biosample.term.name"]
                list_cell = lapply(idx_name,function(y){y})
                names(list_cell) = idx_name
                # if('tissue' %in% input$consensus){
                #     if(x %in% input$tissue){
                #         structure(list_cell, stopened=TRUE, stselected=TRUE)
                #     } else{
                #         structure(list_cell, stopened=FALSE)
                #     }
                # } else {
                    structure(list_cell, stopened=FALSE)
                # }
            })
            names(list_tissue) = idx_tissue
            outList = list('Cell-line consensus' = list_tissue)
            attr(outList[[1]],"stopened")=TRUE
            outList
        } else if(input$peaks=='broad'){
            idx_tissue = tissue[tissue%in%broad_cell$tissueMapped]
            list_tissue = lapply(idx_tissue,function(x){
                idx_name = broad_cell[broad_cell$tissueMapped==x,"Biosample.term.name"]
                list_cell = lapply(idx_name,function(y){y})
                names(list_cell) = idx_name
                # if('tissue' %in% input$consensus){
                #     if(x %in% input$tissue){
                #         structure(list_cell, stopened=TRUE, stselected=TRUE)
                #     } else{
                #         structure(list_cell, stopened=FALSE)
                #     }
                # } else {
                    structure(list_cell, stopened=FALSE)
                # }
            })
            names(list_tissue) = idx_tissue
            outList = list('Cell-line consensus' = list_tissue)
            attr(outList[[1]],"stopened")=TRUE
            outList
        }
    })
    
    # Rendering of composite track
    output$outcomp = renderTnT({
        if(simpleDebug){print("Rendering...")}
        validate(
            if(clickPos){validateRegion(input,T)},
            if(clickGene){validateGenes(input,T)},
            validateInput(input)
        )
        source(file.path("scripts", "linkMenu.R"),  local = TRUE)$value
        
        track()
    })
    
    # observeEvent(c(input$Run,input$searchGenes),ignoreInit = T, {
    #     if(TRUE) {
                #         output$messageMenu <- renderMenu({
                #             dropdownMenu(type = "notifications",.list = v
                #             )
                #         })
    #         
    #     }
    # })    
    
}
