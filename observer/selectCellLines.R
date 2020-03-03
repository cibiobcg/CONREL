observeEvent(input$peaks, {
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
})