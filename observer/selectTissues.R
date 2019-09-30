# Observer that changes the tissue consensus list based on what is selected
observeEvent(c(input$consensus,input$peaks), {
  tissue = c()
  n = unique(narrow_cell$tissueMapped)
  b = unique(broad_cell$tissueMapped)
  
  if (length(input$peaks)>1) {
    # tissue = c(paste0(intersect(n,b)," (Narrow and Broad peak)"),
    #            paste0(setdiff(n,b)," (Narrow peak)"),
    #            paste0(setdiff(b,n)," (Broad peak)"))
    tissue = unique(c(b,n))
  } else {
    if('narrow' %in% input$peaks) {
      # tissue = c(tissue,paste0(n," (Narrow peak)"))
      tissue = unique(n)
    }
    if ('broad' %in% input$peaks) {
      # tissue = c(tissue,paste0(n," (Broad peak)"))
      tissue = unique(b)
    }
  }
  tissue = tissue[tissue %in% tissueOG]
  updateSelectizeInput(session, 'tissue', choices = sort(unique(tissue)), server = TRUE)
})