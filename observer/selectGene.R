# The input gene element is run on the server!
# It increase the speed by at least 50x (especially during first rendering)
# observeEvent(input$selectData, ignoreInit = T,{
#   updateSelectizeInput(session, "genes",
#                        choices = as.character(genesEns$gene_symbol),
#                        server = TRUE,selected = NULL)
# })