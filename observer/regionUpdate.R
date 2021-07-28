### Observed events for zoom/drag/update on track and update the input field
observeEvent(input$positionStart,ignoreInit = F,{
  #Update the input box of region
  updateTextInput(session, "region", value = paste0(input$chr,":",
                                                    prettyNum(input$positionStart,big.mark = ","),"-",
                                                    prettyNum(input$positionEnd,big.mark = ","))
  )
  #Update the URL for DAS server to obtain the sequence
  # url <- a("DNA sequence",
  #          href=paste0(dasServer,"?segment=",elements[1],":",input$positionStart,",",input$positionEnd), target="_blank")
  # updateActionButton(session,inputId = "seq_DAS")
  # output$seq_link <- renderUI({
  #   tagList(url)
  # })
})