observeEvent(c(input$back_2_s1,input$back_2_s1_bottom), ignoreInit = T,{
  updateTabItems(session, "sideBar", "search1")
})

# observeEvent(input$back_2_s1_bottom, {
#   updateTabItems(session, "sideBar", "searchStart")
# })
