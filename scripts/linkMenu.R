output$link <- renderMenu({
  menuItem("Tools", tabName = "link",icon = icon("tools"),startExpanded = T,
           menuSubItem("DNA sequence", tabName = "dna", icon = icon("dna")),
                       # uiOutput("seq_link")
           menuSubItem("Download data", tabName = "dl", icon = icon("download")),
           # This for deselect when DNA or DL are selected. In this way you can click ultiple times on the same Item
           conditionalPanel(condition='input.menu == "hidden"',menuSubItem(NULL, tabName = "empty"))
                       # downloadButton('download',"Download the data",style="color: #fff; background-color: #006502; border-color: #2e6da4")
  )
})