output$outGenome <- renderMenu({
  menuItem("Genome browser", tabName = "genome", icon = icon("globe"))
})