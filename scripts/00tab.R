removeUI(
  ## pass in appropriate div id
  selector = '#tbaTab'
)
removeUI(
  ## pass in appropriate div id
  selector = '#regionTab'
)

if(length(input$tba)>0) {
  insertUI(
    selector = "#placeholder",
    where = "afterEnd",
    ui = tags$div(
      boxPlus(width = 12,title = "TBA info:",collapsible = T,closable = F,
              dataTableOutput("tbaTable")),
      id="tbaTab"
    )
  )
}

insertUI(
  selector = "#placeholder",
  where = "afterEnd",
  ui = tags$div(
    boxPlus(width = 12,collapsible = T,closable = F,title="Info:",
            dataTableOutput("tTable")),
    id="regionTab"
  )
)

