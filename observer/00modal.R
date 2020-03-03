observeEvent(input$dna, {
  # if(input$sideBar == "dna") {
    updateTabItems(session, "sideBar", "genome")
    showModal(
      modalDialog(title = "DNA sequence",easyClose = T,align="center",
                  footer = a("More info about DAS project",
                             href="http://biodas.open-bio.org/wiki/Main_Page", target="_blank"),
                  p("To retrieve the DNA sequence of the region selected:"),
                  # uiOutput("seq_link"),
                  actionButton("SequenceDAS","DNA sequence - DAS server",icon = icon("file-code"),
                               # style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                               onclick = paste0("window.open('",dasServer,"?segment=",elements[1],":",input$positionStart,",",input$positionStop,"', '_blank')")),
                  actionButton("SequenceFasta","DNA sequence - fasta file format",icon = icon("file-alt"),
                               # style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                               onclick = paste0("window.open('http://togows.org/api/ucsc/hg19/",elements[1],":",input$positionStart,"-",input$positionStop,".fasta', '_blank')"))
                  )
    )
    updateTabItems(session,"sideBar",selected = "empty")
  })
observeEvent(input$dl, {
    updateTabItems(session, "sideBar", "genome")
    showModal(
      modalDialog(title = "Download data",easyClose = T,align="center",
                  conditionalPanel(condition = 'input.tooltipTable !== undefined',
                                   downloadButton('downloadTab',"Download the 'Region info' displayed",class = "asdasd"),
                                   hr(),
                                   downloadButton('downloadTBA',"Download the 'TBA info' displayed",class = "asdasd")
                                   ),
                  hr(),
                  disabled(downloadButton('downloadAll',"Download all infos of the region visualized")))
    )
    updateTabItems(session,"sideBar",selected = "empty")
  })