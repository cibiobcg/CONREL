observeEvent(input$sideBar, ignoreInit = T, {
  if(input$sideBar=="genome") {
    output$legendMenu <- renderMenu({
      dropdownBlock(id="legendMenu",badgeStatus = NULL,
                    icon = icon("sliders-h"),title="Colours Legend",
                    fluidRow(h4(if(input$choice_track=='gene'){"Gene track:"}else{"Transcript track:"},style="color:black;font-weight: bold;width:600px"),
                      HTML(paste0('<style type="text/css">
                            .tg  {border-collapse:collapse;border-spacing:0;color:black}
                            .tabglobal {width:300px}

                            .tgInside{text-align:left;vertical-align:top;width:50px}
                            </style>
                            <table class="tabglobal"><tbody><tr><td>
                            <table class="tg">
                            <tbody>',
                        if(input$choice_track=='gene') {unlist(lapply(1:12,function(x){
                          paste0('<tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:',mapColor[legend_idx][x],'"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">',levelsColor[legend_idx][x],'</td></tr>')
                        }))}else{unlist(lapply(1:12,function(x){
                          paste0('<tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:',mapColor_tx[legend_idx_tx][x],'"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">',levelsColor_tx[legend_idx_tx][x],'</td></tr>')
                        }))},'
                            </tbody>
                            </table></td></tr></tbody></table>')
                      ),h4("Consensus track:",style="color:black;font-weight: bold;"),
                      HTML(paste0('<style type="text/css">
                            .tg  {border-collapse:collapse;border-spacing:0;color:black}
                            .tabglobal {width:500px}

                            .tgInside{text-align:left;vertical-align:top;width:50px}
                            </style>
                            <table class="tabglobal"><tbody><tr><td>
                            <table class="tg">
                            <tbody>
                            <tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:#ff8c00"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">Global consensus - Promoter - Narrow Peaks</td></tr>
                            <tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:#8b4500"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">Global consensus - Promoter - Broad Peaks</td></tr>
                            <tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:#4169e1"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">Global consensus - Enhacner - Narrow Peaks</td></tr>
                            <tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:#27408b"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">Global consensus - Enhancer - Broad Peaks</td></tr>
                            <tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:#00ff7f"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">Global consensus - Active Enhancer - Narrow Peaks</td></tr>
                            <tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:#008b45"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">Global consensus - Active Enhancer - Broad Peaks</td></tr>
                            <tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:#333333"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">Tissue consensus</td></tr>
                            <tr><td><table class="tgInside"><tbody><tr style="height:8px"><td style="background-color:#000000"></td></tr></tbody></table></td><td style="padding-left:15px;padding-bottom:5px">Cell-line consensus</td></tr>
                            </tbody>
                            </table></td></tr></tbody></table>')
                      )
                    )
      )
    })
  }
  else {
    output$legendMenu <- NULL
  }
})