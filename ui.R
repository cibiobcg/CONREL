source("./global.R")

tag = span(class = "logo-lg")
# child = tagList(div(img(class = "logo-lg", src="CONREL_LOGO_CORNER.svg",width = "100%", height = "auto",style="float: left;display:inline-block;margin-left:-25px;")),
# div(span(h5(hg_assembly,style="font-weight: 900;"),style="float: left;display:inline-block;margin-left:-20px;")))
child = tagList(div(img(class = "logo-lg", 
                        id='logo-lg', 
                        src="empty_corner.png",
                        width = "100%", 
                        height = "auto",
                        style="float: left;display:inline-block;margin-left:-25px;")),
                img(class='logo-sq',
                    id='logo-sq',
                    src = "empty_square.png"))
tag2 = tagAppendChild(tag, child)

ui <- dashboardPage(title = "CONREL",
  preloader = list(html = tagList(spin_1(), "Loading..."), color = "#3c8dbc"),
  # collapse_sidebar = FALSE,
  # enable_preloader = TRUE,
  # loading_duration = 2,
  md = TRUE,
  sk='green',
  
  ############  
  ### HEAD ###
  ############
  
  dashboardHeader(
    title = tagList(
      # tag = span(class = "logo-lg", "CONREL v1"), 
      tag2),
    leftUi = tagList(dropdownMenuOutput("legendMenu"))
  ),
  
  
  ############
  ### SIDE ###
  ############
  
  dashboardSidebar(
    sidebarMenu(id = "sideBar",
                # menuItemOutput("assembly"),
                menuItem("Change Organism/Assembly", tabName = "changeORG", icon = icon("random")),
                hr(),
                menuItem("Home", tabName = "home", icon = icon("home"),selected = T),
                #menuItem("Search", tabName = "search1", icon = icon("search")),
                #menuItemOutput("outSearch"),
                menuItem("Search", tabName = "search1",icon = icon("search")),
                conditionalPanel(condition='input.menu == "hidden2"',menuItem(NULL, tabName = "search2")),
                conditionalPanel(condition='input.menu == "hiddenError"',menuItem(NULL, tabName = "errorTab")),
                menuItemOutput("outGenome"),
                # hr(),uiOutput("seq_link"),
                # downloadButton('download',"Download the data",style="color: #fff; background-color: #006502; border-color: #2e6da4"),
                menuItemOutput("link"),hr(),
                menuItem("About", tabName = "about", icon = icon("info-circle")),
                menuItem("Help", tabName = "help", icon = icon("question-circle")),
                menuItem("Download", tabName = "singularity", icon = icon("box-open"))
    )
  ),
  
  
  ############
  ### BODY ###
  ############
  
  dashboardBody(
    # setBackgroundImage(src = 'background.png', shinydashboard = TRUE),
    tags$head(tags$style(HTML('
      .custom-ui-errors {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 15px;
        color: red;
        height: auto;
        background-image: url("error.png");
        background-repeat: no-repeat;
        padding-left: 25px;
        display: block;
      }
      .custom-ui-warnings {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 14px;
        color: #AA4400;
        height: auto;
        background-image: url("warning.png");
        background-repeat: no-repeat;
        padding-left: 25px;
        display: block;
      }
      .custom-info {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 12px;
        color: darkgrey;
        height: auto;
        display: block;
      }
      .control-label {
        color: #666 !important;
        font-size: 16px !important;
      }
      .btn {
        border: 2px solid #aaaaaa;
      }
      #boxsearch .box-header {
        padding: 0 !important;
        height: 0 !important;
      }
      #helpBox .box-header {
        background-color: #1ab394;
        color: white;
      }
      
      .styled-table {
        border-collapse: collapse;
        margin: 25px 0;
        font-size: 0.9em;
        font-family: sans-serif;
        min-width: 400px;
        box-shadow: 0 0 20px rgba(0, 0, 0, 0.15);
}
      .styled-table thead tr {
        background-color: #009879;
        color: #ffffff;
        text-align: left;
}
      .styled-table th,
      .styled-table td {
        padding: 12px 15px;
}
      .styled-table tbody tr {
        border-bottom: 1px solid #dddddd;
}
      .styled-table tbody tr:nth-of-type(even) {
        background-color: #f3f3f3;
}
      .styled-table tbody tr:last-of-type {
        border-bottom: 2px solid #009879;
}
      .styled-table tbody tr.active-row {
        font-weight: bold;
        color: #009879;
      }
.borderd-content {
  border: 1px solid #000000;
  border-radius: 4px;
  height: 300px;
  margin-top: 20px;
  position: relative;
}

.borderd-content .title {
  margin: -10px 0 0 10px;
  background: #fff;
  padding: 3px;
  font-size: 1.5em;
  display: inline-block;
  font-weight: bold;
  position: absolute;
}

.borderd-content .content {
  padding: 10px;
}
div.table-loading {
  border: 1px solid #000000;
  border-style: solid none none none;
}
    '))),
    useShinyjs(),
    tabItems(
      tabItem(
        tabName = "home",
        fluidPage(
          fluidRow(width=12,align = "center",
                   #h1(tags$b("Welcome in CONREL")),
                   img(class = "logo-lg", src="CONREL_LOGO_OG.svg",width = "50%", height = "auto"),
                   h3(tags$b("CONsensus Regulatory ELements"))),
          br(),br(),br(),
          fluidRow(column(width=6,offset=3,
                          p("CONREL, a genome browser that allows for the exploration of consensus regulatory
                                  elements at different levels of abstraction. The total binding affinity of transcription
                                  factors on whole consensus region sequences is here fully exploited to characterize and
                                  annotate functional properties of regulatory elements. CONREL can be used to explore genomic
                                  loci, genes or genomic regions of interest across different cell lines and tissues.",style="text-align: justify; text-justify: inter-word;")
                          ,br(),br()))
        )),
      
      
      
      tabItem(
        tabName = "search1",
        fluidRow(column(
          width=6,offset=3,
          actionLink(inputId = "example1",label = "EXAMPLE #1",style='padding:12px;'),
          actionLink(inputId = "example2",label = "EXAMPLE #2")
        )),
        fluidRow(
          column(
            width=6,offset=3,
            box(width = 12,
                title = "SEARCH BY CHROMOSOME POSITION",
                #status = "warning", solidHeader = TRUE,
                textInput("region", "Chromosome region:",value = ''),
                uiOutput("searchRegionError",class="custom-ui-errors"),
                actionButton("searchButtonPos", "Search")
            ),
            box(width = 12,
                title = "SEARCH BY GENE NAME",
                selectizeInput('genes', label = "Gene symbol:",
                               choices = NULL, options = list(
                                 placeholder = 'Select a gene',
                                 maxItems = 1,multiple = F, searchConjunction = 'and')),
                uiOutput("searchGeneError",class="custom-ui-errors"),
                actionButton("searchButtonGene", "Search"),
                conditionalPanel(condition='input.menu == "hidden"',
                                 selectInput(
                                   inputId="selectData",
                                   label=" ", selected = NULL,
                                   choices=c( "title" )
                                 )
                )
            )
          )
        )
      ),
      tabItem(
        tabName = "search2",
        fluidRow(column(
          width=12,offset=0,
          box(width=12,id='boxsearch',
              fluidRow(
                column(
                  width=6,
                  actionButton(width="100%","back_2_s1","Back", icon = icon("long-arrow-alt-left"))
                ),
                column(
                  width=6,
                  actionButton(width="100%","searchSetting", "Search", icon = icon("search"))
                )
              )
          )
        )),
        fluidRow(column(
          width=8,offset=0,
          box(width=12,title = "SETTINGS",closable = F,
              fluidRow(
                column(width=12,offset=0,align = "center",
                       radioGroupButtons(inputId="choice_track", label="select the gene level or the transcript level visualization for the genome browser:", 
                                         choices = c("gene","transcript"),width = "100%")
                )
              ),
              checkboxGroupInput(inputId = "peaks",
                                 label = "Peaks calls format of the regions",
                                 choiceNames = c('Narrow Peaks','Broad Peaks'),
                                 choiceValues = c("narrow","broad")),
              uiOutput("peaksError",class="custom-ui-errors"),
              selectizeInput('regElement',
                             label = "Regulatory element:",
                             choices = NULL,
                             options = list(placeholder = 'select 1 or more',
                                            maxItems = 1500,multiple = T, searchConjunction = 'and')),
              uiOutput("regElementError",class="custom-ui-errors"),
              # checkboxGroupInput(inputId = "consensus",
              #                    label = "Consensus regions track to visualize:",
              #                    choiceNames = c('Global consensus regions','Tissue consensus regions','Cell-lines consensus regions'),
              #                    choiceValues = c("global","tissue","cell-line"),
              #                    selected='global'),
              # conditionalPanel(condition = 'input.consensus.includes("tissue")',
              #                  selectizeInput('tissue', label = "Tissue:",
              #                                 choices = NULL,
              #                                 options = list(placeholder = 'select 1 or more',
              #                                                maxItems = 1500,multiple = T, searchConjunction = 'and'))),
              fluidRow(
                column(width=6,
                       checkboxGroupInput(inputId = "consensus",
                                          label = "Consensus regions track to visualize:",
                                          choiceNames = c('Global consensus regions','Tissue consensus regions','Cell line consensus regions'),
                                          choiceValues = c("global","tissue","cell-line")),
                       uiOutput("consensusError",class="custom-ui-errors")),
                column(width=6,
                       conditionalPanel(condition = 'input.consensus.includes("tissue") | input.consensus.includes("cell-line")',
                                        p("NOTE: The lists of tissues and the tree of cell-lines are built based on the peak format selected. Please select the peak format before any tissue or cell-line selection. Otherwise, by changing the peak format it will reset the corresponding selection.",class="custom-ui-warnings")),
                       conditionalPanel(condition = 'input.consensus.includes("tissue")',
                                        # p("Please select one or more tissues. NOTE: The list of tissues depends on the peak format selected, by changing the peak format, it will reset the tissue list.",class="custom-ui-warnings"),
                                        selectizeInput('tissue', label = "Tissue:",
                                                       choices = NULL,
                                                       options = list(placeholder = 'select 1 or more',
                                                                      maxItems = 1500,multiple = T, searchConjunction = 'and')),
                                        uiOutput("tissueError",class="custom-ui-errors")
                       ))),
              h3("Optional tracks to add:",align="center"),column(width=6,
                                                                  checkboxGroupInput(inputId = "snps",
                                                                                     label = "",
                                                                                     choiceName = c('Add SNPs track - dbSNP v151'),
                                                                                     choiceValues = c('snp'))),column(width=6,
                                                                                                                      checkboxGroupInput(inputId = "tss",
                                                                                                                                         label = "",
                                                                                                                                         choiceName = c('Add SwitchGear TSS track'),
                                                                                                                                         choiceValues = c('tss'))),
              fluidRow(
                column(width=6,
                       checkboxGroupInput(inputId = "tba",
                                          label = "Add TBA info to global and tissue tracks.",
                                          choiceName = c('TBA - Total Binding Affinity'),
                                          choiceValues = c('tba')),
                       conditionalPanel(condition = "input.tba == 'tba'",
                                        uiOutput("tbaError",class="custom-ui-errors"))
                ),
                conditionalPanel(condition = "input.tba == 'tba'",
                                 column(width=6,
                                        # p("Please select one or more p-values cutoff",class="custom-ui-warnings"),
                                        checkboxGroupInput(inputId = "pvalues",label = "p-value cutoff:",selected = "0.00001",
                                                           choices = c("0.00001","0.0001","0.001","0.01","0.05","0.1")))
                )
              ),
              fluidRow(
                column(width=6,
                       conditionalPanel(condition = "input.tba == 'tba'",
                                        h4("TF motifs parameters:"),
                                        numericInput(inputId = "minCount", label = "Motifs count cutoff:",
                                                     50, min = 1, max = NA, step = 1),
                                        numericInput(inputId = "fractionMotifsAssoc", label = "Fraction motifs association in total region cutoff:",
                                                     0.50, min = 0, max = 1, step = 0.01)))
              )
          )
        ),
        column(width=4,
               box(width=12,title="Cell-line Tree selection",closable = F,
                   uiOutput("cellLineError",class="custom-ui-errors"),
                   conditionalPanel(condition = '!input.consensus.includes("cell-line")',
                                    p("*Select 'Cell-lines consensus regions' in setting panel to be able to select and add cell-line tracks level to the genome browser",class="custom-info")),
                   conditionalPanel(condition = 'input.consensus.includes("cell-line")',
                                    shinyTree("tree", checkbox = TRUE,theme = 'proton',themeIcons = F,themeDots = F))
               ))),
        fluidRow(column(
          width=12,offset=0,
          box(width=12,id='boxsearch',
              fluidRow(
                column(
                  width=6,
                  actionButton(width="100%","back_2_s1_bottom","Back", icon = icon("long-arrow-alt-left"))
                ),
                column(
                  width=6,
                  actionButton(width="100%","searchSetting_bottom", "Search", icon = icon("search"))
                )
              )
          )
        ))
      ),
      
      tabItem(
        tabName = "genome",
        fluidPage(
          width = 12,
          TnT::TnTOutput("outcomp",height = 'auto'),hr(),
          tags$div(id = 'placeholder'),
        )
      ),
      
      tabItem(
        tabName = "errorTab",
        fluidPage(
          width = 12,
          p('ERROR')
        )
      ),
      
      tabItem(
        tabName = "spinning",
        fluidPage(
          width = 12,
          actionButton(width="100%","searchSetting_bottom", "Search", icon = icon("search"))
          
        )
      ),
      
      tabItem(
        tabName = "about",
        fluidPage(column(
          width = 12,offset=0,
          # h2("CONREL - CONsensus Regulatory ELements"),
          valueBox(href = "#about_tf",width=3,"5,424","TF DNA-binding sites motifs", color = "yellow", icon = icon("level-down-alt")),
          valueBox(href = "#about_tf",width=3,"1,710","Transcription factor", color = "purple", icon = icon("share")),
          valueBox(href = "#about_cell",width=3,"1,398","ChIP-seq data", color = "green", icon = icon("signal")),
          valueBox(href = "#about_cre",width=3,"1.5 million","genome regulatory elements", color = "blue", icon = icon("align-center")),
          column(width=5,
                 box(width=12,closable = F,title="Datasets sources",id = 'helpBox',
                     status = "danger", solidHeader = F,
                     includeHTML("www/about/about_dataset.html")
                 ),
                 # box(width=12,closable = F,title="liftOver GRCh38/hg38",
                 #     status = "danger", solidHeader = F,
                 #     includeHTML("www/about/about_hg38.html")
                 # )
          ),
          box(width=7,title='Workflow',id = 'helpBox',
              includeHTML("www/about/about.html")),
          
          box(width=12,title='ChIP-seq data',id = 'helpBox',
              includeHTML("www/about/about_cell.html")),
          box(width=12,title='Cis-regulatory elements',id = 'helpBox',
              includeHTML("www/about/about_cre.html")),
          box(width=12,title='TF and TFBS motifs',id = 'helpBox',
              includeHTML("www/about/about_tf.html")),
          box(width=12,title='Comparison resources',id = 'helpBox',
              includeHTML("www/about/about_comparison.html")),
          box(width=12,title='How to cite CONREL',id = 'helpBox',
              includeHTML("www/about/about_cite.html"))
        ))
      ),
      
      tabItem(
        tabName = "help",
        fluidPage(column(
          width = 12,offset=0,
          box(width=12,title='What is CONREL',id = 'helpBox',
              includeHTML("www/help/help.html")),
          box(width=12,title='Starting point',id = 'helpBox',
              includeHTML("www/help/help_start.html")),
          box(width=12,title='Search page',id = 'helpBox',
              includeHTML("www/help/help_search.html")),
          box(width=12,title='Setting page',id = 'helpBox',
              includeHTML("www/help/help_setting.html")),
          box(width=12,title='Browser page',id = 'helpBox',
              includeHTML("www/help/help_browser.html"))
        ))
      ),
      
      tabItem(
        tabName = "singularity",
        fluidPage(
          # column(
          #   width=12,
          #   box(width = 12,closable = F,title="Download singularity image",id = 'helpBox',
          #       status = "danger", solidHeader = F,
          #       includeHTML("www/download/singularity.html"),
          #       box(width=12,closable = F,collapsible = T,collapsed = T,title="How to install and run singularity image",
          #           status = "danger", solidHeader = F,
          #           includeHTML("www/download/singularity_info.html")
          #           
          #       )
          #   ),
          #   box(width = 12,closable = F,title="Download other data",id = 'helpBox',
          #       status = "danger", solidHeader = F,
          #       includeHTML("www/download/download.html"),
          #       box(width=12,closable = F,collapsible = T,collapsed = T,title="Data in the download folder",
          #           status = "danger", solidHeader = F,
          #           includeHTML("www/download/download_info.html")
          #       )
          #   )
          # )
          box(title = 'Download singularity image and data',id = 'helpBox',
              width=12,
              h2('Download resources'),
              includeHTML("www/download/download.html"),
              includeHTML("www/download/download_info.html"),
              hr(),
              h2('Download singularity app'),
              includeHTML("www/download/singularity.html"),
              h2('How to install and use it'),
              includeHTML("www/download/singularity_info.html")
          )
        )
      )
    )
  )
  
)