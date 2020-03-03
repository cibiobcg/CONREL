source("./global.R")

ui <- dashboardPagePlus(
  collapse_sidebar = FALSE,enable_preloader = TRUE,loading_duration = 2,
  md = TRUE,
  skin='green',
  
  ############  
  ### HEAD ###
  ############
  
  dashboardHeaderPlus(
    title = tagList(
      span(class = "logo-lg", "CONREL v1"), 
      img(src = "https://image.flaticon.com/icons/svg/268/268752.svg")),
    left_menu = tagList(dropdownMenuOutput("toolsMenu"))
  ),
  
  
  ############
  ### SIDE ###
  ############
  
  dashboardSidebar(
    sidebarMenu(id = "sideBar",
                menuItem("Home", tabName = "home", icon = icon("home"),selected = T),
                #menuItem("Search", tabName = "search1", icon = icon("search")),
                #menuItemOutput("outSearch"),
                menuItem("Search", tabName = "search1",icon = icon("search")),
                conditionalPanel(condition='input.menu == "hidden2"',menuItem(NULL, tabName = "search2")),
                menuItemOutput("outGenome"),
                # hr(),uiOutput("seq_link"),
                # downloadButton('download',"Download the data",style="color: #fff; background-color: #006502; border-color: #2e6da4"),
                menuItemOutput("link"),hr(),
                menuItem("Help", tabName = "help", icon = icon("info-circle")),
                menuItem("Download", tabName = "singularity", icon = icon("box-open"))
    )
  ),
  
  
  ############
  ### BODY ###
  ############
  
  dashboardBody(tags$head(tags$style(HTML('
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
    '))),
                useShinyjs(),
                tabItems(
                  
                  
                  
                  tabItem(
                    tabName = "home",
                    fluidPage(
                      fluidRow(width=12,align = "center",
                               h1(tags$b("Welcome in CONREL")),
                               h3(tags$b("CONsensus Regulatory ELements"))),
                      br(),br(),br(),
                      fluidRow(column(width=6,offset=3,
                                      p("CONREL, a genome browser that allows for the exploration of consensus regulatory
                                  elements at different levels of abstraction. The total binding affinity of transcription
                                  factors on whole consensus region sequences is here fully exploited to characterize and
                                  annotate functional properties of regulatory elements. CONREL can be used to explore genomic
                                  loci, genes or genomic regions of interest across different cell lines and tissues.")
                                  ))
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
                            textInput("region", "Chromosome region:"),
                            uiOutput("searchRegionError",class="custom-ui-errors"),
                            actionButton("searchButtonPos", "Search by position")
                        ),
                        box(width = 12,
                            title = "SEARCH BY GENE NAME",
                            selectizeInput('genes', label = "Gene symbol:",
                                           choices = NULL, options = list(
                                             placeholder = 'Select a gene',
                                             maxItems = 1,multiple = F, searchConjunction = 'and')),
                            uiOutput("searchGeneError",class="custom-ui-errors"),
                            actionButton("searchButtonGene", "Search by gene"),
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
                      boxPlus(width=12,
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
                      boxPlus(width=12,title = "SETTINGS",closable = F,
                              fluidRow(
                                column(width=12,offset=0,align = "center",
                                       radioGroupButtons(inputId="choice_track", label="Visualization: select the details of visualization for the genome browser:", 
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
                              h3("Optional tracks to add:",align="center"),
                              checkboxGroupInput(inputId = "snps",
                                                 label = "Add SNPs track.",
                                                 choiceName = c('SNPs - dbSNP v151'),
                                                 choiceValues = c('snp')),
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
                                                 p("Please select one or more p-values cutoff",class="custom-ui-warnings"),
                                                        checkboxGroupInput(inputId = "pvalues",label = "p-value cutoff:",
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
                           boxPlus(width=12,title="Cell-line Tree selection",closable = F,
                                   uiOutput("cellLineError",class="custom-ui-errors"),
                                   conditionalPanel(condition = '!input.consensus.includes("cell-line")',
                                                    p("*Select 'Cell-lines consensus regions' in setting panel to be able to select and add cell-line tracks level to the genome browser",class="custom-info")),
                                   conditionalPanel(condition = 'input.consensus.includes("cell-line")',
                                                    shinyTree("tree", checkbox = TRUE,theme = 'proton',themeIcons = F,themeDots = F))
                                   ))),
                    fluidRow(column(
                      width=12,offset=0,
                      boxPlus(width=12,
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
                  
                  # tabItem(
                  #   tabName = "input",
                  #   fluidRow(
                  #     boxPlus(width=5,title = "Input Tab",closable = F,
                  #             status = "info", solidHeader = TRUE
                  #             
                  #             
                  #     ),
                  #     box(width=3,title = "Cell-lines Tree selection",
                  #         status = "info", solidHeader = TRUE,
                  #         conditionalPanel(condition = 'input.consensus.includes("cell-line")',
                  #                          shinyTree("tree", checkbox = TRUE,theme = 'proton',themeIcons = F,themeDots = F)),
                  #         uiOutput("cellWarning",class="custom-ui-warnings"),
                  #     ),
                  #     box(width = 4,title = "Search Tab",
                  #         status = "warning", solidHeader = TRUE,
                  #         textInput("region", "chromosome region", "chr7:139,424,940-141,784,100"),
                  #         actionButton("Run", "Search by position"),hr(),
                  #         selectizeInput('genes', label = "Gene symbol:",
                  #                        choices = NULL, options = list(
                  #                          placeholder = 'Select a gene',
                  #                          maxItems = 1,multiple = F, searchConjunction = 'and')),
                  #         actionButton("searchGenes", "Search by gene"),
                  #         conditionalPanel(condition='input.menu == "hidden"',
                  #                          selectInput(
                  #                            inputId="selectData",
                  #                            label=" ", selected = NULL,
                  #                            choices=c( "title" )
                  #                          )
                  #         )
                  #         
                  #     )
                  #   )
                  # ),
                  
                  
                  
                  tabItem(
                    tabName = "genome",
                    fluidPage(
                      width = 12,
                      TnTOutput("outcomp",height = 'auto'),hr(),
                      tags$div(id = 'placeholder'),
                    )
                  ),
                  
                  
                  
                  tabItem(
                    tabName = "help",
                    fluidPage(column(
                      width = 6,offset=3,
                      includeHTML("www/help.html")
                      ))
                  ),
                  
                  tabItem(
                    tabName = "singularity",
                    fluidPage(
                      width = 12,
                      includeHTML("www/singularity_head.html"),
                      # h1("Download singularity image"),
                      # br(),
                      # h2("Download"),
                      br(),
                      # p("Download a singuarity image to run this shiny app on your local server."),
                      p(HTML('&emsp;'),"Download link:    ",downloadLink('downloadData', 'genomeBrowser_v1.tar')," 17 Oct 2019 (~30GB)"),
                      includeHTML("www/singularity_body.html")
                      # br(),
                      # h2("Instructions"),
                      # br(),
                      # h4("1. Unpack TAR archive"),
                      # box(width=12,"tar -xvf genomeBrowser.tar"),
                      # br(),
                      # h4("2. Prepare shiny-server configuration file"),
                      # p("You will first need to generate a custom configuration for your user, and it will give you instructions for usage:"),
                      # box(width=12,
                      #     p("$ /bin/bash prepare_conf.sh"),
                      #     br(),
                      #     code(
                      #       p("Steps:"),
                      #       p(HTML('&emsp;'),"----------------------------------------------------------------------"),
                      #       p(HTML('&emsp;'),"1. Use this script to prepare your shiny-server.conf (configuration)"),
                      #       
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"/bin/bash prepare_template.sh"),
                      #       
                      #       p(HTML('&emsp;'),"----------------------------------------------------------------------"),
                      #       p(HTML('&emsp;'),"2. If needed, you can provide the following arguments"),
                      #       
                      #       p(HTML('&emsp;'),"Commands:"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"help: show help and exit"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"start: the generation of your config"),
                      #       
                      #       p(HTML('&emsp;'),"Options:"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"--port:  the port for the application (e.g., shiny default is 3737)"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"--user:  the user for the run_as directive in the shiny configuration"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"--base: base folder with applications"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"--logs: temporary folder with write for logs (not required)"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"--disable-index: disable directory indexing"),
                      #       
                      #       p(HTML('&emsp;'),"----------------------------------------------------------------------"),
                      #       p(HTML('&emsp;'),"3. Make sure Singularity is loaded, and run the container using "),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"the commands shown by the template.")
                      #     )
                      # ),
                      # p("When you add 'start' it will do the generation. Here we don't supply any arguments so that they are randomly generated."),
                      # box(width=12,
                      #     p("$ /bin/bash prepare_template.sh start"),
                      #     br(),
                      #     code(
                      #       p("Generating shiny configuration..."),
                      #       p("port: 9870"),
                      #       p("logs: /tmp/shiny-server.gG1X2Z"),
                      #       p("base: /srv/shiny-server/shiny_genomeBrowser"),
                      #       p("Server logging will be in /tmp/shiny-server.gG1X2Z"),
                      #       
                      #       p("To run your server:"),
                      #       
                      #       p(HTML('&emsp;'),"module load singularity/2.4.6"),
                      #       p(HTML('&emsp;'),"singularity run --bind /tmp/shiny-server.gG1X2Z/logs:/var/log/shiny \ "),
                      #       p(HTML('&emsp;'),"--bind /tmp/shiny-server.gG1X2Z/lib:/var/lib/shiny-server \ "),
                      #       p(HTML('&emsp;'),"--bind shiny-server.conf:/etc/shiny-server/shiny-server.conf shiny.simg"),
                      #       p(HTML('&emsp;'),"---------------------------------------------------------------------------"),
                      #       p("For custom applications, also add --bind /srv/shiny-server:/srv/shiny-server"),
                      #       p(HTML('&emsp;'),"To see your applications, open your browser to http://127.0.0.1:9870 or"),
                      #       p(HTML('&emsp;'),"open a ssh connection from your computer to your cluster.")
                      #     )
                      # ),
                      # p("The configuration is generated in your present working directory:"),
                      # box(width=12,
                      #     p("$ cat shiny-server.conf"),
                      #     br(),
                      #     code(
                      #       p("run_as vanessa;"),
                      #       p("server {"),
                      #       p(HTML('&emsp;'),"listen 9098;"),
                      #       
                      #       p(HTML('&emsp;'),"# Define a location at the base URL"),
                      #       p(HTML('&emsp;'),"location / {"),
                      #       
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"# Host the directory of Shiny Apps stored in this directory"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"site_dir /srv/shiny-server;"),
                      #       
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"# Log all Shiny output to files in this directory"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"log_dir /tmp/shiny-server.PtVRXE;"),
                      #       
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"# When a user visits the base URL rather than a particular application,"),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"# an index of the applications available in this directory will be shown."),
                      #       p(HTML('&emsp;'),HTML('&emsp;'),"directory_index on;"),
                      #       p(HTML('&emsp;'),"}"),
                      #       p("}")
                      #     )
                      # ),
                      # p("You can also choose to disable the indexing, meaning that someone that navigates to the root of the server (at the port) won't be able to explore all of your apps."),
                      # code("$ /bin/bash prepare_template.sh --disable-index"),br(),br(),
                      # p("You can also customize the port, temporary folder, 'run_as' user, and base (if somewhere other than /srv/shiny-server)")
                      # ,
                      # br(),
                      # h4("3. Start server"),
                      # p("Once you have that template, follow the instructions to run the container. The temporary folder is already created for you."),
                      # box(width=12,
                      #     p("$ singularity run --bind /tmp/shiny-server.gG1X2Z/logs:/var/log/shiny \\ "),
                      #     p(HTML('&emsp;'),"--bind /tmp/shiny-server.gG1X2Z/lib:/var/lib/shiny-server \\ "),
                      #     p(HTML('&emsp;'),"--bind shiny-server.conf:/etc/shiny-server/shiny-server.conf shiny.simg"),
                      #     br(),
                      #     code(
                      #       p("[2018-04-07T00:14:17.403] [INFO] shiny-server - Shiny Server v1.5.7.890 (Node.js v8.10.0)"),
                      #       p("[2018-04-07T00:14:17.405] [INFO] shiny-server - Using config file '/etc/shiny-server/shiny-server.conf'"),
                      #       p("[2018-04-07T00:14:17.456] [INFO] shiny-server - Starting listener on 0.0.0.0:9870")
                      #     )
                      # ),br(),br(),
                      # h4("Custom application"),
                      # p("When you run the container, if you add a bind to a folder of your own apps in /srv/shiny-server/, you can add your custom applications. The bind would look something like:"),
                      # code("--bind /path/to/apps/folder:/srv/shiny-server")
                      
                      
                      # tags$div(tags$ul(
                      #   tags$li("test1"),
                      #   tags$li("test2"),
                      #   tags$li("test3")),  style = "font-size: 15px")
                    )
                  )
                )
  )
  
)
