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
      span(class = "logo-lg", "Genome browser"), 
      img(src = "https://image.flaticon.com/icons/svg/268/268752.svg")),
    dropdownMenuOutput("messageMenu")
  ),
  
  
  ############
  ### SIDE ###
  ############
  
  dashboardSidebar(
    sidebarMenu(id = "sideBar",
                menuItem("Input data", tabName = "input", icon = icon("folder-open")),
                menuItemOutput("outGenome"),
                # hr(),uiOutput("seq_link"),
                # downloadButton('download',"Download the data",style="color: #fff; background-color: #006502; border-color: #2e6da4"),
                menuItemOutput("link"),hr(),
                menuItem("Download App", tabName = "singularity", icon = icon("box-open"))
    )
  ),
  
  
  ############
  ### BODY ###
  ############
  
  dashboardBody(useShinyjs(),
                tabItems(
                  tabItem(
                    tabName = "input",
                    fluidRow(
                      boxPlus(width=5,title = "Input Tab",closable = F,
                              status = "info", solidHeader = TRUE,
                              fluidRow(
                                column(width=6,boxPad(style = "border-right: 1px solid #BDBDBD",
                                                      checkboxGroupInput(inputId = "peaks",
                                                                         label = "Peaks calls format of the regions",
                                                                         choiceNames = c('Narrow Peaks','Broad Peaks'),
                                                                         choiceValues = c("narrow","broad"),
                                                                         selected='narrow'))),
                                column(width=6,
                                       radioGroupButtons(inputId="choice_track", label="Select gene or transcript track to visualize:", 
                                                         choices = c("gene","transcript"))
                                )),
                              selectizeInput('regElement',
                                             label = "Regulatory element:",
                                             choices = NULL,
                                             options = list(placeholder = 'select 1 or more',
                                                            maxItems = 1500,multiple = T, searchConjunction = 'and')),
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
                                                          choiceNames = c('Global consensus regions','Tissue consensus regions','Cell-lines consensus regions'),
                                                          choiceValues = c("global","tissue","cell-line"),
                                                          selected='global')),
                                column(width=6,
                                       conditionalPanel(condition = 'input.consensus.includes("tissue")',
                                                        selectizeInput('tissue', label = "Tissue:",
                                                                       choices = NULL,
                                                                       options = list(placeholder = 'select 1 or more',
                                                                                      maxItems = 1500,multiple = T, searchConjunction = 'and'))))),
                              checkboxGroupInput(inputId = "snps",
                                                 label = "Add SNPs track.",
                                                 choiceName = c('SNPs - dbSNP v151'),
                                                 choiceValues = c('snp')),
                              fluidRow(
                                column(width=6,
                                       checkboxGroupInput(inputId = "tba",
                                                          label = "Add TBA info to global and tissue tracks.",
                                                          choiceName = c('TBA - Total Binding Affinity'),
                                                          choiceValues = c('tba'))),
                                conditionalPanel(condition = "input.tba == 'tba'",
                                                 column(width=6,
                                                                checkboxGroupInput(inputId = "pvalues",label = "p-value:",
                                                                           choices = c("0.1","0.05","0.01","0.001","0.0001","0.00001"),
                                                                           selected = "0.00001"))
                                )
                              )
                              
                      ),
                      box(width=3,title = "Cell-lines Tree selection",
                          status = "info", solidHeader = TRUE,
                          conditionalPanel(condition = 'input.consensus.includes("cell-line")',
                                           shinyTree("tree", checkbox = TRUE,theme = 'proton',themeIcons = F,themeDots = F))
                      ),
                      box(width = 4,title = "Search Tab",
                          status = "warning", solidHeader = TRUE,
                          textInput("region", "chromosome region", "chr7:139,424,940-141,784,100"),
                          actionButton("Run", "Search by position"),hr(),
                          selectizeInput('genes', label = "Gene symbol:",
                                         choices = NULL, options = list(
                                           placeholder = 'Select a gene',
                                           maxItems = 1,multiple = F, searchConjunction = 'and')),
                          actionButton("searchGenes", "Search by gene"),
                          conditionalPanel(condition='input.menu == "hidden"',
                                           selectInput(
                                             inputId="selectData",
                                             label=" ", selected = NULL,
                                             choices=c( "title" )
                                           )
                          )
                          
                      )
                    )
                  ),
                  tabItem(
                    tabName = "genome",
                    fluidPage(
                      width = 12,
                      TnTOutput("outcomp",height = 'auto'),hr(),
                      tags$div(id = 'placeholder')
                    )
                  ),
                  tabItem(
                    tabName = "singularity",
                    fluidPage(
                      width = 12,
                      h1("Download singularity image"),
                      br(),
                      h2("Download"),
                      br(),
                      p("Download a singuarity image to run this shiny app on your local server."),
                      p("Download link:    ",downloadLink('downloadData', 'genomeBrowser_v1.tar')," 17 Oct 2019 (~30GB)"),
                      br(),
                      h2("Instructions"),
                      br(),
                      h4("1. Unpack TAR archive"),
                      box(width=12,"tar -xvf genomeBrowser.tar"),
                      br(),
                      h4("2. Prepare shiny-server configuration file"),
                      p("You will first need to generate a custom configuration for your user, and it will give you instructions for usage:"),
                      box(width=12,
                          p("$ /bin/bash prepare_conf.sh"),
                          br(),
                          code(
                            p("Steps:"),
                            p(HTML('&emsp;'),"----------------------------------------------------------------------"),
                            p(HTML('&emsp;'),"1. Use this script to prepare your shiny-server.conf (configuration)"),
                            
                            p(HTML('&emsp;'),HTML('&emsp;'),"/bin/bash prepare_template.sh"),
                            
                            p(HTML('&emsp;'),"----------------------------------------------------------------------"),
                            p(HTML('&emsp;'),"2. If needed, you can provide the following arguments"),
                            
                            p(HTML('&emsp;'),"Commands:"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"help: show help and exit"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"start: the generation of your config"),
                            
                            p(HTML('&emsp;'),"Options:"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"--port:  the port for the application (e.g., shiny default is 3737)"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"--user:  the user for the run_as directive in the shiny configuration"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"--base: base folder with applications"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"--logs: temporary folder with write for logs (not required)"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"--disable-index: disable directory indexing"),
                            
                            p(HTML('&emsp;'),"----------------------------------------------------------------------"),
                            p(HTML('&emsp;'),"3. Make sure Singularity is loaded, and run the container using "),
                            p(HTML('&emsp;'),HTML('&emsp;'),"the commands shown by the template.")
                          )
                      ),
                      p("When you add 'start' it will do the generation. Here we don't supply any arguments so that they are randomly generated."),
                      box(width=12,
                          p("$ /bin/bash prepare_template.sh start"),
                          br(),
                          code(
                            p("Generating shiny configuration..."),
                            p("port: 9870"),
                            p("logs: /tmp/shiny-server.gG1X2Z"),
                            p("base: /srv/shiny-server/shiny_genomeBrowser"),
                            p("Server logging will be in /tmp/shiny-server.gG1X2Z"),
                            
                            p("To run your server:"),
                            
                            p(HTML('&emsp;'),"module load singularity/2.4.6"),
                            p(HTML('&emsp;'),"singularity run --bind /tmp/shiny-server.gG1X2Z/logs:/var/log/shiny \ "),
                            p(HTML('&emsp;'),"--bind /tmp/shiny-server.gG1X2Z/lib:/var/lib/shiny-server \ "),
                            p(HTML('&emsp;'),"--bind shiny-server.conf:/etc/shiny-server/shiny-server.conf shiny.simg"),
                            p(HTML('&emsp;'),"---------------------------------------------------------------------------"),
                            p("For custom applications, also add --bind /srv/shiny-server:/srv/shiny-server"),
                            p(HTML('&emsp;'),"To see your applications, open your browser to http://127.0.0.1:9870 or"),
                            p(HTML('&emsp;'),"open a ssh connection from your computer to your cluster.")
                          )
                      ),
                      p("The configuration is generated in your present working directory:"),
                      box(width=12,
                          p("$ cat shiny-server.conf"),
                          br(),
                          code(
                            p("run_as vanessa;"),
                            p("server {"),
                            p(HTML('&emsp;'),"listen 9098;"),
                            
                            p(HTML('&emsp;'),"# Define a location at the base URL"),
                            p(HTML('&emsp;'),"location / {"),
                            
                            p(HTML('&emsp;'),HTML('&emsp;'),"# Host the directory of Shiny Apps stored in this directory"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"site_dir /srv/shiny-server;"),
                            
                            p(HTML('&emsp;'),HTML('&emsp;'),"# Log all Shiny output to files in this directory"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"log_dir /tmp/shiny-server.PtVRXE;"),
                            
                            p(HTML('&emsp;'),HTML('&emsp;'),"# When a user visits the base URL rather than a particular application,"),
                            p(HTML('&emsp;'),HTML('&emsp;'),"# an index of the applications available in this directory will be shown."),
                            p(HTML('&emsp;'),HTML('&emsp;'),"directory_index on;"),
                            p(HTML('&emsp;'),"}"),
                            p("}")
                          )
                      ),
                      p("You can also choose to disable the indexing, meaning that someone that navigates to the root of the server (at the port) won't be able to explore all of your apps."),
                      code("$ /bin/bash prepare_template.sh --disable-index"),br(),br(),
                      p("You can also customize the port, temporary folder, 'run_as' user, and base (if somewhere other than /srv/shiny-server)")
                      ,
                      br(),
                      h4("3. Start server"),
                      p("Once you have that template, follow the instructions to run the container. The temporary folder is already created for you."),
                      box(width=12,
                          p("$ singularity run --bind /tmp/shiny-server.gG1X2Z/logs:/var/log/shiny \\ "),
                          p(HTML('&emsp;'),"--bind /tmp/shiny-server.gG1X2Z/lib:/var/lib/shiny-server \\ "),
                          p(HTML('&emsp;'),"--bind shiny-server.conf:/etc/shiny-server/shiny-server.conf shiny.simg"),
                          br(),
                          code(
                            p("[2018-04-07T00:14:17.403] [INFO] shiny-server - Shiny Server v1.5.7.890 (Node.js v8.10.0)"),
                            p("[2018-04-07T00:14:17.405] [INFO] shiny-server - Using config file '/etc/shiny-server/shiny-server.conf'"),
                            p("[2018-04-07T00:14:17.456] [INFO] shiny-server - Starting listener on 0.0.0.0:9870")
                          )
                      ),br(),br(),
                      h4("Custom application"),
                      p("When you run the container, if you add a bind to a folder of your own apps in /srv/shiny-server/, you can add your custom applications. The bind would look something like:"),
                      code("--bind /path/to/apps/folder:/srv/shiny-server")
                      

                      # tags$div(tags$ul(
                      #   tags$li("test1"),
                      #   tags$li("test2"),
                      #   tags$li("test3")),  style = "font-size: 15px")
                    )
                  )
                )
  )
  
)
