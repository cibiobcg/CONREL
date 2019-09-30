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
                menuItemOutput("link")
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
                      conditionalPanel(condition='input.tooltipTable',
                                       boxPlus(width = 12,title = "Region info:",collapsible = T,closable = F,
                                               dataTableOutput("tTable")
                                       ))
                    )
                  )
                )
  )
  
)
