# Observer that listens to changes in tooltip (click on a regions, Shiny.onInputChange() in javascript)
observeEvent(input$tooltipTable, {
  removeUI(
    ## pass in appropriate div id
    selector = '#tbaTab'
  )
  removeUI(
    ## pass in appropriate div id
    selector = '#regionTab'
  )
  
  tooltip = input$tooltipTable
  df = data.frame(label = tooltip[names(tooltip)=="label"],
                  value = tooltip[names(tooltip)%in%c("value","value1")])
  idx_typeCRE = grep("^typeCRE",df$label)
  #datatable form the list of cell lines used to build the specific consensus
  df_nCells = df[idx_typeCRE,]
  
  if(nrow(df_nCells)>0){
    typeCRE = strsplit(as.character(df_nCells[,2]),"%")[[1]]
    # dt_cells = data.table("title"="Number of cell lines used to build the consensus:",
    #                       "count"=length(list_cells),
    #                       "dt"=data.table(list_cells))
    # dt_cells <- dt_cells[, list(list_cell=list(.SD)), by = list(title,count)]
    dt_nCells <- datatable({
      #configure datatable. Hide row number and cars columns [0,4] and enable details control on plus sign column[1]
      #turn rows into child rows and remove from parent
      if(typeCRE[3]=="global"){
        dt_map = map_Ncells_Global[map_Ncells_Global$peaks==typeCRE[1]&map_Ncells_Global$cre==typeCRE[2],]
        cbind('title' = 'Number of cell lines used to build the global consensus: ',dt_map,' ' = 'Click to expand and see all the cell line used')
      } else if (typeCRE[3]=="tissue") {
        dt_map = map_Ncells_Tissue_table[map_Ncells_Tissue_table$peaks==typeCRE[1]&map_Ncells_Tissue_table$cre==strsplit(typeCRE[2],"\\.")[[1]][2]&map_Ncells_Tissue_table$tissue==strsplit(typeCRE[2],"\\.")[[1]][1],]
        dt_map = dt_map[,-4]
        cbind('title' = 'Number of cell lines used to build the tissue consensus: ',dt_map,' ' = 'Click to expand and see all the cell line used')
      } else if (typeCRE[3]=="cell"){
        dt_map = map_Ncells_Cell_table[map_Ncells_Cell_table$peaks==typeCRE[1]&map_Ncells_Cell_table$cre==strsplit(typeCRE[2],"\\.")[[1]][3]&map_Ncells_Cell_table$cell==strsplit(typeCRE[2],"\\.")[[1]][1],]
        dt_map = dt_map[,-4]
        cbind('title' = 'Number of experiments used to build the cell line consensus: ',dt_map,' ' = 'Click to expand and see all the cell line used')
        # typeCRE[1] narrow
      }
      # cbind(' ' = '&oplus;', mtcars_dt)
    },
    escape = -2,
    options = list(paging = FALSE,
                   searching = FALSE,
                   lengthChange = FALSE,
                   dom = 'tBfrtp',
                   buttons = list('copy',
                                  list(extend='csv',
                                       filename = 'cellLineInfo'),
                                  list(extend='excel',
                                       filename = 'cellLineInfo'),
                                  list(extend='pdf',
                                       filename= 'cellLineInfo')),
                   columnDefs = list(
                     list(visible = FALSE, targets = c(0,2,3,5)),
                     list(orderable = FALSE, className = 'details-control', targets = 6)
                   )
    ),colnames = rep("", 6),
    callback = JS("
                  table.column(1).nodes().to$().css({cursor: 'pointer'});

                  // Format cars object into another table
                  var format = function(d) {
  if(d != null) { 
    var result = ('<table id=\"child_' + d[2] + '_' + d[3] + '_' + d[4] + '\">').replace('.','_') + '<thead><tr>'
    for (var col in d[5][0]) {
      result += '<th>' + col + '</th>'
    }
    result += '</tr></thead></table>'
    return result
  } else {
    return '';
  }
}

var format_datatable = function(d) {
  var dataset = [];
  for (i = 0; i <=  d[5].length-1; i++) {
    var datarow = $.map(d[5][i], function(value, index) {
      return [value];
    });
    dataset.push(datarow);
  }
                  var subtable = $(('table#child_' + d[2] + '_' + d[3] + '_' + d[4]).replace('.','_')).DataTable({
                  'data': dataset,
                  'autoWidth': true, 
                  'deferRender': true, 
                  'info': false, 
                  'lengthChange': false, 
                  'ordering': true, 
                  'paging': false, 
                  'scrollX': false, 
                  'scrollY': false, 
                  'searching': false 
                  });
                  };

                  table.on('click', 'td.details-control', function() {
                  var td = $(this), row = table.row(td.closest('tr'));
                  if (row.child.isShown()) {
                  row.child.hide();
                  td.html('&oplus; Click to expand and see all the cell line used');
                  } else {
                  row.child(format(row.data())).show();
                  td.html('&CircleMinus; Click to collapse');
                  format_datatable(row.data())
                  }
                  });")
    )
    df = df[-idx_typeCRE,]
  }
    
  if(length(input$tba)>0){
    idx_tba = grep("^TBA",df$label)
    df_tba = df[idx_tba,]
    
    if(nrow(df_tba)>0){
      res = lapply(rev(1:nrow(df_tba)),function(i) {
        if(df_tba[i,2]!="%"){
          data.frame(p.value = paste0("< ",formatC(as.numeric(gsub("TBA_pvalue.","",df_tba[i,1])), format = "e", digits = 2)),
                     TF_Symbol_and_Code = strsplit(strsplit(as.character(df_tba[i,2]),"%")[[1]][1],"/")[[1]],
                     ALL_1000GP = strsplit(strsplit(as.character(df_tba[i,2]),"%")[[1]][2],"/")[[1]])
        }
      })
      
      df_tba = do.call(rbind,res)
    } else {
      insertUI(
        selector = "#placeholder",
        where = "beforeBegin",
        ui = tags$div(
          boxPlus(width = 12,collapsible = T,closable = F,title="Gene info:",
                  dataTableOutput("tTable")),
          id="regionTab"
        )
      )
    }
    df_tba = datatable(df_tba,selection = 'none',filter = 'top',rownames = FALSE,extensions = 'Buttons',
                       options = list(paging = TRUE,
                                      searching = TRUE,
                                      dom = 'tBfrtip',
                                      buttons = list('copy',
                                                     list(extend='csv',
                                                          filename = 'tbaInfo'),
                                                     list(extend='excel',
                                                          filename = 'tbaInfo'),
                                                     list(extend='pdf',
                                                          filename= 'tbaInfo'))),
                       )
    tbaDF(df_tba)
    if(length(idx_tba)>0){
      if(length(input$tba)>0) {
        insertUI(
          selector = "#placeholder",
          where = "afterEnd",
          ui = tags$div(
            boxPlus(width = 6,title = "TBA info:",collapsible = T,closable = F,
                    dataTableOutput("tbaTable")),
            id="tbaTab"
          )
        )
        insertUI(
          selector = "#placeholder",
          where = "beforeBegin",
          ui = tags$div(
            boxPlus(width = 6,collapsible = T,closable = F,title="Region info:",
                    dataTableOutput("tTable"),hr(),
                    dataTableOutput("cellTable")),
            id="regionTab"
          )
        )
      }
      df = df[-idx_tba,]
    }
    df = datatable(df,selection = 'none',rownames = FALSE,colnames = c("",""),extensions = 'Buttons',
                   options = list(paging = FALSE,
                                  searching = FALSE,
                                  dom="tBfrtp",
                                  buttons = list('copy',
                                                 list(extend='csv',
                                                             filename = 'geneInfo'),
                                                 list(extend='excel',
                                                      filename = 'geneInfo'),
                                                 list(extend='pdf',
                                                      filename= 'geneInfo'))))
    tooltipDF(df) # set reactiveVal to new value
    if(nrow(df_nCells)>0){cellDF(dt_nCells)}else{cellDF(NULL)}
    # df <- df[, list(cars=list(.SD)), by = list(mpg,cyl)]
  } else {
    if("entrezid"%in%df$label) {
      title = "Gene info:"
    } else {
      title = "Region info:"
    }
    insertUI(
      selector = "#placeholder",
      where = "beforeBegin",
      ui = tags$div(
        boxPlus(width = 12,collapsible = T,closable = F,title=title,
                dataTableOutput("tTable"),hr(),
                dataTableOutput("cellTable")),
        id="regionTab"
      )
    )
    df = datatable(df,selection = 'none',rownames = FALSE,colnames = c("",""),extensions = 'Buttons',
                   options = list(paging = FALSE,
                                  searching = FALSE,
                                  dom="tBfrtp",
                                  buttons = list('copy',
                                                 list(extend='csv',
                                                      filename = 'geneInfo'),
                                                 list(extend='excel',
                                                      filename = 'geneInfo'),
                                                 list(extend='pdf',
                                                      filename= 'geneInfo'))))
    tooltipDF(df) # set reactiveVal to new value
    if(nrow(df_nCells)>0){cellDF(dt_nCells)}else{cellDF(NULL)}
  }
})

# # Observer that listens the click on search by position or genes and reset the ata.frame with tooltip info
# observeEvent(c(input$Run,input$searchGenes), {
#   tooltipDF(NULL) # set reactiveVal to new value
#   tbaDF(NULL)
# })
