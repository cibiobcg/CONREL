# Check the validity of the input region: how it is written and also that END > START position
validateRegion <- function(input,reset=F) {
  output=NULL
  if(!grepl(paste0("^chr([1-9]|1[0-9]|2[0-2]|[MXY]):[0-9,]*-[0-9,]*$"),input$region)){
    output = append(output,"- Please provide a region in valid format.\ne.g. chr7:139424940-141784100")
  } else{
    el = convertPosition(input$region,TRUE)
    chr = elements[1]
    start = as.numeric(el[2])
    end = as.numeric(el[3])
    if (end<=start) {
      output = append(output,"- Please provide a region in valid format.\nEND position need to be > than START position")
    } else if(!(data.table::between(start, 1, hg19[which(hg19$V1==chr),3]) | data.table::between(end, 1, hg19[which(hg19$V1==chr),2]))) {
      output = append(output,"- Please provide a region in valid format.\npositions provided are outside chromosome")
    } else {
      NULL
    }
  }
  if(!is.null(output)){
    output = append("ERROR in input coordinates!",output)
  }
  if(reset){clickPos <<- FALSE}
  return(output)
}

# Check if the input gene is a valid gene and exists into ensembl DB
validateGenes <- function(input,reset=F) {
  output=NULL
  if(input$genes == ""){
    output = append(output,"- Please provide a gene")
  } else if (input$genes %in% genes(EnsDb.Hsapiens.v75)$symbol) {
    NULL
  } else {
    output = append(output,"- Please select a gene from the list")
  }
  if(!is.null(output)){
    output = append("\nERROR in input gene!",output)
  }
  if(reset){clickGene <<- FALSE}
  return(output)
}

# Check the validity of general input:
# type of peak, regulatory element, consensus region, tissue (if selected), are essentials
validateInput <- function(input) {
  output=NULL
  if(length(input$peaks)==0){
    output = append(output,"- Please select at least one type of peaks")
  }
  if(length(input$regElement)==0) {
    output = append(output,"- Please select at least one regulatory element")
  }
  if(length(input$consensus)==0){
    output = append(output,"- Please select at least one type of consensus regions")
  }
  if("tissue" %in% input$consensus & length(input$tissue)==0){
    output = append(output,"- Please select at least one tissue, or deselect 'Tissue consensus regions'")
  }
  if("tba" %in% input$tba & length(input$pvalues)==0){
    output = append(output,"- Please select at least one level of statistical significance, or deselect 'TBA'")
  }
  if(!is.null(output)){
    output = append("\nAt least 1 error has occured in Input data!",output)
  }
  return(output)
}