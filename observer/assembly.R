if("EnsDb.Hsapiens.v75"%in%ip$Package){
  output$assembly <- renderMenu({
    menuItem(text="Switch to GRCh38/hg38", href = "./conrelhg38",icon = icon("dna"))
  })
} else if("EnsDb.Hsapiens.v86"%in%ip$Package){
  output$assembly <- renderMenu({
    menuItem("Switch to GRCh37/hg19", href = "./conrel",icon = icon("dna"))
  })
}