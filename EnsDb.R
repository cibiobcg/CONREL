ip = as.data.frame(installed.packages()[,c(1,3:4)])
ip = ip[is.na(ip$Priority),1:2,drop=FALSE]
if("EnsDb.Hsapiens.v75"%in%ip$Package){
  # library(EnsDb.Hsapiens.v75)
  EnsDb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
  hg19 = fread(paste0(inputFolder,"data/hg19.chrom.bed"),data.table = F)
  hg_assembly="hg19"
  
  tagWorkInProgress = tagList(NULL)
} else if("EnsDb.Hsapiens.v86"%in%ip$Package){
  # library(EnsDb.Hsapiens.v86)
  EnsDb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  base_path = "/shares/CIBIO-Storage/BCGLAB/CONRELhg38/"
  encodeFolder = "/shares/CIBIO-Storage/BCGLAB/CONRELhg38/consensus/"
  inputFolder = "/shares/CIBIO-Storage/BCGLAB/CONRELhg38/"
  TBA_folder = "/shares/CIBIO-Storage/BCGLAB/CONRELhg38/"
  hg19 = fread(paste0(inputFolder,"data/hg38.chrom.bed"),data.table = F)
  hg_assembly="hg38"
  
  tagWorkInProgress = tagList(div(style="border:2px red dotted;",
                                  p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
                                    "ATTENTION! WORK IN PROGESS"),
                                  p(style="font-size:16px;font-weight:bold;text-align: justify; text-justify: inter-word;",
                                    "The GRCh38 version of CONREL is online but not fully functional. Some TBA annotations are still on processing, then the TBA informations will not show for those CREs.")))
  tagWorkInProgress = tagList(NULL)
}

# Ensembl data
seqlevelsStyle(EnsDb) <- "UCSC"

genesEns = as.data.table(genes(EnsDb)$symbol)
allEns = genes(EnsDb,columns="gene_biotype",return.type="DataFrame")$gene_biotype
allEns_tx = genes(EnsDb,columns="tx_biotype",return.type="DataFrame")$tx_biotype
levelsColor = unique(allEns)
levelsColor_tx = unique(allEns_tx)
#mapColor <- c(grDevices::rainbow(10),grDevices::rainbow(10,start = 0.2),grDevices::rainbow(length(levelsColor)-20,start = 0.4))
mapColor <- grDevices::rainbow(length(levelsColor))
c25 <- c(
  "#1c86ee", "#E31A1C", "#008b00", "#6A3D9A","#FF7F00", "#d4d4d4", "#ffd700", "#7ec0ee", "#FB9A99", "#90ee90", "#CAB2D6","#FDBF6F",
  "#b3b3b3", "#eee685", "#b03060", "#ff83fa", "#ff1493", "#0000ff", "#36648b", "#00ced1", "#00ff00", "#8b8b00", "#cdcd00", "#8b4500", "#a52a2a"
)
mapColor = c(c25,c25)[1:length(levelsColor)]
mapColor_tx = c(c25,c25)[1:length(levelsColor_tx)]

tabEns = table(allEns)[match(levelsColor,names(table(allEns)))]
tabEns_tx = table(allEns_tx)[match(levelsColor_tx,names(table(allEns_tx)))]
legend_idx = order(tabEns,decreasing = T)
legend_idx_tx = order(tabEns_tx,decreasing = T)

colnames(genesEns) = c("gene_symbol")
chrFilter = AnnotationFilterList(SeqNameFilter("chr1"),SeqNameFilter("chr2"),SeqNameFilter("chr3"),SeqNameFilter("chr4"),
                                 SeqNameFilter("chr5"),SeqNameFilter("chr6"),SeqNameFilter("chr7"),SeqNameFilter("chr8"),
                                 SeqNameFilter("chr9"),SeqNameFilter("chr10"),SeqNameFilter("chr11"),SeqNameFilter("chr12"),
                                 SeqNameFilter("chr13"),SeqNameFilter("chr14"),SeqNameFilter("chr15"),SeqNameFilter("chr16"),
                                 SeqNameFilter("chr17"),SeqNameFilter("chr18"),SeqNameFilter("chr19"),SeqNameFilter("chr20"),
                                 SeqNameFilter("chr21"),SeqNameFilter("chr22"),SeqNameFilter("chrX"),SeqNameFilter("chrY"),
                                 SeqNameFilter("chrM"),logicOp = rep(c("|"),24))