# Check the availability of tissue consensus for narrow and broad peak
tissueAvailability <- function(tissue,peaks,regElement) {
  r = lapply(peaks,function(p){
    regElement=gsub("active enhancer","active",regElement)
    listFiles = list.files(paste0(encodeFolder,p,"Peaks/bgz_tissueConsensus"),
                           pattern = paste(paste0("(",paste(tissue,collapse = "|"),")"),
                                           paste0("(",paste(regElement,collapse = "|"),")"),
                                           "union",
                                           sep = "."))
    listFiles = listFiles[-grep(".tbi",listFiles)]
    if(length(listFiles)>0){
      data.frame(p,listFiles)
    }
  })
  return(as.data.frame(do.call(rbind,r)))
}

# The function returns the name of the cell-line selected in the tree
extractCellline <- function(cellList){
  unlist(lapply(cellList,function(x){
    if(list.depth(x)==3) {
      names(x[[1]][[1]])
    }
  }))
}

# Check the availability of tissue consensus for narrow and broad peak
celllineAvailability <- function(cellList,peaks,regElement) {
  cellline = extractCellline(cellList)
  r = lapply(peaks,function(p){
    regElement=gsub("active enhancer","active",regElement)
    listFiles = list.files(paste0(encodeFolder,p,"Peaks/bgz_consensus"),
                           pattern = paste(paste0("(",paste(cellline,collapse = "|"),")"),
                                           "intersect",
                                           paste0("(",paste(regElement,collapse = "|"),")"),
                                           sep = "."))
    listFiles = listFiles[-grep(".tbi",listFiles)]
    if(length(listFiles)>0){
      data.frame(p,listFiles)
    }
  })
  return(as.data.frame(do.call(rbind,r)))
}

# The function returns the files (combination of, in a DF) for global consensus selected by user
mappingInput <- function(consensus,peaks) {
  consensus = gsub("active enhancer","active",consensus)
  combFolder = expand.grid(peaks,consensus)
  return(combFolder)
}


# TISSUES consensus for all ENCODE (the list in the bcglab DB)
mapping = fread(paste0(encodeFolder,"mapping_cellLines.csv"),data.table = F)
mapping$name = gsub("'","",mapping$name)
mapping$name = gsub(" ","_",mapping$name)
mapping$name = gsub("/","_",mapping$name)
mapping$tissue = gsub(" ","_",mapping$tissue)

tissue = sort(unique(mapping$tissue))
tissueOG = tissue[-grep("/",tissue)]

# CELL-LINES consensus (list of cell-lines available) for narrow and broad peaks
narrow_cell = fread(paste0(encodeFolder,"narrowPeaks/dataTableNoReplicates.tsv"),data.table = F)
narrow_cell$Biosample.term.name = gsub(" ","_",narrow_cell$Biosample.term.name)
narrow_cell$tissueMapped = mapping[match(narrow_cell$Biosample.term.name,mapping$name),"tissue"]
narrow_cell = unique(narrow_cell[,c("Biosample.term.name","tissueMapped")])
idx=grep("/",narrow_cell$tissueMapped)
split_multiTissue = strsplit(narrow_cell[idx,2],"/")
res=lapply(1:length(split_multiTissue),function(x){
  data.frame(Biosample.term.name=narrow_cell[idx[x],1],tissueMapped=split_multiTissue[[x]])
})
narrow_cell = rbind(narrow_cell[-idx,],do.call(rbind,res))

broad_cell = fread(paste0(encodeFolder,"broadPeaks/dataTableNoReplicates.tsv"),data.table = F)
broad_cell$Biosample.term.name = gsub(" ","_",broad_cell$Biosample.term.name)
broad_cell$tissueMapped = mapping[match(broad_cell$Biosample.term.name,mapping$name),"tissue"]
broad_cell = unique(broad_cell[,c("Biosample.term.name","tissueMapped")])
idx=grep("/",broad_cell$tissueMapped)
split_multiTissue = strsplit(broad_cell[idx,2],"/")
res=lapply(1:length(split_multiTissue),function(x){
  data.frame(Biosample.term.name=broad_cell[idx[x],1],tissueMapped=split_multiTissue[[x]])
})
broad_cell = rbind(broad_cell[-idx,],do.call(rbind,res))
