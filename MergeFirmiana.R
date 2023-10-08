#' @title MergeFirmiana



#' @function MergeFirmiana.pro
#' @description merge proteomic firmiana results
#' @param inpath  the path where firmiana result storaged
#' @param selectFileName select a part of files in the path to merge  
#' @param file_pattern select the pattern of file names in the path to merge 
#' @param columns the coloumns in each file to merge: Symbol Area iBAQ FoT(1e-6)
#' @param ID the gene or peptide ID
#' @param outpath where to save the result. the default is 0 for not saving the result
#' @import stringr
#' @import dplyr

MergeFirmiana.pro <- function(inpath, selectFileNames=c(), file_pattern = '^Exp(\\d{6})(.*)(\\.txt)$', 
                                    columns=c('Symbol','iBAQ'), ID = 'Symbol', outpath = 0, unique.pep = 0) {
  require(stringr)
  require(dplyr)
  #set dataframe
  if (length(selectFileNames)){
    FilePath <- sapply(selectFileNames, function(x) paste(inpath,paste(x,'.txt',sep=''),sep = '\\'))
  } else {
    File <-list.files(inpath, pattern = file_pattern) 
    selectFileNames <- sapply(File, function(x) str_split(x,'\\.',simplify = T)[1])
    FilePath <- sapply(File,function(x) paste(inpath,x,sep = '\\'))
  }
  FileFrame=data.frame('Name'= selectFileNames,'Path'=FilePath)
  #cycle
  for (i in seq(nrow(FileFrame))) {
    if (i == 1) {
      merge_data <- read.delim(FileFrame$Path[i], sep = '\t')
      merge_data <- merge_data %>% dplyr::filter(Unique.Peptide.Num >= unique.pep)
      merge_data <- as.data.frame(merge_data)[, columns]
    } else {
      single_data <- read.delim(FileFrame$Path[i], sep = '\t')
      single_data <- single_data %>% dplyr::filter(Unique.Peptide.Num >= unique.pep)
      single_data <- single_data[,columns]
      
      merge_data <- merge(merge_data, single_data, by = ID, all = T)
    }
    colnames(merge_data)[ncol(merge_data)] <- FileFrame$Name[i]
  }
  #output
  if(outpath!=0) {
    write.table(merge_data, outpath, na = "NA", row.names = F, quote = F, sep = '\t')
  }
  return(merge_data)
}

MergeFirmiana.pep <- function(inpath, selectFileNames=c(), file_pattern = '^Exp(\\d{6})(.*)(\\.txt)$', 
                              columns=c('Symbol','iBAQ'), ID = 'Symbol', f = NA, outpath = 0) {
  require(stringr)
  require(dplyr)
  pb <- txtProgressBar(style = 3)
  start_time <- Sys.time()
  #set dataframe
  if (length(selectFileNames)){
    FilePath <- sapply(selectFileNames, function(x) paste(inpath,paste(x,'.txt',sep=''),sep = '\\'))
  } else {
    File <-list.files(inpath, pattern = file_pattern) 
    selectFileNames <- sapply(File, function(x) str_split(x,'\\.',simplify = T)[1])
    FilePath <- sapply(File,function(x) paste(inpath,x,sep = '\\'))
  }
  FileFrame=data.frame('Name'= selectFileNames,'Path'=FilePath)
  #cycle
  for (i in seq(nrow(FileFrame))) {
    setTxtProgressBar(pb, i/nrow(FileFrame))
    if (i == 1) {
      merge_data <- read.csv(FileFrame$Path[i])
      merge_data <- merge_data[, columns]
      if (typeof(f) == 'closure') merge_data <- f(merge_data)
    } else {
      single_data <- read.csv(FileFrame$Path[i])
      single_data <- single_data[,columns]
      if (typeof(f) == 'closure') single_data <- f(single_data)
      merge_data <- merge(merge_data, single_data, by = ID, all = T)
    }
    colnames(merge_data)[ncol(merge_data)] <- FileFrame$Name[i]
  }
  close(pb)
  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
  #output
  if(outpath!=0) {
    write.table(merge_data, outpath, na = "NA", row.names = F, quote = F, sep = '\t')
  }
  return(merge_data)
}

# data <- read.csv("D:/Projects/Pan-sarcoma/Data/SS_PhoFirmianaResult_processed/result1/Exp114756_1.csv")
# data <- data[,c('Aligned.sequence','Abundance')]
MergeSameRow_pep <- function(data) {
  require(stringr)
  data <- data[-which(is.na(data$Aligned.sequence)),]
  
  dup_line_mul <- which(str_detect(data[,'Aligned.sequence'],'/')) #example ADEDJWEON/NCWOECEW
  dup_seq_mul <- unique(data[dup_line_mul,'Aligned.sequence'])
  if (length(dup_seq_mul) != 0) {
      for (i in dup_seq_mul) {
        index <- which(data[,'Aligned.sequence'] == i)
        data[index,'Aligned.sequence'] <- as.character(str_split(i, pattern = '/', simplify = T))
      }
    }
  
  dup_line <- which(duplicated(data[,'Aligned.sequence']))
  dup_seq_same <- unique(data[dup_line, 'Aligned.sequence'])
  if (length(dup_seq_same) != 0) {
    for (j in dup_seq_same) {
      index <- which(data[,'Aligned.sequence'] == j)
      data[index[1],'Abundance'] <-  sum(data[index,'Abundance'], na.rm = T)
    }
      data <- data[-dup_line,]
  }
  return(data)
}

# MergeFirmiana.pho <- function(phospho_files_dir,phospho_files_pattern = '(\\.txt)$',  selectFileName = c(), protein_sequence_fasta_path, outpath) {
#   # save the proteins in protein sequence database as a list
#   require(stringr)
#   fasta <- read.delim(protein_sequence_fasta_path, header = F)
#   protein_sequence_data <- list()
#   for (i in seq(nrow(fasta))) {
#     if (str_sub(fasta[i],1,1) == '>') gi_name <- str_sub(fasta[i], 2)
#     else protein_sequence_data[[gi_name]] <- fasta[i]
#   }
#   
#   # create a new dir to save the result
#   if (!dir.exists(outpath)) dir.create(outpath)
#   
#   # a cylce to match the peptides with the protein 
#   phospho_files_dir <- normalizePath(phospho_files_dir)
#   phospho_files <- list.files(phospho_files_dir,pattern = phospho_files_pattern)
#   if (length(selectFileName) != 0) pho_files <- selectFileName
#   for (i in phospho_files) {
#     phospho_file <- read.table(paste(phospho_files_dir,i,sep = '/'), sep = '\t',header = T)
#     for (i in seq(nrow(phospho_file))) {
#       
#     }
#   }
#   
# }
