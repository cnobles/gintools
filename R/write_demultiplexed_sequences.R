#' Write demultiplexed sequences to fasta/q files in an output director
#' 
#' \code{write_demultiplexed_sequences} will write fasta/q files to an output 
#' directory when supplied with grouping information (`multiplexed_data`).
#' 
#' @description This function will read in a sequence file, parse it given 
#' `multiplexed_data` information and write the parsed files to an output 
#' directory. Specify the `type` as "fasta" or "fastq" and give a sequence file
#' path (`read_file_path`) and the output folder path (`out_folder`). Data can 
#' either be comressed (into .gz format) or uncompressed. 
#' 
#' @param read_file_path path to sequence file (fasta/q) to be demultiplexed.
#' @param multiplexed_data
#' @param out_folder path to output directory.
#' @param type character identifying the read class (R1, I1, I2, R2).
#' @param read_name_pattern regular expression to be applied to the sequence 
#' names for identification. Default: "[\\w\\:\\-\\+]+".
#' @param compress logical indicating if the output files should be gzip 
#' compressed (.gz, default) or uncompressed.
#' 
#' @author Christopher Nobles, Ph.D.
#' 

write_demultiplexed_sequences <- function(read_file_path, multiplexed_data, 
                                          out_folder, type, 
                                          read_name_pattern = "[\\w\\:\\-\\+]+", 
                                          compress = TRUE){
  # Require packages for parallel processing
  dependencies <- c("stringr", "ShortRead", "Biostrings")
  loaded <- sapply(dependencies, require, character.only = TRUE)
  stopifnot(all(loaded))
  # Load read sequences and sequence names then write to file
  reads <- ShortRead::readFastq(read_file_path)
  reads@id <- Biostrings::BStringSet(
    stringr::str_extract(as.character(ShortRead::id(reads)), read_name_pattern))
  reads <- reads[multiplexed_data$index]
  reads <- split(reads, multiplexed_data$sampleName)
  null <- lapply(1:length(reads), function(i, reads, type, out_folder, compress){
    if(compress){  
      file_path <- file.path(
        out_folder, paste0(names(reads[i]), ".", type, ".fastq.gz"))
    }else{
      filePath <- file.path(
        out_folder, paste0(names(reads[i]), ".", type, ".fastq"))
    }
    if(file.exists(file_path)) unlink(file_path)
    ShortRead::writeFastq(reads[[i]], file = file_path, compress = compress)
    message(
      paste0("\nWrote ", length(reads[[i]]), " reads to:\n", file_path, "."))
  }, reads = reads, type = type, out_folder = out_folder, compress = compress)
  return(list(read_file_path, type, out_folder))
}