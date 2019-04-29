#### There are some additional functions we will used in this package

#' DT.replace.NA
#'
#' replace NA with 0 in data.table
#'
#' @importFrom data.table set
DT.replace.NA  <- function(x, r = 0) {
  for (j in seq_len(ncol(x))) data.table::set(x, which(is.na(x[[j]])), j, r)
}


#' @name usePackage
#' @title usePackage
#' loading all packages simultaneously
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dependencies = TRUE)
  require(p, character.only = TRUE)
}



#' @name useBiocPackage
#' @title useBiocPackage
#' @import BiocInstaller
#' @param p The packages list that you want to install 
#' loading all bioconductor packages simultaneously
useBiocPackage <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p, ask = FALSE)
  }
  require(p, character.only = TRUE)
}






#' Calculate PSI based on STAR output
#'
#' Function for calculate PSI based on splice junction reads. We used to build
#' a \code{SCASSet-class}, or you can use it directly.
#'
#' @importFrom  parallel mclapply
#' @import data.table
#' @param path The directory of STAR output *SJ.out.tab files, you shou put
#' all *SJ.out.tab files in same directory.
#' @param pattern The filename suffix of SJ.out.tab file.
#' @param minSJ The minimum mapped reads of every SJs
#' @param SSMinSJ Minimum junction reads of a same splice site
#' @param minSJs The minimum of row Sums of every SJs
#' @param minSamps The minimum of samples(cells) which SJ reads more than \code{minSJ} of every SJ.
#' @param uniqueMapOnly logical value indicates whether only the unique value
#' is counted for subsequent operations.
#' @param NT The number of cores in multithreaded parallel computing.
#' @return A list with 4 tables, count for SJ reads, sj.info for STAR
#' annotation of all SJ, start.psi for PSI of all same start SJ, and
#' end.psi for PSI of all same end SJ.
MakeFromSJs <- function(path,
                        pattern = "SJ.out.tab",
                        minSJ = 2,
                        minSJs = 100,
                        SSMinSJ = 5,
                        minSamps = 2,
                        uniqueMapOnly = TRUE,
                        NT = 1){
  if(length(list.files(path, pattern, full.names = TRUE)) == 0){
    stop(paste("There are no SJ.out.tab in", path))
  }

  if( !is.numeric(minSJ) ) {
    stop(paste0("Invalid minSJ", minSJ, "and minSJ must be a numeric value"))
  }
  if( !is.numeric(minSJs) ) {
    stop(paste0("Invalid minSJs", minSJs, "and minSJs must be a numeric value"))
  }
  if( !is.numeric(SSMinSJ) ) {
    stop(paste0("Invalid SSMinSJ", SSMinSJ, "and SSMinSJ must be a numeric value"))
  }
  if( !is.numeric(minSamps) ) {
    stop(paste0("Invalid minSamps", minSamps, "and minSamps must be a numeric value"))
  }

  if( !is.logical(uniqueMapOnly) ) {
    stop(paste0("Invalid uniqueMapOnly", uniqueMapOnly, "and uniqueMapOnly must be logical value"))
  }

  ## read SJ files
  filenames <- list.files(path, pattern, full.names=TRUE)
  message("Reading Splice Junction files...")
  datalist <- parallel::mclapply(filenames, FUN = function(x){
    tmp <- unique(data.table::fread(x, header=FALSE, select = 1:8))
    if(as.logical(uniqueMapOnly)) tmp$Reads = tmp$V7 else tmp$Reads = tmp$V7 + tmp$V8
    tmp <- tmp[Reads >= minSJ, ]
    data.table::setnames(tmp, "Reads", tail(strsplit(x, "/")[[1]], n = 1))
    return(tmp)}, mc.cores = NT)

  ## STAR annotation information
  info <- Reduce(function(x, y) {rbind(x, y, fill=TRUE)}, datalist)
  reads <- rowSums(info[, -c(1:8)], na.rm = TRUE)
  info <- info[, 1:6]
  colnames(info) <- c("seqname", "start", "end", "strand", "motif", "annotation")
  info$reads <- reads
  info <- info[, .(readsSum = sum(reads), N = .N), by = c("seqname", "start", "end", "strand", "motif", "annotation")]
  info <- info[readsSum >= minSJs & N >= minSamps, ]
  data.table::setkey(info, seqname, start, end, strand)

  ## raw splice junction reads
  message("Merging raw Splice Junction counts...")
  datalist <- lapply(datalist, function(x) {
    tmp <- x[, c(1:4, 9)]
    colnames(tmp)[1:4] <- c("seqname", "start", "end", "strand")
    data.table::setkey(tmp, seqname, start, end, strand)
    tmp[info[, .(seqname, start, end, strand)],]
  })
  mat.filter <- cbind(info[, .(seqname, start, end, strand)], Reduce(cbind, parallel::mclapply(datalist, function(x) x[, 5], mc.cores = NT)))
  DT.replace.NA(mat.filter)
  names(mat.filter)[-c(1:4)] <- rownames(SampleInfo)
  mat.filter[, strand:=as.character(strand)]
  mat.filter[strand == "0", strand:="*"]
  mat.filter[strand == "1", strand:="+"]
  mat.filter[strand == "2", strand:="-"]
  data.table::setkey(mat.filter, seqname, start, end, strand)

  ## STAR annotation information
  info[, strand:=as.character(strand)]
  info[strand == "0", strand:="*"]
  info[strand == "1", strand:="+"]
  info[strand == "2", strand:="-"]
  data.table::setkey(info, seqname, start, end, strand)
  # sj.info <- info[mat.filter[, 1:4], ]

  ## same start PSI
  message("Calculating same start PSI...")
  data.table::setkey(mat.filter, seqname, start, strand)
  start.sj <- mat.filter[, lapply(.SD, function(x) {if(as.numeric(sum(x, na.rm = T)) >= SSMinSJ) as.integer(x) else rep(NA_integer_, length(x))}), .SDcols = -c(1:4), by = list(seqname, start, strand)]
  start.psi <- start.sj[, lapply(.SD, function(x) x/sum(x, na.rm = T)), .SDcols = -c(1:3), by = list(seqname, start, strand)]
  data.table::setkey(start.psi, seqname, start, strand)
  start.psi$end <- mat.filter$end
  data.table::setcolorder(start.psi, c("seqname", "start", "end", "strand"))
  data.table::setkey(start.psi, seqname, start, end, strand)

  ## same end PSI
  message("Calculating same end PSI...")
  data.table::setkey(mat.filter, seqname, end, strand)
  end.sj <- mat.filter[, lapply(.SD, function(x) {if(as.numeric(sum(x, na.rm = T)) >= SSMinSJ) as.integer(x) else rep(NA_integer_, length(x))}), .SDcols = -c(1:4), by = list(seqname, end, strand)]
  end.psi <- end.sj[, lapply(.SD, function(x) x/sum(x, na.rm = T)), .SDcols = -c(1:3), by = list(seqname, end, strand)]
  data.table::setkey(end.psi, seqname, end, strand)
  end.psi$start <- mat.filter$start
  data.table::setcolorder(end.psi, c("seqname", "start", "end", "strand"))
  data.table::setkey(end.psi, seqname, start, end, strand)

  message(" Saving result...")
  data.table::setkey(mat.filter, seqname, start, end, strand)
  return(list(count = mat.filter, sj.info = info, start.psi = start.psi, end.psi = end.psi))
}


#' GTF file index and AS type classify
#'
#' Get all annotated alternative splicing sites from given GTF file.
#' We strongly recommend using the GTF file used when STAR alignment.
#'
#' @import data.table
#' @import GenomicFeatures
#' @import SummarizedExperiment
#' @import S4Vectors
#' @importFrom  stringr str_replace_all
#'
#' @param sjInfo data.frame or data.table of splice junction. Information such as
#' seqnames, start, end and strand are essential.
#' @param gtfFile The file.path of GTF file when you used in STAR allignment.
#' @param NT The number of cores in multithreaded parallel computing.
#' @export
GTFIndexAndAStype <- function(sjInfo, gtfFile, NT = 1){
  if(!is(sjInfo, "data.table")) setDT(sjInfo)
  if(!all(is.element(c("seqnames", "start", "end", "strand"), colnames(sjInfo))))
    stop(paste0("You must have seqnames, start, end and strand columns in your sjInfo file."))
  if(!is.integer(sjInfo$start))
    stop(paste0("The type of start must be integer."))
  if(!is.integer(sjInfo$end))
    stop(paste0("The type of end must be integer."))
  if(!file.exists(gtfFile))
    stop("The GTF file does not exist! And you should provide the GTF file you used in STAR alignment.")

  # GTF index ===============================
  # AF/LE -----------------------------------
  txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtfFile, format = "gtf")
  txintron <- GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE)
  txintron <- txintron[S4Vectors::elementNROWS(txintron) > 0]
  s_s <- as.character(unlist(unique(GenomicRanges::strand(txintron))))
  s_l <- S4Vectors::elementNROWS(txintron)
  txintron_gr <- unlist(txintron)
  my.rank <- function(s, l) ifelse(s == "+", return(1:l), return(-(l:1)))
  GenomicRanges::mcols(txintron_gr)$intron_rank <- do.call(c, lapply(seq_along(s_s), function(i) my.rank(s_s[i], s_l[i])))
  GenomicRanges::mcols(txintron_gr)$tx_id <- names(txintron_gr)
  txintron_tab_rank <- data.table::as.data.table(txintron_gr)

  first_intron_tab <- txintron_tab_rank[abs(intron_rank) == 1, ]
  first_intron_tab[, junc:=paste0(seqnames, ":", start, "-", end, ":", strand)]
  setkey(first_intron_tab, tx_id)

  last_intron_tab <- txintron_tab_rank[ , .SD[abs(intron_rank) == max(abs(intron_rank)), ], by = "tx_id"]
  last_intron_tab[, junc:=paste0(seqnames, ":", start, "-", end, ":", strand)]
  data.table::setkey(last_intron_tab, tx_id)

  # ASS -------------------------------------
  exons <- data.table::as.data.table(unique(GenomicFeatures::exons(txdb, columns = c("gene_id", "EXONNAME"))))

  # AS type identification ==================
  # Alternative SJ --------------------------
  same.start.sj <- sjInfo[, .SD[.N >= 2, ], by = list(seqnames, start, strand)]
  same.start.sj.id <- same.start.sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "start", "strand")][, V1]
  same.start.sj$ID <- rep(same.start.sj.id, same.start.sj[, .N, by = c("seqnames", "start", "strand")][, N])
  same.start.2.sj <- same.start.sj[, .SD[.N == 2, ], by = list(seqnames, start, strand)]
  same.end.sj <- sjInfo[, .SD[.N >= 2, ], by = list(seqnames, end, strand)]
  same.end.sj.id <- same.end.sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "end", "strand")][, V1]
  same.end.sj$ID <- rep(same.end.sj.id, same.end.sj[, .N, by = c("seqnames", "end", "strand")][, N])
  same.end.2.sj <- same.end.sj[, .SD[.N == 2, ], by = list(seqnames, end, strand)]
  same.start.2.sj2 <- same.start.2.sj[, .(end1 = min(end), end2 = max(end)), by = list(seqnames, start, strand)]
  same.end.2.sj2 <- same.end.2.sj[, .(start1 = max(start), start2 = min(start)), by = list(seqnames, end, strand)]
  data.table::setkey(same.start.2.sj2, seqnames, start, end1, end2, strand)
  data.table::setkey(same.end.2.sj2, seqnames, end, start2, start1, strand)


  # SE & MXE SJ -----------------------------

  match_nage <- same.start.2.sj2[strand=="-", same.end.2.sj2[strand!="+",], by = c("seqnames", "start", "end1", "end2", "strand")]
  match_posi <- same.start.2.sj2[strand=="+", same.end.2.sj2[strand!="-",], by = c("seqnames", "start", "end1", "end2", "strand")]
  start_end_match <- rbind(match_posi, match_nage)
  colnames(start_end_match)[colnames(start_end_match) == "seqnames"] <- c("seqnames", "seqnames1")
  start_end_match <- start_end_match[as.logical(start_end_match[, seqnames]==start_end_match[, seqnames1]), ]

  SE <- start_end_match[start==start2 & end==end2 & end1 < start1,]
  SE[, ID:=paste0(seqnames, ":", start, "-", end1, ":", strand, "@", seqnames, ":", start1, "-", end, ":", strand, "@", seqnames, ":", start, "-", end2, ":", strand)]
  SE_SJ1 <- SE[, .(seqnames, start, end1, strand, ID)]
  colnames(SE_SJ1) <- c("seqnames", "start", "end", "strand", "ID")
  SE_SJ2 <- SE[, .(seqnames, start1, end, strand, ID)]
  colnames(SE_SJ2) <- c("seqnames", "start", "end", "strand", "ID")
  SE_SJ3 <- SE[, .(seqnames, start, end2, strand, ID)]
  colnames(SE_SJ3) <- c("seqnames", "start", "end", "strand", "ID")
  SE_SJ <- rbind(SE_SJ1, SE_SJ2, SE_SJ3)
  SE_SJ[, seqnames:=as.character(seqnames)]
  SE_SJ[, start:=as.integer(start)]
  SE_SJ[, end:=as.integer(end)]
  SE_SJ[, strand:=as.character(strand)]
  SE_SJ[, ID:=as.character(ID)]
  data.table::setkey(SE_SJ, seqnames, start, end, strand)

  MXE <- start_end_match[start < end1 & end1 < start2 & start2 < end2 & end2 < start1 & start1 < end, ]
  MXE <- MXE[!MXE[, paste0(seqnames, ":", start2, "-", end2, ":", strand)] %in% sjInfo[, paste0(seqnames, ":", start, "-", end, ":", strand)], ]
  MXE[, ID:=paste0(seqnames, ":", start, "-", end1, ":", strand, "@", seqnames, ":", start2, "-", end, ":", strand, "@", seqnames, ":", start, "-", end2, ":", strand, "@", seqnames, ":", start1, "-", end, ":", strand)]
  MXE_SJ1 <- MXE[, .(seqnames, start, end1, strand, ID)]; colnames(MXE_SJ1) <- c("seqnames", "start", "end", "strand", "ID")
  MXE_SJ2 <- MXE[, .(seqnames, start2, end, strand, ID)]; colnames(MXE_SJ2) <- c("seqnames", "start", "end", "strand", "ID")
  MXE_SJ3 <- MXE[, .(seqnames, start, end2, strand, ID)]; colnames(MXE_SJ3) <- c("seqnames", "start", "end", "strand", "ID")
  MXE_SJ4 <- MXE[, .(seqnames, start1, end, strand, ID)]; colnames(MXE_SJ4) <- c("seqnames", "start", "end", "strand", "ID")
  MXE_SJ <- rbind(MXE_SJ1, MXE_SJ2, MXE_SJ3, MXE_SJ4)
  MXE_SJ[, seqnames:=as.character(seqnames)]
  MXE_SJ[, start:=as.integer(start)]
  MXE_SJ[, end:=as.integer(end)]
  MXE_SJ[, strand:=as.character(strand)]
  MXE_SJ[, ID:=as.character(ID)]
  data.table::setkey(MXE_SJ, seqnames, start, end, strand)

  #### AF/LE SJ
  ## AF/LE of same start
  tmp1 <- merge(same.start.sj, unique(first_intron_tab[strand == "-", c("seqnames", "start", "end", "strand", "junc")]), by = c("seqnames", "start", "end", "strand"), all.x = T)
  same_start_afe <- same.start.sj[!ID %in% tmp1[is.na(junc), ID], ]
  tmp1 <- merge(same.start.sj, unique(last_intron_tab[strand == "+", c("seqnames", "start", "end", "strand", "junc")]), by = c("seqnames", "start", "end", "strand"), all.x = T)
  same_start_ale <- same.start.sj[!ID %in% tmp1[is.na(junc), ID], ]
  tmp1 <- merge(same.end.sj, unique(first_intron_tab[strand == "+", c("seqnames", "start", "end", "strand", "junc")]), by = c("seqnames", "start", "end", "strand"), all.x = T)
  same_end_afe <- same.end.sj[!ID %in% tmp1[is.na(junc), ID], ]
  tmp1 <- merge(same.end.sj, unique(last_intron_tab[strand == "-", c("seqnames", "start", "end", "strand", "junc")]), by = c("seqnames", "start", "end", "strand"), all.x = T)
  same_end_ale <- same.end.sj[!ID %in% tmp1[is.na(junc), ID], ]
  AFE_SJ <- rbind(same_start_afe, same_end_afe)[, c("seqnames", "start", "end", "strand", "ID")]
  ALE_SJ <- rbind(same_start_ale, same_end_ale)[, c("seqnames", "start", "end", "strand", "ID")]

  AFE_SJ[, seqnames:=as.character(seqnames)]
  AFE_SJ[, start:=as.integer(start)]
  AFE_SJ[, end:=as.integer(end)]
  AFE_SJ[, strand:=as.character(strand)]
  AFE_SJ[, ID:=as.character(ID)]
  data.table::setkey(AFE_SJ, seqnames, start, end, strand)
  ALE_SJ[, seqnames:=as.character(seqnames)]
  ALE_SJ[, start:=as.integer(start)]
  ALE_SJ[, end:=as.integer(end)]
  ALE_SJ[, strand:=as.character(strand)]
  ALE_SJ[, ID:=as.character(ID)]
  data.table::setkey(ALE_SJ, seqnames, start, end, strand)

  #### ASS SJ

  same.start.2.sj2[, ID:=paste0(seqnames, ":", start, "-", end1, ":", strand, "@", seqnames, ":", start, "-", end2, ":", strand)]
  same.end.2.sj2[, ID:=paste0(seqnames, ":", start1, "-", end, ":", strand, "@", seqnames, ":", start2, "-", end, ":", strand)]
  same_start_junc2_gr <- GenomicRanges::GRanges(seqnames = as.character(same.start.2.sj2[, seqnames]),
                                                ranges = IRanges::IRanges(start = as.integer(same.start.2.sj2[, end1]) + 1, end = as.integer(same.start.2.sj2[, end2])),
                                                strand = as.character(same.start.2.sj2[, strand]),
                                                ID = same.start.2.sj2[, ID])
  same_end_junc2_gr <- GenomicRanges::GRanges(seqnames = as.character(same.end.2.sj2[, seqnames]),
                                              ranges = IRanges::IRanges(start = as.integer(same.end.2.sj2[, start2]), end = as.integer(same.end.2.sj2[, start1]) - 1),
                                              strand = as.character(same.end.2.sj2[, strand]),
                                              ID = same.end.2.sj2[, ID])
  tmp1 <- data.table::as.data.table(same_start_junc2_gr)
  tmp3 <- data.table::as.data.table(same_end_junc2_gr)

  same_start_ass <- merge(tmp1, exons, by = c("seqnames", "start", "strand"))[end.x < end.y, ]
  same_end_ass <- merge(tmp3, exons, by = c("seqnames", "end", "strand"))[start.x > start.y, ]
  A3SS_junc <- unique(c(same_start_ass[strand == "+", ID], same_end_ass[strand == "-", ID]))
  A5SS_junc <- unique(c(same_start_ass[strand == "-", ID], same_end_ass[strand == "+", ID]))

  A3SS_SJ <- lapply(strsplit(stringr::str_replace_all(A3SS_junc, "-(?=[[:digit:]])", ":"), "[@:]"), function(x){
    matrix(x, ncol = 4, byrow = T, dimnames = list(NULL, c("seqnames", "start", "end", "strand")))
  })
  A3SS_SJ <- data.table::as.data.table(do.call(rbind, A3SS_SJ))
  A3SS_SJ$ID <- rep(A3SS_junc, each = 2)

  A5SS_SJ <- lapply(strsplit(stringr::str_replace_all(A5SS_junc, "-(?=[[:digit:]])", ":"), "[@:]"), function(x){
    matrix(x, ncol = 4, byrow = T, dimnames = list(NULL, c("seqnames", "start", "end", "strand")))
  })
  A5SS_SJ <- data.table::as.data.table(do.call(rbind, A5SS_SJ))
  A5SS_SJ$ID <- rep(A5SS_junc, each = 2)

  A3SS_SJ[, seqnames:=as.character(seqnames)]
  A3SS_SJ[, start:=as.integer(start)]
  A3SS_SJ[, end:=as.integer(end)]
  A3SS_SJ[, strand:=as.character(strand)]
  A3SS_SJ[, ID:=as.character(ID)]
  data.table::setkey(A3SS_SJ, seqnames, start, end, strand)
  A5SS_SJ[, seqnames:=as.character(seqnames)]
  A5SS_SJ[, start:=as.integer(start)]
  A5SS_SJ[, end:=as.integer(end)]
  A5SS_SJ[, strand:=as.character(strand)]
  A5SS_SJ[, ID:=as.character(ID)]
  data.table::setkey(A5SS_SJ, seqnames, start, end, strand)

  SE_SJ[, AS:="SE"]
  MXE_SJ[, AS:="MXE"]
  AFE_SJ[, AS:="AFE"]
  ALE_SJ[, AS:="ALE"]
  A3SS_SJ[, AS:="A3SS"]
  A5SS_SJ[, AS:="A5SS"]
  AS_SJ <- rbind(SE_SJ, MXE_SJ, AFE_SJ, ALE_SJ, A3SS_SJ, A5SS_SJ)

  return(S4Vectors::DataFrame(merge.data.frame(sjInfo, AS_SJ, by = c("seqnames", "start", "end", "strand"))))
}


#' Find host gene of SJs based on GTF file.
#'
#' We used \code{GenomicFeatures} to find the host gene of given splice
#' junctions.
#'
#' @references GenomicFeatures
#' @title HostGeneofSJ
#' @rdname HostGeneofSJ
#'
#' @param sjInfo data.frame or data.table of splice junction. Information such as
#' seqnames, start, end and strand are essential.
#' @param gtfFile The file.path of GTF file when you used in STAR allignment.
#' @param NT The number of cores in multithreaded parallel computing.
#'
#' @importFrom rtracklayer readGFF
#' @import GenomicFeatures
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors nrow
#' @importFrom S4Vectors DataFrame
#'
#' @export

HostGeneofSJ <- function(sjInfo, gtfFile, NT = 1){

  if(!is(sjInfo, "data.table")) setDT(sjInfo)
  if(!all(is.element(c("seqnames", "start", "end", "strand"), colnames(sjInfo))))
    stop(paste0("You must have seqnames, start, end and strand columns in your sjInfo file."))
  if(!is.integer(sjInfo$start))
    stop(paste0("The type of start must be integer."))
  if(!is.integer(sjInfo$end))
    stop(paste0("The type of end must be integer."))
  if(!file.exists(gtfFile))
    stop("The GTF file does not exist! And you should provide the GTF file you used in STAR alignment.")

  gtf <- data.table::as.data.table(rtracklayer::readGFF(gtfFile))
  protein_coding_gene <- gtf[gene_biotype == "protein_coding" & type == "gene", .(seqid, start, end, strand, gene_id)]

  # gene id ---------------------------------
  txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtfFile, format = "gtf")
  intron <- GenomicFeatures::intronicParts(txdb, linked.to.single.gene.only = FALSE)
  intron_tab <- data.table::as.data.table(intron)
  intron_tab[, names:=as.character(intron)]

  # exactly match ----
  junc_tu <- merge(sjInfo, intron_tab[, .(seqnames, start, end, strand, gene_id)], by = c("seqnames", "start", "end", "strand"), all.x = TRUE)
  junc_tu$gene_id <- parallel::mclapply(junc_tu$gene_id, FUN = function(x) {
    if(length(x) <= 1) {
      x
    } else {
      x[x %in% protein_coding_gene[, gene_id]]
    }
  }, mc.cores = NT)
  sj_gene_gf <- junc_tu[mapply(FUN = function(x) length(x) == 1, junc_tu$gene_id), ]
  sj_gene_gf$gene_id <- unlist(sj_gene_gf$gene_id)

  # mapRangesToIds ----
  tmp <- junc_tu[mapply(FUN = function(x) length(x) == 0, junc_tu$gene_id), ]

  parallel::mclapply(seq_len(nrow(tmp)), function(i){
    return(GenomicRanges::GRanges(seqnames = as.character(tmp[i, seqnames]),
                                  ranges = IRanges::IRanges(start = as.integer(tmp[i, start]), end = as.integer(tmp[i, end])),
                                  strand = as.character(tmp[i, strand])))
  }, mc.cores = NT) -> res_s

  grl_s <- GenomicRanges::GRangesList(res_s)
  names(grl_s) <- tmp[, paste0(seqnames, ":", start, "-", end, ":", strand)]
  res_s_gene <- GenomicFeatures::mapRangesToIds(txdb, grl_s, "gene")
  res_s_gene_cd <- parallel::mclapply(res_s_gene, FUN = function(x) {
    if(S4Vectors::nrow(x) == 1) {
      x
    } else {
      y = x$gene_id[x$gene_id %in% protein_coding_gene[, gene_id]]
      S4Vectors::DataFrame(gene_id = y)
    }
  }, mc.cores = NT)
  res_s_gene_cd <- res_s_gene_cd[mapply(S4Vectors::nrow, res_s_gene_cd) == 1]
  sj_gene_rang <- data.table::as.data.table(data.frame(seqnames = as.character(mapply(function(x) x[1], strsplit(names(res_s_gene_cd), ":"))),
                                                       start = as.integer(mapply(function(x) x[2], strsplit(names(res_s_gene_cd), "[:-]"))),
                                                       end = as.integer(mapply(function(x) x[3], strsplit(names(res_s_gene_cd), "[:-]"))),
                                                       strand = as.character(mapply(function(x) x[3], strsplit(names(res_s_gene_cd), ":"))),
                                                       gene_id = do.call(rbind, res_s_gene_cd)$gene_id))
  sj_gene_rang <- merge(sjInfo, sj_gene_rang, by = c("seqnames", "start", "end", "strand"))

  sj_host_gene <- rbind(sj_gene_rang, sj_gene_gf)

  # gene_name -------------------------------
  sj_host_gene <- merge(sj_host_gene, unique(gtf[, .(gene_id, gene_name, gene_biotype)]), by = "gene_id")
  HostGene <- merge(sjInfo, sj_host_gene, by = c("seqnames", "start", "end", "strand", "width", "motif", "annotation"), all.x = T)
  # return host gene ------------------------
  return(HostGene)
}



#' AS type and host gene of given SJs.
#'
#' We used \code{GenomicFeatures} to find the host gene of given splice
#' junctions.
#'
#' @references GenomicFeatures
#' @title AStypeAndHostGene
#' @rdname AStypeAndHostGene
#'
#' @param sjInfo data.frame or data.table of splice junction. Information such as
#' seqnames, start, end and strand are essential.
#' @param gtfFile The file.path of GTF file when you used in STAR allignment.
#' @param NT The number of cores in multithreaded parallel computing.
#'
#' @importFrom rtracklayer readGFF
#' @import GenomicFeatures
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors nrow
#' @importFrom S4Vectors DataFrame

AStypeAndHostGene <- function(sjInfo, gtfFile, NT = 1){
  if(!is(sjInfo, "data.table")) setDT(sjInfo)
  if(!all(is.element(c("seqnames", "start", "end", "strand"), colnames(sjInfo))))
    stop(paste0("You must have seqnames, start, end and strand columns in your sjInfo file."))
  if(!is.integer(sjInfo$start))
    stop(paste0("The type of start must be integer."))
  if(!is.integer(sjInfo$end))
    stop(paste0("The type of end must be integer."))
  if(!file.exists(gtfFile))
    stop("The GTF file does not exist! And you should provide the GTF file you used in STAR alignment.")

  message("Indexing GTF file...")
  # GTF index ===============================
  # AF/LE -----------------------------------
  suppressMessages(txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtfFile, format = "gtf"))
  txintron <- GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE)
  txintron <- txintron[S4Vectors::elementNROWS(txintron) > 0]
  s_s <- as.character(unlist(unique(GenomicRanges::strand(txintron))))
  s_l <- S4Vectors::elementNROWS(txintron)
  txintron_gr <- unlist(txintron)
  my.rank <- function(s, l) ifelse(s == "+", return(1:l), return(-(l:1)))
  GenomicRanges::mcols(txintron_gr)$intron_rank <- do.call(c, lapply(seq_along(s_s), function(i) my.rank(s_s[i], s_l[i])))
  GenomicRanges::mcols(txintron_gr)$tx_id <- names(txintron_gr)
  txintron_tab_rank <- data.table::as.data.table(txintron_gr)

  first_intron_tab <- txintron_tab_rank[abs(intron_rank) == 1, ]
  first_intron_tab[, junc:=paste0(seqnames, ":", start, "-", end, ":", strand)]
  setkey(first_intron_tab, tx_id)

  last_intron_tab <- txintron_tab_rank[ , .SD[abs(intron_rank) == max(abs(intron_rank)), ], by = "tx_id"]
  last_intron_tab[, junc:=paste0(seqnames, ":", start, "-", end, ":", strand)]
  data.table::setkey(last_intron_tab, tx_id)

  # ASS -------------------------------------
  exons <- data.table::as.data.table(unique(GenomicFeatures::exons(txdb, columns = c("gene_id", "EXONNAME"))))

  # AS type identification ==================
  message("Identifying AS type...")
  # Alternative SJ --------------------------
  same.start.sj <- sjInfo[, .SD[.N >= 2, ], by = list(seqnames, start, strand)]
  same.start.sj.id <- same.start.sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "start", "strand")][, V1]
  same.start.sj$ID <- rep(same.start.sj.id, same.start.sj[, .N, by = c("seqnames", "start", "strand")][, N])
  same.start.2.sj <- same.start.sj[, .SD[.N == 2, ], by = list(seqnames, start, strand)]
  same.end.sj <- sjInfo[, .SD[.N >= 2, ], by = list(seqnames, end, strand)]
  same.end.sj.id <- same.end.sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "end", "strand")][, V1]
  same.end.sj$ID <- rep(same.end.sj.id, same.end.sj[, .N, by = c("seqnames", "end", "strand")][, N])
  same.end.2.sj <- same.end.sj[, .SD[.N == 2, ], by = list(seqnames, end, strand)]
  same.start.2.sj2 <- same.start.2.sj[, .(end1 = min(end), end2 = max(end)), by = list(seqnames, start, strand)]
  same.end.2.sj2 <- same.end.2.sj[, .(start1 = max(start), start2 = min(start)), by = list(seqnames, end, strand)]
  data.table::setkey(same.start.2.sj2, seqnames, start, end1, end2, strand)
  data.table::setkey(same.end.2.sj2, seqnames, end, start2, start1, strand)


  # SE & MXE SJ -----------------------------

  match_nage <- same.start.2.sj2[strand=="-", same.end.2.sj2[strand!="+",], by = c("seqnames", "start", "end1", "end2", "strand")]
  match_posi <- same.start.2.sj2[strand=="+", same.end.2.sj2[strand!="-",], by = c("seqnames", "start", "end1", "end2", "strand")]
  start_end_match <- rbind(match_posi, match_nage)
  colnames(start_end_match)[colnames(start_end_match) == "seqnames"] <- c("seqnames", "seqnames1")
  start_end_match <- start_end_match[as.logical(start_end_match[, seqnames]==start_end_match[, seqnames1]), ]

  SE <- start_end_match[start==start2 & end==end2 & end1 < start1,]
  SE[, ID:=paste0(seqnames, ":", start, "-", end1, ":", strand, "@", seqnames, ":", start1, "-", end, ":", strand, "@", seqnames, ":", start, "-", end2, ":", strand)]
  SE_SJ1 <- SE[, .(seqnames, start, end1, strand, ID)]
  colnames(SE_SJ1) <- c("seqnames", "start", "end", "strand", "ID")
  SE_SJ2 <- SE[, .(seqnames, start1, end, strand, ID)]
  colnames(SE_SJ2) <- c("seqnames", "start", "end", "strand", "ID")
  SE_SJ3 <- SE[, .(seqnames, start, end2, strand, ID)]
  colnames(SE_SJ3) <- c("seqnames", "start", "end", "strand", "ID")
  SE_SJ <- rbind(SE_SJ1, SE_SJ2, SE_SJ3)
  SE_SJ[, seqnames:=as.character(seqnames)]
  SE_SJ[, start:=as.integer(start)]
  SE_SJ[, end:=as.integer(end)]
  SE_SJ[, strand:=as.character(strand)]
  SE_SJ[, ID:=as.character(ID)]
  data.table::setkey(SE_SJ, seqnames, start, end, strand)

  MXE <- start_end_match[start < end1 & end1 < start2 & start2 < end2 & end2 < start1 & start1 < end, ]
  MXE <- MXE[!MXE[, paste0(seqnames, ":", start2, "-", end2, ":", strand)] %in% sjInfo[, paste0(seqnames, ":", start, "-", end, ":", strand)], ]
  MXE[, ID:=paste0(seqnames, ":", start, "-", end1, ":", strand, "@", seqnames, ":", start2, "-", end, ":", strand, "@", seqnames, ":", start, "-", end2, ":", strand, "@", seqnames, ":", start1, "-", end, ":", strand)]
  MXE_SJ1 <- MXE[, .(seqnames, start, end1, strand, ID)]; colnames(MXE_SJ1) <- c("seqnames", "start", "end", "strand", "ID")
  MXE_SJ2 <- MXE[, .(seqnames, start2, end, strand, ID)]; colnames(MXE_SJ2) <- c("seqnames", "start", "end", "strand", "ID")
  MXE_SJ3 <- MXE[, .(seqnames, start, end2, strand, ID)]; colnames(MXE_SJ3) <- c("seqnames", "start", "end", "strand", "ID")
  MXE_SJ4 <- MXE[, .(seqnames, start1, end, strand, ID)]; colnames(MXE_SJ4) <- c("seqnames", "start", "end", "strand", "ID")
  MXE_SJ <- rbind(MXE_SJ1, MXE_SJ2, MXE_SJ3, MXE_SJ4)
  MXE_SJ[, seqnames:=as.character(seqnames)]
  MXE_SJ[, start:=as.integer(start)]
  MXE_SJ[, end:=as.integer(end)]
  MXE_SJ[, strand:=as.character(strand)]
  MXE_SJ[, ID:=as.character(ID)]
  data.table::setkey(MXE_SJ, seqnames, start, end, strand)

  #### AF/LE SJ
  ## AF/LE of same start
  tmp1 <- merge(same.start.sj, unique(first_intron_tab[strand == "-", c("seqnames", "start", "end", "strand", "junc")]), by = c("seqnames", "start", "end", "strand"), all.x = T)
  same_start_afe <- same.start.sj[!ID %in% tmp1[is.na(junc), ID], ]
  tmp1 <- merge(same.start.sj, unique(last_intron_tab[strand == "+", c("seqnames", "start", "end", "strand", "junc")]), by = c("seqnames", "start", "end", "strand"), all.x = T)
  same_start_ale <- same.start.sj[!ID %in% tmp1[is.na(junc), ID], ]
  tmp1 <- merge(same.end.sj, unique(first_intron_tab[strand == "+", c("seqnames", "start", "end", "strand", "junc")]), by = c("seqnames", "start", "end", "strand"), all.x = T)
  same_end_afe <- same.end.sj[!ID %in% tmp1[is.na(junc), ID], ]
  tmp1 <- merge(same.end.sj, unique(last_intron_tab[strand == "-", c("seqnames", "start", "end", "strand", "junc")]), by = c("seqnames", "start", "end", "strand"), all.x = T)
  same_end_ale <- same.end.sj[!ID %in% tmp1[is.na(junc), ID], ]
  AFE_SJ <- rbind(same_start_afe, same_end_afe)[, c("seqnames", "start", "end", "strand", "ID")]
  ALE_SJ <- rbind(same_start_ale, same_end_ale)[, c("seqnames", "start", "end", "strand", "ID")]

  AFE_SJ[, seqnames:=as.character(seqnames)]
  AFE_SJ[, start:=as.integer(start)]
  AFE_SJ[, end:=as.integer(end)]
  AFE_SJ[, strand:=as.character(strand)]
  AFE_SJ[, ID:=as.character(ID)]
  data.table::setkey(AFE_SJ, seqnames, start, end, strand)
  ALE_SJ[, seqnames:=as.character(seqnames)]
  ALE_SJ[, start:=as.integer(start)]
  ALE_SJ[, end:=as.integer(end)]
  ALE_SJ[, strand:=as.character(strand)]
  ALE_SJ[, ID:=as.character(ID)]
  data.table::setkey(ALE_SJ, seqnames, start, end, strand)

  #### ASS SJ

  same.start.2.sj2[, ID:=paste0(seqnames, ":", start, "-", end1, ":", strand, "@", seqnames, ":", start, "-", end2, ":", strand)]
  same.end.2.sj2[, ID:=paste0(seqnames, ":", start1, "-", end, ":", strand, "@", seqnames, ":", start2, "-", end, ":", strand)]
  same_start_junc2_gr <- GenomicRanges::GRanges(seqnames = as.character(same.start.2.sj2[, seqnames]),
                                                ranges = IRanges::IRanges(start = as.integer(same.start.2.sj2[, end1]) + 1, end = as.integer(same.start.2.sj2[, end2])),
                                                strand = as.character(same.start.2.sj2[, strand]),
                                                ID = same.start.2.sj2[, ID])
  same_end_junc2_gr <- GenomicRanges::GRanges(seqnames = as.character(same.end.2.sj2[, seqnames]),
                                              ranges = IRanges::IRanges(start = as.integer(same.end.2.sj2[, start2]), end = as.integer(same.end.2.sj2[, start1]) - 1),
                                              strand = as.character(same.end.2.sj2[, strand]),
                                              ID = same.end.2.sj2[, ID])
  tmp1 <- data.table::as.data.table(same_start_junc2_gr)
  tmp3 <- data.table::as.data.table(same_end_junc2_gr)

  same_start_ass <- merge(tmp1, exons, by = c("seqnames", "start", "strand"))[end.x < end.y, ]
  same_end_ass <- merge(tmp3, exons, by = c("seqnames", "end", "strand"))[start.x > start.y, ]
  A3SS_junc <- unique(c(same_start_ass[strand == "+", ID], same_end_ass[strand == "-", ID]))
  A5SS_junc <- unique(c(same_start_ass[strand == "-", ID], same_end_ass[strand == "+", ID]))

  A3SS_SJ <- lapply(strsplit(stringr::str_replace_all(A3SS_junc, "-(?=[[:digit:]])", ":"), "[@:]"), function(x){
    matrix(x, ncol = 4, byrow = T, dimnames = list(NULL, c("seqnames", "start", "end", "strand")))
  })
  A3SS_SJ <- data.table::as.data.table(do.call(rbind, A3SS_SJ))
  A3SS_SJ$ID <- rep(A3SS_junc, each = 2)

  A5SS_SJ <- lapply(strsplit(stringr::str_replace_all(A5SS_junc, "-(?=[[:digit:]])", ":"), "[@:]"), function(x){
    matrix(x, ncol = 4, byrow = T, dimnames = list(NULL, c("seqnames", "start", "end", "strand")))
  })
  A5SS_SJ <- data.table::as.data.table(do.call(rbind, A5SS_SJ))
  A5SS_SJ$ID <- rep(A5SS_junc, each = 2)

  A3SS_SJ[, seqnames:=as.character(seqnames)]
  A3SS_SJ[, start:=as.integer(start)]
  A3SS_SJ[, end:=as.integer(end)]
  A3SS_SJ[, strand:=as.character(strand)]
  A3SS_SJ[, ID:=as.character(ID)]
  data.table::setkey(A3SS_SJ, seqnames, start, end, strand)
  A5SS_SJ[, seqnames:=as.character(seqnames)]
  A5SS_SJ[, start:=as.integer(start)]
  A5SS_SJ[, end:=as.integer(end)]
  A5SS_SJ[, strand:=as.character(strand)]
  A5SS_SJ[, ID:=as.character(ID)]
  data.table::setkey(A5SS_SJ, seqnames, start, end, strand)

  SE_SJ[, AS:="SE"]
  MXE_SJ[, AS:="MXE"]
  AFE_SJ[, AS:="AFE"]
  ALE_SJ[, AS:="ALE"]
  A3SS_SJ[, AS:="A3SS"]
  A5SS_SJ[, AS:="A5SS"]
  AS_SJ <- rbind(SE_SJ, MXE_SJ, AFE_SJ, ALE_SJ, A3SS_SJ, A5SS_SJ)

  AS_Type <- S4Vectors::DataFrame(merge.data.frame(sjInfo, AS_SJ, by = c("seqnames", "start", "end", "strand")))

  #### Host gene ======================================
  message("Identifying host gene...")
  gtf <- data.table::as.data.table(rtracklayer::readGFF(gtfFile))
  protein_coding_gene <- gtf[gene_biotype == "protein_coding" & type == "gene", .(seqid, start, end, strand, gene_id)]

  # gene id ---------------------------------
  #   txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtfFile, format = "gtf")
  intron <- GenomicFeatures::intronicParts(txdb, linked.to.single.gene.only = FALSE)
  intron_tab <- data.table::as.data.table(intron)
  intron_tab[, names:=as.character(intron)]

  # exactly match ----
  junc_tu <- merge(sjInfo, intron_tab[, .(seqnames, start, end, strand, gene_id)], by = c("seqnames", "start", "end", "strand"), all.x = TRUE)
  junc_tu$gene_id <- parallel::mclapply(junc_tu$gene_id, FUN = function(x) {
    if(length(x) <= 1) {
      x
    } else {
      x[x %in% protein_coding_gene[, gene_id]]
    }
  }, mc.cores = NT)
  sj_gene_gf <- junc_tu[mapply(FUN = function(x) length(x) == 1, junc_tu$gene_id), ]
  sj_gene_gf$gene_id <- unlist(sj_gene_gf$gene_id)

  # mapRangesToIds ----
  tmp <- junc_tu[mapply(FUN = function(x) length(x) == 0, junc_tu$gene_id), ]

  parallel::mclapply(seq_len(nrow(tmp)), function(i){
    return(GenomicRanges::GRanges(seqnames = as.character(tmp[i, seqnames]),
                                  ranges = IRanges::IRanges(start = as.integer(tmp[i, start]), end = as.integer(tmp[i, end])),
                                  strand = as.character(tmp[i, strand])))
  }, mc.cores = NT) -> res_s

  grl_s <- GenomicRanges::GRangesList(res_s)
  names(grl_s) <- tmp[, paste0(seqnames, ":", start, "-", end, ":", strand)]
  res_s_gene <- GenomicFeatures::mapRangesToIds(txdb, grl_s, "gene")
  res_s_gene_cd <- parallel::mclapply(res_s_gene, FUN = function(x) {
    if(S4Vectors::nrow(x) == 1) {
      x
    } else {
      y = x$gene_id[x$gene_id %in% protein_coding_gene[, gene_id]]
      S4Vectors::DataFrame(gene_id = y)
    }
  }, mc.cores = NT)
  res_s_gene_cd <- res_s_gene_cd[mapply(S4Vectors::nrow, res_s_gene_cd) == 1]
  sj_gene_rang <- data.table::as.data.table(data.frame(seqnames = as.character(mapply(function(x) x[1], strsplit(names(res_s_gene_cd), ":"))),
                                                       start = as.integer(mapply(function(x) x[2], strsplit(names(res_s_gene_cd), "[:-]"))),
                                                       end = as.integer(mapply(function(x) x[3], strsplit(names(res_s_gene_cd), "[:-]"))),
                                                       strand = as.character(mapply(function(x) x[3], strsplit(names(res_s_gene_cd), ":"))),
                                                       gene_id = do.call(rbind, res_s_gene_cd)$gene_id))
  sj_gene_rang <- merge(sjInfo, sj_gene_rang, by = c("seqnames", "start", "end", "strand"))

  sj_host_gene <- rbind(sj_gene_rang, sj_gene_gf)

  # gene_name -------------------------------
  sj_host_gene <- merge(sj_host_gene, unique(gtf[, .(gene_id, gene_name, gene_biotype)]), by = "gene_id")
  HostGene <- merge(sjInfo, sj_host_gene, by = c("seqnames", "start", "end", "strand", "width", "motif", "annotation"), all.x = T)
  # return AS_Type & SJ host gene -----------

  return(list(AS_Type = AS_Type, HostGene = S4Vectors::DataFrame(HostGene)))
}



##obtain the junc.as from junction count table
##example junc.as = ASSite(junc,NT=3L)
ASSite <- function(junc, NT = 10){
  
  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")
  usePackage("plyr")
  
  # Identify alternative splicing events ----
  juncs <- junc$junctions
  junc.names <- data.table(chr = as.character(mapply(function(x) x[1], strsplit(juncs, split = ":"))),
                           start = as.integer(mapply(function(x) x[2], strsplit(juncs, split = "[:-]"))),
                           end = as.integer(mapply(function(x) x[3], strsplit(juncs, split = "[:-]"))),
                           strand = as.character(mapply(function(x) x[3], strsplit(juncs, split = ":"))))
  #rownames(junc.names) <- juncs
  junc.names$names <- juncs
  junc.names <- junc.names[order(chr, start, end),][chr%in%c(1:22,'X','Y')]
  #### add the exon-intron boundary site (different for + and - strand)
  #### for - strand (end_base)
  junc.names[, first_boundary:=paste0(chr,':',end,':-')]
  ### using boundary:= failed !
  ##### for + strand (start base)
  junc.names[strand=='+'] = transform(junc.names[strand=='+'], first_boundary=paste0(chr,':',start,':+'))
  junc.names <- as.data.frame(junc.names)
  
  junc.as <- do.call(c, mclapply(unique(junc.names$chr), function(chr) {
    junc.chr <- junc.names[junc.names$chr == chr, ]
    #### notice: when applu the dlply to data.table, will return eorror results !!!!
    same.start <- plyr::dlply(junc.chr, c("start"), function(x) x)
    same.start <- lapply(same.start,function(x){ x$share = 'same_start'; return(x)})
    same.end <- plyr::dlply(junc.chr, c("end"), function(x) {
      if(nrow(x)>1) {
        return(x)
      } else {
        return(NULL)
      }
    })
    ### remove the NULL values 
    same.end = same.end[!sapply(same.end, is.null)]
    same.end <- lapply(same.end,function(x){ x$share = 'same_end'; return(x)})
    return(c(same.start, same.end))
  }, mc.cores = NT))
  ## junction just has two alternative starts/ends:
  junc.as <- junc.as[sapply(junc.as, nrow) > 1]
  names(junc.as) <- sapply(junc.as, function(x) paste(x$names, collapse="@"))
  return(junc.as)
}
##### create indexed gtf (extract the first exon/intron) ------

##GTF index outputs with intron exon 5utr promoter regions
#' @export
GTFIndex <- function(gtfFile, NT = 1){
  if(!file.exists(gtfFile)) {
    stop("The GTF file does not exist! And you should provide the GTF file you used in STAR alignment.")
  }
  
  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")
  useBiocPackage('GenomicFeatures')
  #### contrust the database from GTF -----------
  txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtfFile, format = "gtf")
  #### the intron regions to map gene id -------------
  intron <- intronicParts(txdb, linked.to.single.gene.only = TRUE)
  intron_tab <- as.data.table(intron)
  intron_tab[, names:=as.character(intron)]
  # AF/LE -----------------------------------
  txintron <- intronsByTranscript(txdb, use.names=TRUE)
  txintron <- txintron[S4Vectors::elementNROWS(txintron) > 0]
  s_s <- as.character(unlist(unique(strand(txintron))))
  s_l <- S4Vectors::elementNROWS(txintron)
  txintron_gr <- unlist(txintron)
  my.rank <- function(s, l) ifelse(s == "+", return(1:l), return(-(l:1)))
  GenomicRanges::mcols(txintron_gr)$intron_rank <- do.call(c, lapply(seq_along(s_s), function(i) my.rank(s_s[i], s_l[i])))
  GenomicRanges::mcols(txintron_gr)$tx_id <- names(txintron_gr)
  txintron_tab_rank <- as.data.table(txintron_gr)
  
  first_intron_tab <- txintron_tab_rank[abs(intron_rank) == 1, ]
  first_intron_tab[, junc:=paste0(seqnames, ":", start, "-", end, ":", strand)]
  setkey(first_intron_tab, tx_id)
  
  last_intron_tab <- txintron_tab_rank[ , .SD[abs(intron_rank) == max(abs(intron_rank)), ], by = "tx_id"]
  last_intron_tab[, junc:=paste0(seqnames, ":", start, "-", end, ":", strand)]
  setkey(last_intron_tab, tx_id)
  
  # ASS -------------------------------------
  exons <- as.data.table(unique(exons(txdb, columns = c("gene_id", "EXONNAME"))))
  ### find first exon
  exon_tab <- as.data.table(unique(GenomicFeatures::exonicParts(txdb)))
  exon_tab <- exon_tab[seqnames %in% c(1:22,'X','Y')][,exon_rank:=as.numeric(as.character(exon_rank))]
  first_exon_tab <- exon_tab[abs(exon_rank) == 1, ][,gene_id:=as.character(gene_id)][order(gene_id)]
  first_exon_tab = transform(first_exon_tab, seqnames = as.character(seqnames), start=as.numeric(start), end = as.numeric(end))
  ### obtain the first exon boundary site
  ### end+1 for + strand
  first_exon_tab[,exon_boundary:=paste0(seqnames,':',end+1,':+')]
  ### start-1 for - strand
  first_exon_tab[strand=='-'] = transform(first_exon_tab[strand=='-'], exon_boundary=paste0(seqnames,':',start-1,':-'))
  ### we retain the gene with multiple transcripts
  retained_gene = first_exon_tab[duplicated(gene_id)]$gene_id
  first_exon_tab = first_exon_tab[gene_id%in%retained_gene]
  setkey(first_exon_tab, gene_id)
  
  #### contrust the database from GTF -----------
  txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtfFile, format = "gtf")
  #### obtain the 5 UTR information
  fiveUTR=fiveUTRsByTranscript(txdb)
  ###Creat the RLe of 5' region
  fiveUTRs_grl = GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  # Create Rle of the tx_names
  fiveUTRs_txname_rle = S4Vectors::Rle(names(fiveUTRs_grl), S4Vectors::elementNROWS(fiveUTRs_grl))
  fiveUTRs_txname_vec = as.character(fiveUTRs_txname_rle)
  # Unlist and add the tx_names
  fiveUTRs_gr = unlist(fiveUTRs_grl, use.names = FALSE)
  GenomicRanges::mcols(fiveUTRs_gr)$tx_name = fiveUTRs_txname_vec
  # Add Entrez ID, symbol, and type
  # NOTE: here we match on the tx_name because the tx_id is not given
  #GenomicRanges::mcols(fiveUTRs_gr)$gene_id = id_maps[match(GenomicRanges::mcols(fiveUTRs_gr)$tx_name, id_maps$TXNAME), 'GENEID']
  #GenomicRanges::mcols(fiveUTRs_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(fiveUTRs_gr)$gene_id, eg2symbol$symbol), 'symbol']
  #GenomicRanges::mcols(fiveUTRs_gr)$entrez = eg2symbol[match(GenomicRanges::mcols(fiveUTRs_gr)$gene_id, eg2symbol$symbol),'gene_id']
  GenomicRanges::mcols(fiveUTRs_gr)$type = sprintf('5UTRs')
  GenomicRanges::mcols(fiveUTRs_gr)$id = paste0('5UTR:ExonRank', fiveUTRs_gr$exon_rank)
  pos=fiveUTRs_gr[strand(fiveUTRs_gr)=='+',]
  neg=fiveUTRs_gr[strand(fiveUTRs_gr)=='-',]
  ###Create Rle of the promoter region
  promoters_gr = GenomicFeatures::promoters(txdb, upstream = 1000, downstream = 200)
  # Add Entrez ID, symbol, and type
  #GenomicRanges::mcols(promoters_gr)$gene_id = id_maps[match(GenomicRanges::mcols(promoters_gr)$tx_id, id_maps$TXID), 'GENEID']
  #GenomicRanges::mcols(promoters_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(promoters_gr)$gene_id, eg2symbol$symbol), 'symbol']
  #GenomicRanges::mcols(promoters_gr)$type = sprintf('%s_genes_promoters', genome)
  GenomicRanges::mcols(promoters_gr)$id = paste0('promoter:', seq_along(promoters_gr))
  #GenomicRanges::mcols(promoters_gr) = GenomicRanges::mcols(promoters_gr)[, c('id','tx_name','gene_id','symbol','type')]
  #colnames(GenomicRanges::mcols(promoters_gr)) = c('id','tx_id','gene_id','symbol','type')
  
  
  
  #### the intron regions to map gene id -------------
  intron <- intronicParts(txdb, linked.to.single.gene.only = TRUE)
  intron_tab <- as.data.table(intron)
  intron_tab[, names:=as.character(intron)]
  
  ### classify the distinct exons and
  ### tandem exons in + strand (shared same start, with different ends)
  #first_exon_tab_pos = first_exon_tab[strand=='+',.(startl = length(start), startu = length(unique(start)), endl=length(end), endu = length(unique(end))),by=gene_id]
  #first_exon_tab_pos[,exon_type:=ifelse(startl == startu & endu < endl, 'Tandem_exon', 'Distinct_exon')] 
  ### tandem exons in - strand (shared same end, with different starts)
  #first_exon_tab_neg = first_exon_tab[strand=='-',.(startl = length(start), startu = length(unique(start)), endl=length(end), endu = length(unique(end))),by=gene_id]
  # first_exon_tab_pos[,exon_type:=ifelse(startu < startl & endu == endl, 'Tandem_exon', 'Distinct_exon')] 
  
  return(list(intron_tab = intron_tab, first_intron_tab = first_intron_tab, last_intron_tab = last_intron_tab,
              fiveUTRs_gr = fiveUTRs_gr, exons = exons, first_exon_tab = first_exon_tab,promoters_gr=promoters_gr))
}
####




# normally sj = unique(rownames())
#### find the host gene for AS events using ID------
#' @export
map_HostGene_for_ID <- function(sj, gtfFile, NT = 1){
  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")
  useBiocPackage('GenomicFeatures')
  useBiocPackage('rtracklayer')
  
  print(paste("Starting map host genes at ", Sys.time(), sep = ""))
  if(!file.exists(gtfFile)) {
    stop("The GTF file does not exist! And you should provide the GTF file you used in STAR alignment.")
  }
  
  
  # gene id ---------------------------------
  sj_tu <- unique(do.call(c, lapply(sj, function(x) unlist(strsplit(x, "@")))))
  junc_tu <- as.data.table(data.frame(names = sj_tu,
                                      chr = mapply(function(x) as.character(x[1]), strsplit(sj_tu, ":")),
                                      start = mapply(function(x) as.integer(x[2]), strsplit(sj_tu, "[:-]")),
                                      end = mapply(function(x) as.integer(x[3]), strsplit(sj_tu, "[:-]")),
                                      strand = mapply(function(x) as.character(x[3]), strsplit(sj_tu, "[:]"))))
  junc_tu <- merge(junc_tu, intron_tab[, .(names, gene_id)], by = "names", all.x = TRUE)
  
  ### treat the junctions not mapped to intergenic (find the ranking genes) --------
  tmp <- junc_tu[is.na(gene_id), ]
  mclapply(seq_len(nrow(tmp)), function(i){
    return(GRanges(seqnames = tmp[i, chr], ranges = IRanges(start = tmp[i, start], end = tmp[i, end]), strand = tmp[i, strand]))
  }, mc.cores = NT) -> res_s
  grl_s <- GRangesList(res_s); names(grl_s) <- tmp[, names]
  res_s_gene <- mapRangesToIds(txdb, grl_s, "gene")
  range_gene <- as.data.table(data.frame(gene_id = mcmapply(function(x) paste(x$gene_id, collapse = "|"), res_s_gene, mc.cores = NT)), keep.rownames = "junction")
  tmp[range_gene$junction, ]$gene_id <- range_gene$gene_id
  
  junc_tu_gene <- rbind(na.omit(junc_tu), tmp)
  setkey(junc_tu_gene, names)
  
  # gene_name -------------------------------
  tu_gene_id <- as.list(junc_tu_gene$gene_id); names(tu_gene_id) <- junc_tu[, names]
  gtf <- rtracklayer::readGFF(gtfFile)
  gtf_infor <- as.data.table(gtf[gtf$type == "gene", c("gene_id", "gene_name")])
  setkey(gtf_infor, gene_id)
  junc_tu_gene$gene_name <- do.call(c, lapply(tu_gene_id, function(x) paste(gtf_infor[unlist(strsplit(x, "\\|")), gene_name], collapse = "|") ))
  
  # add host gene ---------------------------
  do.call(rbind, mclapply(sj, function(x) {
    as.data.frame(t(apply(junc_tu_gene[unlist(strsplit(x, "@")), .(names, gene_id, gene_name)], 2, function(y) paste0(y, collapse = "@"))))
  }, mc.cores = NT)) -> psi_id_gene
  setDT(psi_id_gene, key = "names")
  psi_id_gene[, names:=as.character(names)]
  colnames(psi_id_gene)[1] = 'id'
  print(paste("Finished at ", Sys.time(), sep = ""))
  return(psi_id_gene)
}
#### optimization of identifying AFE ---------
#tandem first exon (with different length)
#A5SS & AFE in + strand
#A3SS & AFE in - strand 
####The canonical AFE:
#### the alternative exon are located in different transcripts from the sam gene 
#### map the first exon's end wiith junction (+);  the first exon's start wiith junction (-)
#' @export
identify_AStype <- function(cell.psi, sj_list, indexedGTF, NT = 1){
  ### cell_psi: the result table from betabinomial test
  ### sj_list: the list of AS junctions
  ### indexedGTF: index files' list from GTF
  # required package ------------------------
  usePackage("data.table")
  usePackage("parallel")
  
  junc_as2 <- as.list(cell.psi[cell.psi$event == 2, 'id'])
  junc_as <- do.call(rbind, lapply(junc_as2, function(x) {
    unlist(strsplit(stringr::str_replace_all(x, "-(?=[[:digit:]])", "@"), "[:@]"))
  }))
  junc_as <- as.data.table(junc_as)
  colnames(junc_as) <- paste0(c("chr", "start", "end", "strand"), rep(1:2, c(4,4)))
  junc_as[, names:=cell.psi[cell.psi$event == 2, 'id']]
  
  same_start_junc2 <- junc_as[start1 == start2, ]
  same_end_junc2 <- junc_as[end1 == end2, ]
  same_start_junc2 <- same_start_junc2[, .(names, chr1, start1, end1, end2, strand1)]
  colnames(same_start_junc2) <- c("names", "chr", "start", "end1", "end2", "strand")
  same_end_junc2 <- same_end_junc2[, .(names, chr1, start1, start2, end1, strand1)]
  colnames(same_end_junc2) <- c("names", "chr", "start1", "start2", "end", "strand")
  
  
  
  match_nega <- same_start_junc2[strand=="-", same_end_junc2[strand=="-",], by = c("names", "chr", "start", "end1", "end2", "strand")]
  match_posi <- same_start_junc2[strand=="+", same_end_junc2[strand=="+",], by = c("names", "chr", "start", "end1", "end2", "strand")]
  start_end_match <- rbind(match_posi, match_nega)
  colnames(start_end_match)[colnames(start_end_match) == "seqnames"] <- c("seqnames", "seqnames1")
  start_end_match <- start_end_match[as.logical(start_end_match[, seqnames]==start_end_match[, seqnames1]), ]
  
  # SE --------------------------------------
  SE <- start_end_match[start==start2 & end==end2 & end1 < start1,]
  SE_junc <- unique(c(SE[[1]], SE[[7]]))
  
  # MXE -------------------------------------
  MXE <- start_end_match[start < end1 & end1 < start2 & start2 < end2 & end2 < start1 & start1 < end, ]
  MXE_junc <- unique(c(MXE[[1]], MXE[[7]]))
  
  # AF/LE (special we need use sj_list)-----------------------------------
  first_intron_tab <- indexedGTF$first_intron_tab
  last_intron_tab <- indexedGTF$last_intron_tab
  first_exon_tab <- indexedGTF$first_exon_tab
  
  afle_index <- mcmapply(sj_list, FUN = function(x){
    if(sum(x$names %in% first_intron_tab[, junc]) > 0 ){
      res = "Pseudo_AFE"
      if(sum(x$first_boundary %in% first_exon_tab$exon_boundary)>1){
        if(x$strand == '+' & x$share == 'same_end'){
          res = 'Real_AFE' 
        } else if(x$strand == '-' & x$share == 'same_start') {
          res = 'Real_AFE' 
        }
      }
    } else if (sum(x$names %in% last_intron_tab[, junc]) > 0 ) {
      res = "ALE"
    } else {
      res = "Non"
    }
    return(res)
  }, mc.cores = NT)
  #AFE_junc <- cell.psi[, as.character(id)][grep('AFE',afle_index)]
  R_AFE_junc <- cell.psi[, 'id'][which(afle_index == "Real_AFE")]
  P_AFE_junc <- cell.psi[, 'id'][which(afle_index == "Pesudo_AFE")]
  ALE_junc <- cell.psi[, 'id'][which(afle_index == "ALE")]
  
  
  # A3/5SS ----------------------------------
  same_start_junc2_gr <- GRanges(seqnames = same_start_junc2[, chr],
                                 ranges = IRanges(start = as.integer(same_start_junc2[, end1]) + 1, end = as.integer(same_start_junc2[, end2])),
                                 strand = same_start_junc2[, strand],
                                 names = same_start_junc2[, names])
  same_end_junc2_gr <- GRanges(seqnames = same_end_junc2[, chr],
                               ranges = IRanges(start = as.integer(same_end_junc2[, start2]), end = as.integer(same_end_junc2[, start1]) - 1),
                               strand = same_end_junc2[, strand],
                               names = same_end_junc2[, names])
  tmp1 <- as.data.table(same_start_junc2_gr)
  tmp2 <- indexedGTF$exons
  tmp3 <- as.data.table(same_end_junc2_gr)
  
  same_start_ass <- merge(tmp1, tmp2, by = c("seqnames", "start", "strand"))[end.x < end.y, ]
  same_end_ass <- merge(tmp3, tmp2, by = c("seqnames", "end", "strand"))[start.x > start.y, ]
  A3SS_junc <- unique(c(same_start_ass[strand == "+", names], same_end_ass[strand == "-", names]))
  A5SS_junc <- unique(c(same_start_ass[strand == "-", names], same_end_ass[strand == "+", names]))
  
  AS_list <- list(SE = SE_junc, A3SS = A3SS_junc, A5SS = A5SS_junc,
                  Real_AFE = R_AFE_junc, Pseudo_AFE = P_AFE_junc, ALE = ALE_junc, MXE = MXE_junc)
  
  # Add AS type ------------------------
  as.type <-  function(x) {
    if(sum(mapply(function(y) is.element(x, y), AS_list)) > 0){
      AS = paste0(names(AS_list)[mapply(function(y) is.element(x, y), AS_list)], collapse = "|")
    } else {
      AS = "Other"
    }
    return(AS)
  }
  
  cell.psi[!is.na(event), AS:=mcmapply(as.type, cell.psi[!is.na(event), id], mc.cores = NT)]
  ### assign #A5SS & AFE in + and A3SS & AFE in - to tandem AFE
  
  cell.psi[grep('+',cell.psi$id) & cell.psi$AS=='A5SS|Real_AFE','AS'] = 'Tandem_AFE'
  cell.psi[grep('+',cell.psi$id,invert = T) & cell.psi$AS=='A5SS|Real_AFE','AS'] = 'Tandem_AFE'
  return(cell.psi)
}


















