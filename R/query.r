#' Query UKB-PPP downloaded files
#' 
#' The files are arranged according to https://www.synapse.org/#!Synapse:syn51364943/files/.
#' The \code{metadata/} directory has a mapping of rsid and chr/pos in build hg19 and hg38.
#' The \code{UKB-PPP pGWAS summary statistics} directory has a directory for each ancestry
#' with complete summary statistics for ~3k proteins in each directory. Each GWAS result is 
#' split into separate chromosome files and then tarballed.
#' 
#' This function takes a data frame of variants and extracts them from a specified tarball file
#' 
#' @param tarfile Path to a tarball e.g. \code{/path/to/African/A1BG_P04217_OID30771_v1_Inflammation_II.tar}
#' @param variants Output from \code{get_snpid_from_chrpos} or \code{get_snpid_from_rsid}
#' 
#' @return data frame
#' @export
query_pgwas <- function(tarfile, variants) {
    td <- tempdir()
    cmd <- glue::glue("tar xvf {tarfile} -C {td}")
    system(cmd)
    fd <- file.path(td, gsub(".tar", "", basename(tarfile)))
    fnut <- list.files(fd, full.names=TRUE)
    l <- list()
    for(ch in unique(variants$chr))
    {
        message(ch)
        fn <- grep(paste0("chr", ch, "_"), fnut, value=T)
        if(!file.exists(fn))
        {
            message("missing")
            next
        }
        d <- lookup_txt(fn, subset(variants, chr == ch)$ID)
        l[[ch]] <- d
    }
    l <- bind_rows(l)
    con <- gzfile(fn)
    names(l) <- scan(con, nlines=1, what=character())
    close(con)
    l$tarfile <- tarfile
    system(glue::glue("rm -r {td}"))
    return(l)
}

#' Get variant information based on chr/pos lookup
#' 
#' Provide a path to the \code{Metadata/SNP RSID maps} directory and lookup the snp IDs by chromosome position
#' 
#' @param chrpos vector of variants by chr/pos e.g. \code{c("10:60494", "10:60515", "10:60523")}
#' @param build build of chrpos being queries - 19 or 38
#' @param map_files path to the \code{Metadata/SNP RSID maps} directory
#' 
#' @return data frame
#' @export
get_snpid_from_chrpos <- function(chrpos, build=c("19", "38")[1], map_files=get_mapfiles()) {
    td <- tempdir(check=TRUE)
    chrpos <- strsplit(chrpos, ":")
    chrpos <- dplyr::tibble(chr=sapply(chrpos, \(x) x[1]), pos=sapply(chrpos, \(x) x[2]))
    l <- list()
    for(ch in unique(chrpos$chr)) {
        message(ch)
        data.table::fwrite(subset(chrpos, chr == ch)$pos, file=file.path(td, "pos"), row.names=FALSE, col.names=FALSE, quote=FALSE)
        cmd <- glue::glue("zgrep -wf {file.path(td, 'pos')} {grep(paste0('chr', ch), map_files, value=TRUE)} > {file.path(td, 'out')}")
        system(cmd)
        l[[ch]] <- data.table::fread(file.path(td, 'out'))
    }
    a <- bind_rows(l)
    names(a) <- c("ID", "REF", "ALT", "rsid", "POS19", "POS38")
    if(nrow(a) > 0) {
        a$chr <- strsplit(a$ID, ":") %>% sapply(., \(x) x[1])
        if(build == "19") {
            a <- subset(a, paste(chr, POS19) %in% paste(chrpos$chr, chrpos$pos))
        } else {
            a <- subset(a, paste(chr, POS38) %in% paste(chrpos$chr, chrpos$pos))
        }
    }
    return(a)    
}

#' Get variant information based on rsid lookup
#' 
#' Provide a path to the \code{Metadata/SNP RSID maps} directory and lookup the snp IDs by rsid
#' 
#' @param rsid vector of variants by rsid
#' @param map_files path to the \code{Metadata/SNP RSID maps} directory
#' 
#' @return data frame
#' @export
get_snpid_from_rsid <- function(rsid, map_files=get_mapfiles()) {
    td <- tempdir(check=TRUE)
    data.table::fwrite(rsid, file=file.path(td, "rsid"), row.names=FALSE, col.names=FALSE, quote=FALSE)
    map_files <- paste0()
    cmd <- glue::glue("zgrep -wf {file.path(td, 'rsid')} {paste(map_files, collapse=' ')} > {file.path(td, 'out')}")
    system(cmd)
    a <- data.table::fread(file.path(td, 'out'))
    names(a) <- c("ID", "REF", "ALT", "rsid", "POS19", "POS38")
    if(nrow(a) > 0) {
        a$ID <- strsplit(a$ID, ".gz:") %>% sapply(., \(x) x[2])
        a$chr <- strsplit(a$ID, ":") %>% sapply(., \(x) x[1])
    }
    return(a)
}


lookup_txt <- function(fn, pos) {
    tf <- tempfile()
    tf2 <- tempfile()
    data.table::fwrite(unique(pos), file=tf, row.names=FALSE, col.names=FALSE, quote=FALSE)
    cmd <- glue::glue("zgrep -wf {tf} {paste(fn, collapse=' ')} > {tf2}")
    system(cmd)
    data.table::fread(tf2)
}

#' Get map file paths
#' 
#' @param d Directory (defaults to bc4 location)
#' 
#' @return array of paths for each map file chromosome
#' @export
get_mapfiles <- function(d=options()$ukbpppdir) {
    file.path(d, "Metadata", "SNP RSID maps") %>% list.files(full.names=TRUE)
}

#' Get protein information
#' 
#' Adds the tarball file name column too
#' 
#' @param d base directory (defaults to bc4 location)
#' 
#' @return data frame
#' @export
protein_info <- function(d=options()$ukbpppdir) {
    file.path(d, "Metadata", "Protein annotation", "olink_protein_map_3k_v1.tsv") %>%
        data.table::fread(.) %>%
        mutate(tarfile = paste0(gsub(":", "_", UKBPPP_ProteinID), "_", Panel, ".tar"))
}
