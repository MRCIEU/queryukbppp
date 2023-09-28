.onAttach <- function(libname, pkgname) {
    bc4dir <- "/mnt/storage/private/mrcieu/research/scratch/IGD/data/dev/UKB-PPP"
    if(file.exists(bc4dir)) {
        packageStartupMessage("Found bc4 directory: ", bc4)
        options(ukbpppdir=bc4dir)
    }
}

