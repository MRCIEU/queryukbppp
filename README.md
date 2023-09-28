# Convenience functions for extracting UKB-PPP data

Assumed that the data is synced from here https://www.synapse.org/#!Synapse:syn51364943/files/

e.g. on bc4 the path is

```
/mnt/storage/private/mrcieu/research/scratch/IGD/data/dev/UKB-PPP
```

## Installation

```r
install.packages("remotes")
remotes::install_github("mrcieu/queryukbppp")
```

## Example

To extract a set of variants from e.g. all African summary stat files on bc4

### 1. Get variant IDs

```r
library(queryukbppp)
rsid <- c("rs568182971", "rs566854205", "rs112920234", "rs552582111", "rs569167217")
v <- get_snpid_from_rsid(rsid)
```

### 2. Extract from one protein

```r
prot <- protein_info()
tarfile <- file.path(options()$ukbpppdir, "UKB-PPP pGWAS summary statistics", "African", prot$tarfile[1])
out <- query_pgwas(tarfile, v)
```

### 3. Extract variants from all proteins

```r
out <- lapply(
    file.path(options()$ukbpppdir, "UKB-PPP pGWAS summary statistics", "African", prot$tarfile), 
    \(x) {
        message(x)
        query_pgwas(x, v)
    }    
)
```

