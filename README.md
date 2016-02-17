# pqsfinder package

This package provides functions for identification of potential intramolecular quadruplex patterns in DNA sequence. The main functionality is to detect the positions of subsequences capable of folding into an intramolecular G4-quadruplex in a much larger sequence.

## How to Install

Download the latest develompment version of R:
```bash
wget https://stat.ethz.ch/R/daily/R-devel.tar.gz
```

Unpack:
```bash
tar -xzf R-devel.tar.gz
```

Configure with `--enable-R-shlib` option in case you want to run the newest R inside RStudio. It is possible, you will be required to install additional libraries. The corresponding package names are usually in the form `libXXX-dev`, where `XXX` is the name of the library.
```bash
cd R-devel
./configure --enable-R-shlib=yes
```

**Only** compile. There is no reason to run `make install` as the resulting binary can be run directly from the compilation folder without messing up your system.
```bash
make
```

Run:
```bash
./bin/R
```

Install Bioconductor and `pqsfinder` package dependencies:
```R
source("https://bioconductor.org/biocLite.R")
biocLite()

biocLite(c("Rcpp", "GenomicRanges", "Biostrings", "flowCore", "devtools"))
```

Install `pqsfinder` package directly from GitHub:
```R
library("devtools")
install_github("jirihon/pqsfinder")
```
