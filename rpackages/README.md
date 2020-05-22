# R packages

Required [R](https://www.r-project.org/) packages `adegenet`, `adegraphics`, `StAMPP` and `vcfR` and their dependencies for R version used **must** be installed here.

Installation of source packages in Linux requires basic compilation tools and developmental files for GDAL, GEOS and UDUNITS-2.

Within `R` command line started in `radseq` directory use e.g. command

```R
install.packages(pkgs=c("adegenet", "adegraphics", "memuse", "BH", "pegas",
"permute", "plogr", "sf", "spData", "StAMPP", "testthat", "vcfR"),
lib="rpackages", repos="https://mirrors.nic.cz/R/", dependencies="Imports")
```

to install needed packages.

