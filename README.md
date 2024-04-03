# MSstatsBig

MSstatsBig provides support for very large data files which either do not fit 
into a standard computer's memory, or take a very long time to process. These 
data files are generally the output of spectral processing tools, such as 
Spectronaut or FragPipe. This workflow was designed with DIA experiments in mind.

The package includes converters which replace the standard MSstats converters 
(such as SpectronauttoMSstatsFormat). After the converters are run, the 
standard MSstats workflow can be followed.

## Installation 

This package is currently available on Bioconductor and Github.

It can be installed from Bioconductor as follows:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsBig")
```

It can be installed from Github as follows:

```
devtools::install_github("Vitek-Lab/MSstatsBig", build_vignettes = TRUE)
```

Important functions:

```
BigFragPipetoMSstatsFormat()
BigSpectronauttoMSstatsFormat()
MSstatsPreprocessBig()
```

## License

[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0)
