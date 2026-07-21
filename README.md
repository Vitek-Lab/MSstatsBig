# MSstatsBig

<!-- badges: start -->
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/MSstatsBig.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MSstatsBig/)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstatsBig/branch/devel/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstatsBig/branch/devel)
[![Bioconductor Downloads Rank](https://bioconductor.org/shields/downloads/release/MSstatsBig.svg)](https://bioconductor.org/packages/stats/bioc/MSstatsBig/)
[![Years in Bioconductor](https://bioconductor.org/shields/years-in-bioc/MSstatsBig.svg)](https://bioconductor.org/packages/release/bioc/html/MSstatsBig.html#since)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->

MSstatsBig is an R/Bioconductor package for preprocessing mass spectrometry-based
proteomics experiments that are larger than a standard computer's memory. Modern
acquisition protocols — most often data-independent acquisition (DIA) with a
large number of MS runs — can produce quantitative reports that the standard
MSstats converters cannot load into RAM. MSstatsBig provides drop-in replacement
converters that filter each tool's PSM report down to the data required for
differential analysis using out-of-memory backends (Arrow / Spark), so the
result fits in memory. Once a converter has run, the standard MSstats workflow
(`dataProcess()`, `groupComparison()`) can be followed as usual.

MSstatsBig is part of the [MSstats](https://github.com/Vitek-Lab/MSstats)
family of packages, developed and maintained by the
[Vitek Lab](https://olga-vitek-lab.khoury.northeastern.edu/) at Northeastern
University. The package and its documentation are also available at
[msstats.org](http://msstats.org).

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsBig")
```

The development version can be installed directly from this repository:

```r
BiocManager::install("Vitek-Lab/MSstatsBig", ref = "devel")
```

## Quick Start

```r
library(MSstatsBig)

# Example FragPipe msstats.csv output bundled with the package
input <- system.file("extdata", "fgexample.csv", package = "MSstatsBig")

# Convert out-of-memory using the Arrow backend; writes to output_file.csv
converted <- bigFragPipetoMSstatsFormat(input,
                                        "output_file.csv",
                                        backend = "arrow",
                                        max_feature_count = 20)

# Collect the Arrow object into an in-memory data.frame
converted <- as.data.frame(dplyr::collect(converted))

# Continue with the standard MSstats workflow
library(MSstats)
summarized <- dataProcess(converted, use_log_file = FALSE)
```

Users can also feed data into the underlying `MSstatsPreprocessBig()` function
directly — either from a tool's native export format or by converting raw data
chunk by chunk — to support other tools such as DIA-NN. See the workflow
vignette for details.

## Supported Converters

MSstatsBig provides larger-than-memory replacements for the standard MSstats
converters. After a converter runs, its output is used with
[MSstats](https://github.com/Vitek-Lab/MSstats):

| Search tool / format | Converter function |
| --- | --- |
| FragPipe | `bigFragPipetoMSstatsFormat()` |
| Spectronaut | `bigSpectronauttoMSstatsFormat()` |
| DIA-NN | `bigDIANNtoMSstatsFormat()` |
| Generic MSstats-format data | `MSstatsPreprocessBig()` |

The `MSstatsAddAnnotationBig()` helper attaches experimental-design annotation
to converted data. See the [MSstatsBig workflow vignette](vignettes/MSstatsBig_Workflow.Rmd)
for the required input files and options for each converter.

## Documentation

- [MSstatsBig workflow](vignettes/MSstatsBig_Workflow.Rmd) — full worked example from raw data to results
- [Official website: msstats.org](http://msstats.org)
- [Bioconductor package page and reference manual](https://bioconductor.org/packages/MSstatsBig)

## Getting Help / Reporting Bugs

- **Questions about usage, statistical methods, or troubleshooting:** please
  post to the [MSstats Google Group](https://groups.google.com/forum/#!forum/msstats).
  This is monitored by the development team and searchable, so it's the fastest
  way to get help and to see if your question has already been answered.
- **Bug reports and feature requests for this repository:** please open a
  [GitHub issue](https://github.com/Vitek-Lab/MSstatsBig/issues).

## References

MSstatsBig prepares data for the standard MSstats workflow. If you use it, please
cite the MSstats publication(s):

1. Kohler D, Staniak M, Tsai TH, Huang T, Shulman N, Bernhardt OM, MacLean BX,
   Nesvizhskii AI, Reiter L, Sabido E, Choi M, Vitek O. **MSstats Version 4.0:
   Statistical Analyses of Quantitative Mass Spectrometry-Based Proteomic
   Experiments with Chromatography-Based Quantification at Scale.**
   *J Proteome Res*. 2023;22(5):1466-1482.
   [DOI: 10.1021/acs.jproteome.2c00834](https://doi.org/10.1021/acs.jproteome.2c00834)
2. Kohler D, Vitek O, et al. **An MSstats workflow for detecting differentially
   abundant proteins in large-scale data-independent acquisition mass
   spectrometry experiments with FragPipe processing.** *Nat Protoc*.
   2024;19:2915-2938.
   [DOI: 10.1038/s41596-024-01000-3](https://doi.org/10.1038/s41596-024-01000-3)

## Funding

MSstats development has been supported by the Chan Zuckerberg Initiative's
[Essential Open Source Software for Science](https://chanzuckerberg.com/eoss/proposals/).

## License

MSstatsBig is released under the [Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
