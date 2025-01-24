# TSO500R(eader)

[![R-CMD-check](https://github.com/Boehringer-Ingelheim/TSO500R/actions/workflows/test.yml/badge.svg)](https://github.com/Boehringer-Ingelheim/TSO500R/actions/workflows/test.yml)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/Boehringer-Ingelheim/TSO500R/blob/master/LICENSE)

#### Handle TSO500 data using R.

---

*TSO500R(eader)* is an R package developed for handling Illumina [TruSight Oncology 500](https://emea.illumina.com/products/by-type/clinical-research-products/trusight-oncology-500.html) data. It can be used for importing and processing of files produced by the Illumina [TSO500 DRAGEN analysis pipeline](https://support-docs.illumina.com/SW/DRAGEN_TSO500_v2.1/Content/SW/FrontPages/DRAGENTSO500_v2.1.htm), [TSO500 DRAGEN ctDNA analysis pipeline](https://support-docs.illumina.com/SW/DRAGEN_TSO500_ctDNA_v2.1/Content/SW/FrontPages/DRAGENTSO500_ctDNA_v2.1.htm) and the [LocalApp](https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/trusight/trusight-oncology-500/trusight-oncology-500-local-app-v2.2-user-guide-1000000137777-01.pdf).

The package provides different functions for parsing the various output files produced by the Illumina pipelines. This includes quality control files as well as the analysis outputs. Besides, it offers functionality to integrate the different result types, e.g. small variants and amplifications.

Other features include functions for basic plotting and export functionality for writing `RData` objects or to generate DRAGEN analysis pipeline samplesheets.

## Installation

Clone the repository:
```bash
git clone https://github.com/Boehringer-Ingelheim/TSO500R.git
```

Install from path using R:
```r
# Install TS500R from local clone
install.packages("/path/to/TSO500R", repos = NULL)
```

## Usage

TSO500R provides an API for interaction with most of the output files produced by the TSO500 analysis pipeline. In addition the data in each section of the quality control file (`MetricsOutput`) and the main output file of the TSO500 analysis, the `CombinedVariantOutput.tsv` file, can be read in separately.

Once installed, you can load the TSO500R package. 

```r
# load the tso500R(eader) package
library(tso500R)
```

Afterwards, you can, e.g., load all files of type `MetricsOutput` in a given directory.

```r
# load the MetricsOutput.tsv file in a specific folder
tso500_data_dir <- "/path/to/your/TSO500/analysis/output/data/"
qmo_data <- read_qmo_data(tso500_data_dir)
```

For detailed instructions and examples, please refer to the [documentation](https://boehringer-ingelheim.github.io/tso500R/).

## Contributions & Support 

Suggestions for improvements and new features, as well as contributions are welcome. Please use the [issue tracker](https://github.com/Boehringer-Ingelheim/TSO500R/issues) to request features and for bug reports.

Also feel free to get in touch directly by contacting [@christopher-mohr](https://github.com/christopher-mohr).

## Credits

TODO