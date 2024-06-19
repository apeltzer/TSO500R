# TSO500R(eader)

## Installation

Clone the repository and run:

```r
# Install TS500R from local clone
install.packages("/path/to/TSO500R", repos = NULL)
```

## Overview

*TSO500R(eader)* is an R package developed for Illumina TSO500 data. It can be used for importing and processing of files produced by the Illumina [TSO500 DRAGEN analysis pipeline](https://support-docs.illumina.com/SW/DRAGEN_TSO500_v2.1/Content/SW/FrontPages/DRAGENTSO500_v2.1.htm), [TSO500 DRAGEN ctDNA analysis pipeline](https://support-docs.illumina.com/SW/DRAGEN_TSO500_ctDNA_v2.1/Content/SW/FrontPages/DRAGENTSO500_ctDNA_v2.1.htm) and the [LocalApp](https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/trusight/trusight-oncology-500/trusight-oncology-500-local-app-v2.2-user-guide-1000000137777-01.pdf).

The package provides different functions for parsing the various output files produced by the Illumina pipelines. This includes the quality control files as well as the analysis outputs. Besides, it offers functionality to integrate the different result types, e.g. small variants and amplifications.

Other features include functions for basic plotting and export functionality for writing `RData` objects or to generate DRAGEN analysis pipeline samplesheets.

## Usage

TSO500R provides an API for interaction with most of the output files produced by the TSO500 analysis pipeline. In addition the data in each section of the quality control file (`MetricsOutput`) and the main output file of the TSO500 analysis, the `CombinedVariantOutput.tsv` file, can be read in separately.

```r
# load the tso500R(eader) package
library(tso500R)
```

- [Loading MetricsOutput Data](#loading-metricsoutput-data)
- [Loading CombinedVariantOutput Datas](#loading-combinedvariantoutput-data)
  - [Accessing Analysis Details](#accessing-analysis-details)
  - [Accessing TMB and MSI Metrics](#accessing-tmb-and-msi-metrics)
  - [Accessing Small Variants](#accessing-small-variants)
- [Integration of Data](#integration-of-data)
- [Plotting](#plotting)
- [Other Convenience Functions](#other-convenience-functions)


### Loading MetricsOutput Data

In order to load the `MetricsOutput.tsv` file or multiple such output files of a DRAGEN analysis run, that contains all quality control metrics of a run and all included samples, use the following code:

```r
# load the MetricsOutput.tsv file in a specific folder
tso500_data_dir <- "/path/to/your/TSO500/analysis/output/data/"
qmo_data <- read_qmo_data(test_path)
```

In case you want to read in data that resulted from a `LocalApp` run, you need to set `local_app=TRUE` due to the (slight) differences between the file formats:

```r
# load the MetricsOutput.tsv file of a LocalApp run in a specific folder
tso500_data_dir <- "/path/to/your/TSO500/analysis/output/data/"
qmo_data <- read_qmo_data(test_path, local_app=TRUE)
```

In case you want to read in data that resulted from a `ctDNA analysis` run, you need to set `ctdna=TRUE` due to the (slight) differences between the file formats:

```r
# load the MetricsOutput.tsv file of a LocalApp run in a specific folder
tso500_data_dir <- "/path/to/your/TSO500/analysis/output/data/"
qmo_data <- read_qmo_data(test_path, ctdna=TRUE)
```

This will produce a list of `combined.quality.metrics.output` objects, which can be further processed and contain all the information as given in the parsed files. Each section in the `CombinedVariantOutput` file can be accessed with standard list indexing. The `combined.quality.metrics.output` object in the last can be accessed via the corresponding file name (without extension). This will therfore probably be `MetricsOutput` in case you read in a metrics file of all samples of a run and will include sample identifiers (`sample1_MetricsOutput`) if you want to read in individual metrics file.

```r
# access the header information of a specific MetricsOutput file
qmo_data$filename$header

# access the qc metrics on DNA of a specific MetricsOutput file
qmo_data$filename$dna_qc_metrics

# access the qc metrics on DNA MSI of a specific MetricsOutput file
qmo_data$filename$dna_qc_metrics_msi
```

The available sections that can be accessed are the following: `header`, `notes`, `run_qc_metrics`, `analysis_status`, `dna_qc_metrics`, `dna_qc_metrics_snvtmb`, `dna_qc_metrics_msi`, `dna_qc_metrics_cnv`, `rna_qc_metrics`, `dna_expanded_metrics`, and `rna_expanded_metrics` 

*TSO500R(eader)* also provides convenience functions to access the data from all individual tables given in a `MetricsOutput` file:

```r
# get a data frame with all qc metrics on CNVs for all samples and the given recommended thresholds
read_dna_qc_metrics_cnv(qmo_data)

# get a data frame with the expanded dna qc metrics for all samples and the given recommended thresholds
read_dna_expanded_metrics(qmo_data)
```
These data frames provide the metrics for one sample per row.

Please refer to the [manual](https://github.com/apeltzer/TSO500R/tree/main/man) for a complete list of functions.

### Loading CombinedVariantOutput Data

To load all `CombinedVariantOutput.tsv` files of a DRAGEN analysis run in a specified folder, run the following command on the directory containing the files:

```r
# load all CombinedVariantOutput.tsv files of a LocalApp run in a specific folder
cvo_data_dir <- "/path/to/your/TSO500/analysis/output/data/"
cvo_data <- read_cvo_data(tso500_data_dir)
```

In case you want to read in data that resulted from a `LocalApp` run, you need to set `local_app=TRUE` due to the (slight) differences between the file formats:

```r
# load all CombinedVariantOutput.tsv files in a specific folder
cvo_data_dir <- "/path/to/your/TSO500/analysis/output/data/"
cvo_data <- read_cvo_data(tso500_data_dir, local_app=TRUE)
```

In case you want to read in data that resulted from a `ctdna` pipeline run, you need to set `ctdna=TRUE` due to the (slight) differences between the file formats:

```r
# load all CombinedVariantOutput.tsv files in a specific folder
cvo_data_dir <- "/path/to/your/TSO500/analysis/output/data/"
cvo_data <- read_cvo_data(tso500_data_dir, ctdna=TRUE)
```

This will produce a list of `combined.variant.output` objects, which can be further processed and contain all the information as given in the parsed files. Each section in the `CombinedVariantOutput` file can be accessed with standard list indexing. The individual `combined.variant.output` objects in the list can be accessed via the `PairID` as given in the `CombinedVariantOutput.tsv` file.

```r
# access analysis detail of sample1
cvo_data$sample1$sequencing_run_details

# access small variants of sample2
cvo_data$sample2$small_variants
```

The available sections are the following: `analysis_details`, `sequencing_run_details`, `tmb`, `msi`, `fusions`, `gene_amplifications`, `splice_variants`, and `small_variants`

#### Accessing Analysis Details

You can also use a convenience function to retrieve the analysis details for every samples the `CombinedVariantOutput.tsv` file has been parsed for:

```r 
# get data frame with analysis details for every sample in list of CVO objects 
get_analysis_details_df(cvo_data)
```

The same can be done for the sequencing run details using the corresponding function `get_sequencing_run_details_df`.

#### Accessing TMB and MSI Metrics

You can use the following function to extract TMB/MSI data from the list of `combined.variant.output` objects:

```r
# get TMB and MSI metrics from the CombinedVariantOutput.tsv files of all samples
tmb_msi_df <- get_metrics_df(cvo_data)
```

This will produce a dataframe with TMB/MSI metrics given in `CombinedVariantOutput.tsv` files and the corresponding sample identifiers. The metrics given as columns are the following: `total_tmb`, `coding_region_size_in_megabases`, `number_of_passing_eligible_variants`, `usable_msi_sites`, `total_msi_sites_unstable`, and `percent_unstable_msi_sites`

The TSO500 analysis using `LocalApp` or the DRAGEN pipeline also provides additional files with details on the TMB. We also provide parser functions for them.

```r
tmb_directory <- "/path/to/your/TSO500/analysis/output/data/"

# read in a batch of TMB_trace.tsv files into a data frame
tmb_trace_data <- read_tmb_trace_data(tmb_directory)

# read in a batch of *tmb.json files into a data frame
tmb_details <- read_tmb_details_data(tmb_directory)
```

#### Accessing Small Variants 

You can also extract all reported small variants across all of the `combined.variant.output` objects:

```r
# get small variants from the CombinedVariantOutput.tsv files of all samples
small_variant_df <- read_small_variants(cvo_data)
```

This will produce a dataframe with small variants given in `CombinedVariantOutput.tsv` files and the corresponding sample identifiers.

The same is possible for genetic rearrangements:

```r
# get data frame with amplifications from the CombinedVariantOutput.tsv files of all samples
read_gene_amplifications(cvo_data)

# get data frame with fusions from the CombinedVariantOutput.tsv files of all samples
read_fusions(cvo_data)

# get data frame with splice variants from the CombinedVariantOutput.tsv files of all samples
read_splice_variants(cvo_data)
```

### Integration of Data

We provide different functions to enrich the data frame containing small variant data information with other types of data.

```r 
# add information given in Illumina TMB trace table to small variant data frame
small_variants_with_tmb_info <- small_variant_df %>%
  add_tmb_variant_data(tmb_trace_data)

# adds amplifications to small variant data frame
small_variants_with_amps <- small_variant_df %>%
  add_amplification_data(read_gene_amplifications(cvo_data))
```

### Filtering and Processing

We also provide processing/filtering functions. If you want to do the following filtering steps, please use the function `process_small_variant_data` on your data frame with small variants.

- filters for variants that
  - have a consequence in a pre-defined list
  - are present with depth > 0
- extracts NP ID and protein information from the P Dot-notation column
- adds columns to faciliate addition of reference data

The following variant consequences are included:

- frameshift_variant
- inframe_deletion
- inframe_insertion
- missense_variant
- missense_variant:splice_region_variant
- splice_acceptor_variant
- splice_donor_variant
- splice_donor_variant:intron_variant
- start_lost
- stop_gained
- stop_gained:splice_region_variant
- stop_lost

```r
# filter small variant data frame as described above
small_variants_with_filtered <- small_variant_df %>%
  process_small_variant_data()
```

The package also provides filtering functions for filterting out or keeping variants with specified consequences:

```r
# filter out (remove) variants from data frame that have any of the specified consequences
# in this case remove all intronic variants
consequences <- c("intron_variant")
filter_consequences(small_variant_df, consequences, "consequence_s")

# keep variants from data frame that have any of the specified consequences
# in this case keep all missense variants
consequences <- c("missense_variant")
keep_consequences(small_variant_df, consequences, "consequence_s")
```

If your consequences are in a different column (other than `consequence_s`) please specify the name as a third argument.

You can also apply filters based on other properties such as the read depth of the determined variant:

```r
# keep all small variants that have a depth HIGHER than 100
small_variants_with_filtered <- small_variant_df %>%
  filter_depth(100)
```

The other filterting functions that can be applied on your small variant data frame if it has the corresponding columns are the following:

- filter_germline_db 
- filter_germline_proxi
- filter_for_cosmic_id
- filter_for_Included_in_TMB

### Plotting

*TSO500R(eader)* also provides some (basic) plotting functions. The following functions exist in the current version:

- plot allele frequency per variant consequence
- plot allele frequency kernel density estimate (KDE) per sample
- plot allele frequeny histogram per sample

The basic usage is as follows:

```r
plot_af_per_variant_consequence(small_variant_df)
```

You can also add commonly-used theme elements, if you prefer:

```r
plot_af_per_variant_consequence(small_variant_df) %>%
  add_common_theme_elements()
```

We also provide functionality for generating basic OncoPrint plots. If you need specific more complex annotations, we recommend to use the OncoPrint functionality directly.

```r
# prepare matrix based on small variant data frame
df_for_oncoprint <- prepare_dataframe_for_oncoprint(small_variant_df_filtered, "sample_id", "gene", "consequence_s")

# call function to generate OncoPrint plot
# alter_list: list with alter_graphics functions for variant types
# variant_colors: list of defined colors for variant types
# heatmap_legend: list containing legend specs (title, at, labels)
plot_onco_print(variant_matrix, "your favorite column title", alter_list, variant_colors, heatmap_legend)
```

Please refer to the [OncoPrint documentation](https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html)  for details.

### Other Convenience Functions

If we need to extract the counts of the different variant types from a list of `combined.variant.output` objects, you can use the following function which will return a data frame of counts per sample.

```r
# get count data frame on variant type for list of `combined.variant.output` objects
get_count_df <- function(cvo_data)
```

If you want to use your TSO500 data for generating an [OncoPrint](https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html) plot, you can use the following function to prepare the data accordingly.

```r
# make sure that your small variant data frame does not have empty fields/ NAs in the needed columns
# you do not have to use the same function as in this example
small_variant_df_filtered <- read_small_variants(cvo_data) %>%
  process_and_filter_small_variant_data()

# run function on small variant data frame
# if the needed columns have the same name as here you do not have to define them explicitly since these are the default values
# the columns are id (e.g. sample), gene, and variant type 
df_for_oncoprint <- prepare_dataframe_for_oncoprint(small_variant_df_filtered, "sample_id", "gene", "consequence_s")
```
