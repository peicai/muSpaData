# muSpaData
Multi-sample multi-group spatially resolved transcriptomic data

The `muSpaData` package includes datasets for use in the `DESpace`
package's examples and vignettes. It provides access to a publicly available Stereo-seq spatial dataset with complex experimental designs. This dataset, containing multiple samples (e.g., serial sections) measured under various experimental conditions (e.g., time points), is formatted as `SpatialExperiment` (SPE) Bioconductor objects. All provided SPEs include  processed data and meta data obtained from the original data source.

Currently available datasets:

- `Wei22_full`:  
 A comprehensive spatial transcriptomics dataset at single-cell resolution, capturing axolotl telencephalon regeneration stages. The dataset includes axolotl brain tissues collected from multiple sections across various regeneration stages: 2 (3 sections), 5 (3 sections), 10 (3 sections), 15 (4 sections), and 20 (3 sections) days post injury (DPI).  [[source](https://db.cngb.org/stomics/artista/), [reference](https://doi.org/10.1126/science.abp9444)]

- `Wei22_example`:  
  A subset of the Wei22_full dataset, focusing on fewer regeneration stages. This dataset includes 2 DPI (2 sections), 10 DPI (2 sections), and 20 DPI (2 sections) brain tissues. [[source](https://db.cngb.org/stomics/artista/), [reference](https://doi.org/10.1126/science.abp9444)]
