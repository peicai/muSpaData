meta <- data.frame(
    stringsAsFactors = FALSE,
    Title = c(
      "Wei22_full",
      "Wei22_example"
      ),
    Description = c(paste(
        "Single-cell Stereo-seq spatial transcriptomics data with serial sections",
        "along the rostral-caudal axis of axolotl brain tissues,",
        "collected at 2 (3 sections), 5 (3 sections), 10 (3 sections), ",
        "15 (4 sections) and 20 (3 sections) days post injury (DPI)",
        "after removal of a portion of the dorsal pallium of",
        "the 11 cm length axolotl.",
        "Represented as a SpatialExperiment;",
        "derived from https://db.cngb.org/stomics/artista/download/."),
        paste(
          "A subset of the Wei22_full dataset,",
          "focusing on fewer regeneration stages.",
          "It includes 2 (2 sections), 10 (2 sections), ",
          "and 20 (2 sections) days post injury (DPI).")
    ),
    BiocVersion = "3.20",
    Genome = NA,
    SourceType = "rds",
    SourceUrl = "https://db.cngb.org/stomics/artista/download/",
    SourceVersion = "Aug 1 2021",
    Species = "Ambystoma mexicanum",
    TaxonomyId = "8296",
    Coordinate_1_based  = NA,
    DataProvider = "Spatial Transcript Omics DataBase (STOmics DB)",
    Maintainer = "Peiying Cai <peiying.cai@uzh.ch>",
    RDataClass = "SpatialExperiment",
    DispatchClass = "Rda",
    RDataPath = file.path(
      "muSpaData",
      "Data",
      c(
        "Wei22_full.rda",
        "Wei22_example.rda"
      )
    )
)

# write to .csv
write.csv(meta, file = "./inst/extdata/metadata.csv", row.names = FALSE)
