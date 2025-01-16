\name{Wei22_full}
\title{\code{Wei22_full} dataset}
\alias{Wei22_full}

\description{
\itemize{
\item{Single-cell Stereo-seq spatial transcriptomics data includes axolotl brain tissues
collected from multiple sections across five regeneration stages 
(16 samples in total)}
\item{The dataset contains 13890 features across 147432 cells, 
and preprocessing has been applied.}
\item{The original data is deposited in the Spatial Transcript Omics DataBase (STOmics DB)
under dataset id STDS0000056.}}}

\usage{Wei22_full(metadata = FALSE)}

\arguments{
\item{metadata}{
    A \code{logical} value indicating whether to return only 
    the ExperimentHub metadata, which describes the overall dataset, 
    or to load the entire dataset. Defaults to \code{FALSE}.}}

\details{
\describe{
    \item{link to reference:}{
        https://www.science.org/doi/10.1126/science.abp9444}
    \item{link to raw data:}{
        https://db.cngb.org/stomics/artista/download/}}}

\value{an object of class \code{SpatialExperiment}.}

\examples{
# load metadata only
Wei22_full(metadata = TRUE)

# load SpatialExperiment
Wei22_full(metadata = FALSE)
}

\author{Peiying Cai \email{peiying.cai@uzh.ch}} 

\references{
Wei et al. (2022). 
"Single-cell Stereo-seq reveals induced progenitor cells involved in axolotl brain regeneration", 
\emph{Science}, \bold{377}(6610).
}
