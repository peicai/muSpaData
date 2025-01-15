\name{Wei22_example}
\title{\code{Wei22_example} dataset}
\alias{Wei22}
\alias{Wei22_example}

\description{
\itemize{
\item{A subset of the Wei22_full dataset, focusing on fewer regeneration stages.}
\item{It includes 2 (2 sections), 10 (2 sections), and 20 (2 sections) days 
post injury (DPI).}
\item{The original data is deposited in the Spatial Transcript Omics DataBase (STOmics DB)
under dataset id STDS0000056.}}}

\usage{Wei22_example(metadata = FALSE)}

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
Wei22_example(metadata = TRUE)

# load SpatialExperiment
Wei22_example(metadata = FALSE)
}

\author{Peiying Cai \email{peiying.cai@uzh.ch}} 

\references{
Wei et al. (2022). 
"Single-cell Stereo-seq reveals induced progenitor cells involved in axolotl brain regeneration", 
\emph{Science}, \bold{377}(6610).
}
