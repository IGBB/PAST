\name{merge_data}
\alias{merge_data}
\title{Merge statistics and effects data}
\usage{
merge("/path/to/statistics", "/path/to/effects")
}
\description{
Merges statistics data and effects data, filtering out non-biallelic data.
}
\examples{
merge_data("../example/GrainColorGLMStatsHiMAF.txt", "../example/GrainColorGLMAllEffHiMAF.txt")
}
