\name{parse_SNP}
\alias{parse_SNP}
\title{parses LD file by filtering for SNPs that meet a cutoff value and detecting linkage}
\usage{
  parse_SNP(all_data, LD, gff_file, window_size, r_squared_cutoff)
}
\description{
  Parses LD file by filtering for SNPs that meet a cutoff value and detecting linkage and then detecting genes near the SNPs based on window size
}
\examples{
  parse_SNP(all_data, LD, "../example/example.gff.gz", 1000, 0.8)
}
