\name{HSE}
\alias{HSE}
\title{Homogeneity of Stratum Effects Test}
\description{
To test the homogeneity of stratum effects in stratified paired binary data.
}
\usage{
HSE(Tab)
}
\arguments{
  \item{Tab}{the table combining k 2*2 tables by rows, k is 1,2,...,n}
}
\value{
  \item{T}{the statistic following the chi-squared distribution}
  \item{p}{the p value, if p<0.05, there was a statistically significant difference among the stratum}
}
\note{
Please feel free to contact us, if you have any advice and find any bug!

Reference:

1. Yan D. Zhao, Dewi Rahardja, De-Hui Wang & Haili Shen (2014) Testing Homogeneity of Stratum Effects in Stratified Paired Binary Data, Journal of Biopharmaceutical Statistics, 24:3, 600-607, DOI: 10.1080/10543406.2014.888440

Update:

Version 0.1.0: The first version.

}
\examples{
table1 <- matrix(c(30,4,3,20),ncol=2)
table2 <- matrix(c(25,2,3,20),ncol=2)
table <- rbind(table1,table2)
HSE(table)
}
