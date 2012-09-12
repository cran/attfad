\name{expressionStats}
\alias{expressionStats}

\title{Expression data analysis}

\description{
Calculates several statistics from an expression data matrix, a gene regulatory network, and, optionally, a vector of replicate groups.
}

\usage{
expressionStats(expression.data, network, exp.group = NA, permute.genes = FALSE)
}

\arguments{
  \item{expression.data}{
numeric matrix of expression values with microarray experiments as rows and genes as columns.
}
  \item{network}{
character matrix or data.frame with transcription factors in the first column and target genes in the second column.
}
  \item{exp.group}{
numeric vector containing seperate integers for each replicate group.
}
  \item{permute.genes}{
logical indicating if genes should be permuted in each experiment. Useful to compute background distributions. Default is false.
}
}

\value{
A list of statistics, most in the form of numeric vectors:
\item{Replicate.noise.distribution}{Pearson correlation coefficients of all pairs of experiments that belong to the same replicate group.}
\item{Absolute.intensity.distribution}{Simply the values from expression.data.}
\item{Range.of.gene.expression}{The expression range for each gene is estimated by taking the range from the 5\% quantile to the 95\% quantile and dividing by the median of these ranges across all genes.}
\item{Silhouette.coefficient.genes}{The distribution of silhouette coefficients at each node of the dendrogram resulting from hierarchical average linkage clustering of genes. High sihouette coefficients indicate small within-cluster-distances and large between-cluster-distances.}
\item{Silhouette.coefficient.microarrays}{The distribution of silhouette coefficients at each node of the dendrogram resulting from hierarchical average linkage clustering of microarrays. High sihouette coefficients indicate small within-cluster-distances and large between-cluster-distances.}
\item{TF.TG.correlation}{All pairwise Pearson correlation coefficients of transcription factors and corresponding target genes.}
\item{TF.TG.correlation.enrichment}{Indicates how the correlation coefficients of TF-TG pairs are distributed compared to the background of all gene pairs. For each of several bins between -1 and 1, the normalized number of all pairs falling into this bin is subtracted from the number of pairs which are in a TF-TG relation.} 
\item{TG.TG.correlation}{All pairwise Pearson correlation coefficients of target genes which are regulated by the same set of transcription factors.}
\item{TG.TG.correlation.enrichment}{Equivalent to TF.TG.correlation.enrichment.}
\item{TF.activity.distribution}{} 

}


\author{
Robert Maier
}


\seealso{
\code{plotComparison}
\code{plotAll}
\code{getAllScores}
}

\examples{
data(GNW.s.expression.data)
data(GNW.s.network)
data(GNW.s.exp.group)
results <- expressionStats(GNW.s.expression.data, GNW.s.network, GNW.s.exp.group)
str(results)
}

