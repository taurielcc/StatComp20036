% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{kNN}
\alias{kNN}
\title{A KNN function using Rcpp}
\usage{
kNN(X, y, Z, n_classes, k)
}
\arguments{
\item{X}{matrix of training set cases}

\item{y}{label of true classification of training set}

\item{Z}{matrix of test set cases}

\item{n_classes}{number of classes in y}

\item{k}{number of neighbours considerd}
}
\value{
dist Euclidean distance matrix

results label of classification of test set
}
\description{
A KNN(K-Nearest Neighbor) classification for test set from training set using Rcpp. For each row of the test set, the k nearest (in Euclidean distance) training set vectors are found, and the classification is decided by majority vote.
}
\examples{
\dontrun{
X<-matrix(c(1:12),ncol=2)
y<-c(0,0,0,1,1,1)
Z<-matrix(c(2,6,7,11),ncol=2)
res<-KNN(X,y,Z,2,k=3)
}
}
