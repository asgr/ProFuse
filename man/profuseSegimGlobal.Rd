\name{profuseSegimGlobal}
\alias{profuseSegimGlobal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Global Segmentation Map
}
\description{
Based on a list, this will make a maximal global segmentation map. This means the minimal cutout will definitely contain the largest segment in any of the provided pre-band segmentation maps in segim_list. This is mostly used internally by \code{\link{profuseMultiBandFound2Fit}}, but uses can use this more directly if they have a complicated use-case.
}
\usage{
profuseSegimGlobal(segim_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim_list}{
List; required. The per-band segmentation map to maximally combine.
}
}
\details{
Uses \code{imager} and \code{foreach}, so these must be installed.
}
\value{
Returns an integer matrix of the maximal segmentation map.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profuseMultiBandFound2Fit}}
}
\examples{
##Nothing yet.
}
\concept{segmentation}

