\name{profuseRegenPSF}
\alias{profuseRegenPSF}
\alias{profuseRegenPSF_F2F}
\alias{profuseRegenPSF_MF2F}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regenerate PSF Convolvers
}
\description{
Simple functions to rebuild the session dependent convolver kernels should the Found2Fit objects be loaded into a new session.
}
\usage{
profuseRegenPSF_F2F(F2F)
profuseRegenPSF_MF2F(MF2F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{F2F}{
List; output of \code{\link{profuseFound2Fit}}.
}
  \item{MF2F}{
List; output of \code{\link{profuseMultiBandFound2Fit}}.
}
}
\details{
Not much to say here.
}
\value{
Returns a modified (with convolvers added) version \option{F2F} or \option{MF2F} as relevant.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profuseFound2Fit}}, \code{\link{profuseMultiBandFound2Fit}}
}
\examples{
# In rough code you would have
# F2F = profuseFound2Fit(args here)
# Then in a new session you would run:
# F2F = profuseRegenPSF_F2F(F2F)
# Or
# MF2F = profuseMultiBandFound2Fit(args here)
# MF2F = profuseRegenPSF_MF2F(MF2F)
}
