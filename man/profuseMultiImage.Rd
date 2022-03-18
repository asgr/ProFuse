\name{profuseMultiImageFound2Fit}
\alias{profuseMultiImageFound2Fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProFound to ProFit Profile Fitting for Multiple Images
}
\description{
Functions to easily take users from multiple image pixel data all the way to fully fitted profiles. The key reason to use this function is when you have multiple images (i.e. exposures) of a galaxy in the same band. It is better to fit these simulataneously than try to stack and fit a single image. This is because stacking is prone to introducing pixel covariance and producing less well defined PSFs.
}
\usage{
profuseMultiImageFound2Fit(image_list, segim_list = NULL, mask_list = NULL,
  Ncomp = 2, loc = NULL, cutbox = NULL, psf_list = NULL, magzero = NULL, gain = NULL,
  resamp = NULL, offset_list = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_list}{
List; required, images to pass into \code{\link{profuseFound2Fit}} \option{image}.
}
  \item{segim_list}{
List; required, segims to pass into \code{\link{profuseFound2Fit}} \option{segim}.
}
  \item{mask_list}{
List; optional, masks to pass into \code{\link{profuseFound2Fit}} \option{mask}.
}
  \item{Ncomp}{
Integer scalar; optional, number of components to fit to galaxy. Options for galaxies are "1" (single Sersic profile); "1.5" (disk Sersic profile with a central PSF to model the bulge); "2" (double Sersic profile, the default); and "3" (triple Sersic profile with one bulge and two [thin/thick] disks). In these cases an input \option{psf} should be supplied. To fit a star with a Moffat profile set Ncomp to "0.5". In this case no input \option{psf} should be supplied. To fit a star with a provided PSF set Ncomp to "0".
}
  \item{loc}{
Numeric vector; optional, the [X,Y] location which will be targetted and/or cut out. The default is NULL, which targets the galaxy at the centre of the first \option{image_list}.
}
  \item{cutbox}{
Numeric vector; optional, the dimensions of the box to cut out from \option{image} centred on \option{loc}. This is ignored if \option{loc} is NULL (default). This can be useful if the target \option{image} is very large since a part of this process is to run \code{profoundProFound}, where we only really care about the solution near our object of interest.
}
  \item{psf_list}{
List; optional, psfs to pass into \code{\link{profuseFound2Fit}} \option{psf}.
}
  \item{magzero}{
Numeric vector; optional, magnitude zeropoints to pass into \code{\link{profuseFound2Fit}} \option{magzero}. This vector should exactly match the length of \option{image_list}.
}
  \item{gain}{
Numeric vector; optional, gains to pass into \code{\link{profuseFound2Fit}} \option{gain}. This vector should exactly match the length of \option{image_list}.
}
  \item{resamp}{
Numeric scalar; optional, resampling factor to pass into \code{\link{profuseFound2Fit}} \option{resamp}.
}
  \item{offset_list}{
List; optional, offsets to pass into \code{\link{profuseFound2Fit}} \option{psf}.
}
  \item{\dots}{
Other arguments to pass into \code{\link{profuseFound2Fit}}.
}
}
\details{
\code{profitMultiImageFound2Fit} is a high-level and simplified multi-image focussed version of \code{profitFound2Fit}. For many users who wish to fit multiple images in the same band this is probably the easiest entry point unless you really need to lower level control provided by \code{profitFound2Fit}.
}
\value{
A list of length the number of images provided by \option{image_list}, where each list component is \option{Data} structure created by \code{\link{profitSetupData}}.
}
\references{
None yet (paper submitted).
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profuseFound2Fit}}
}
\examples{
library(Rfits)
library(ProFit)

\dontrun{
im1 = Rfits_read_image(system.file("extdata/Offset/Zim1.fits",package="ProFit"),
  ext=2)$imDat
im2 = Rfits_read_image(system.file("extdata/Offset/Zim2.fits",package="ProFit"),
  ext=2)$imDat

F2F = profuseMultiImageFound2Fit(image_list = list(im1, im2),
  offset_list = list(c(0,0), c(-23.112,-9.697)), magzero = 30, tolerance=20)

profitLikeModel(F2F[[1]]$init, Data=F2F, makeplots=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
