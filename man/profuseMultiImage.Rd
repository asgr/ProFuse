\name{profuseMultiImageFound2Fit}
\alias{profuseMultiImage}
\alias{profuseMultiImageFound2Fit}
\alias{profuseMultiImageDoFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProFound to ProFit Profile Fitting for Multiple Images
}
\description{
Functions to easily take users from multiple image pixel data all the way to fully fitted profiles. The key reason to use this function is when you have multiple images (i.e. exposures) of a galaxy in the same band. It is better to fit these simulataneously than try to stack and fit a single image. This is because stacking is prone to introducing pixel covariance and producing less well defined PSFs.
}
\usage{
profuseMultiImageFound2Fit(image_list, segim_list = NULL, mask_list = NULL,
  Ncomp = 2, loc = NULL, cutbox = NULL, psf_list = NULL, magzero = NULL,
  gain = NULL, resamp = NULL, tightcrop = TRUE, offset_list = NULL, ...)

profuseMultiImageDoFit(image_list, F2F = NULL, Ncomp = 2, psf_list = NULL,
  magzero = NULL, psf_dim = c(51, 51), plot = FALSE, seed = 666,
  optim_iters = 5, Niters = c(200, 200), NfinalMCMC = 1000, walltime = Inf,
  keepall = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_list}{
List; required, images to pass into \code{\link{profuseFound2Fit}} \option{image}.
}
  \item{F2F}{
List; output of \code{profuseMultiImageFound2Fit}. If this is provided then the \code{profuseMultiImageFound2Fit} step of \code{profuseMultiImageDoFit} is skipped and this object is used directly for inference. Note you may need to run \code{\link{profuseRegenPSF_MF2F}} on the \option{F2F} object if you are using across different session to regenerate the PSF convolution pointer.
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
  \item{tightcrop}{
Logical; optional, should the image be cut down to the tightest fitting box that fully contains the segmentation map? This might be a good idea if the parent image is very large (cuts down on memory allocation etc), but might be a bad idea if you want to know the exact position of your source with the frame, e.g. \option{xcen} and \option{ycen} in the \option{modellist} are relative to the new tightly cropped image not the original cropped image.
}
  \item{offset_list}{
List; optional, offsets to pass into \code{\link{profuseFound2Fit}} \option{psf}.
}
  \item{psf_dim}{
Integer vector; optional, the dimensions of the output PSF generated when fitting a star with \option{Ncomp} = 0.5. The PSF will be centrally located in the image. The dimensions will be forced to be odd (adding 1 if required) to ensure the PSF mode is at the centre of the middle pixel.
}
  \item{plot}{
Logical; optional, should a before and after modelling plot be created? This is useful when testing, but should be set to FALSE (default) for big runs.
}
  \item{seed}{
Integer scalar; random seed to start the \code{Highlander} function with.
}
  \item{optim_iters}{
Integer scalar; number of CMA / LD loops. See \code{\link{Highlander}}. The default of 5 works pretty well in practice.
}
  \item{Niters}{
Integer vector; number of iterations per CMA and LD respectively (must be length 2). See \code{\link{Highlander}}. The default of c(200,200) works pretty well in practice.
}
  \item{NfinalMCMC}{
Integer scalar; number of iterations to run for the final MCMC run. See \code{\link{Highlander}}.
}
  \item{walltime}{
Numeric scalar; the maximum allowed CPU time for \code{\link{Highlander}} in minutes. The function will stop with the best solution to date once the walltime has been exceeded (although it does not stop mid CMA or MCMC, so the walltime will usually be exceeded a bit).  
}
  \item{keepall}{
Logical; if FALSE (default) then does nothing. If TRUE then output objects \option{CMA_all} and \option{LD_all} will be concatenated lists of all CMA and LD iterations respectively.
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
