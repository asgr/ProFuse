\name{profuseAllStar}
\alias{profuseAllStarFound2Fit}
\alias{profuseAllStarDoFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Automated Star Fitting
}
\description{
Functions to identify stars from \code{ProFound} and fit them with \code{ProFit} using \code{Highlander}.
}
\usage{
profuseAllStarFound2Fit(image, sigma = NULL, locs = NULL, segim = NULL, mask = NULL,
  magzero = 0, resamp = NULL, psf_dim = c(51, 51), star_con = 2, star_con_fit = TRUE,
  star_circ = TRUE, rough = TRUE, star_dom_mag = NULL, Nstar = 4, autoclip = TRUE,
  roughpedestal = TRUE, ...)

profuseAllStarDoFit(image, sigma = NULL, locs = NULL, magzero = 0, psf_dim = c(51,51),
  rough = TRUE, plot = FALSE, seed = 666, optim_iters = 2, Niters = c(100,100),
  NfinalMCMC = 100, walltime = Inf, keepall = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the galaxy image we want to fit a galaxy of PSF/star model to. For convenience this can also be an Rfits_pointer object for very large on disk images. Critically, this input \option{image} needs to be well background subtracted (this is not done again internally).
}
  \item{sigma}{
Numeric matrix; optional, the sigma map of the \option{image} data we want to fit a model to. If not provided this will be estimated by running \code{profoundProFound} (to find the sky and skyRMS), \code{profoundGainEst} (to find the gain) and \code{profoundMakeSigma} to make a decent effort at creating an appropriate sigma map. Users should usually try to provide this themselves for critical fitting, but the automatic sigma map is usually good enough for testing purposes. For convenience this can also be an Rfits_pointer object for very large on disk sigma maps.
}
  \item{segim}{
Integer matrix; optional, users can pass in their own segmentation maps, which should be pixel matched to the input \option{image}. This will be used instead of the internal \code{ProFound} run, although the segmentation statistics will still be recomputed.
}
 \item{mask}{
Boolean matrix or integer scalar (1,0); optional, parts of the image to mask out (i.e. ignore). If a matrix is provided, this matrix *must* be the same dimensions as \option{image} where 1 means mask out and 0 means use for analysis. If a scalar is provided it indicates the exact \option{image} values that should be treated as masked (e.g. by setting masked pixels to 0 or -999). The latter achieves the same effect as setting masked \option{image} pixels to NA, but allows for the fact not all programs can produce \code{R} legal NA values.
}
  \item{locs}{
Numeric matrix; optional, the [X,Y] locations of target stars on the original \option{image} for \code{AllStar} fitting to create out PSF model. The idea is we pass in a list of multiple potential stars, and internally this is limited to just those star segments that do not touch another source or hit the \option{image} border. If this is not provided then \code{profitAllStarFound2Fit} will attempt to guess which objects might be good stars for modelling (but this should not be relied on generally).
}
  \item{magzero}{
Numeric scalar/vector; optional, the magnitude zero point, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)). For \option{profoundMultiBandFound2Fit} this vector should exactly match the length of \option{sky_list}, \option{skyRMS_list}, \option{image_list}, \option{psf_list} if these are supplied.
}
  \item{resamp}{
Numeric scalar; optional, reampling factors for the \option{image}. If you have warped natively low resolution data (say, FIR) to match higher resolution optical data the resampling factor for the FIR will be larger than 1 (less than 1 is ignored) which modifies the measure sky RMS. In simple terms, by over-sampling lower resolution data and conserving flux you reduce the estimated RMS likelihood (your RMS per new pixel will be a factor \option{resamp} larger, but you have \option{resamp}^2 fewer pixels).
}
  \item{star_con}{
Numeric scalar; optional, start value (or fixed value if \option{star_con_fit} = FALSE) of the Moffat profile concentration index for star fitting.
}
  \item{star_con_fit}{
Logical; optional, should the Moffat profile concentration index be optimised during fitting? If FALSE then is will be fixed at \option{star_con}.
}
  \item{star_circ}{
Logical; optional, should the star be forced to be circular? If TRUE this means the \option{axrat} is fixed to be 1.
}
  \item{psf_dim}{
Integer vector; optional, the dimensions of the output PSF generated when fitting a star with \option{Ncomp} = 0.5. The PSF will be centrally located in the image. The dimensions will be forced to be odd (adding 1 if required) to ensure the PSF mode is at the centre of the middle pixel.
}
  \item{rough}{
Logical; optional, should ProFit be run using an approximate fit? This is less accurate, and potentially much less for compact sources, but it also runs much faster. Potentially if the target object is very extended then this will be accurate enough and speed up fitting by a factor of a few usually. It is also useful for quick testing. For publishable results this should usually be left as the default of FALSE. For \code{AllStar} setting this to FALSE is perhaps a good option if you do no care about the model parameters themselves (e.g. FWHM) and just want a decent PSF image.
}
  \item{star_dom_mag}{
Numeric scalar; optional, specify at what magnitude (and brighter) stars start to dominate source counts. This helps the code to guess which objects are stars if \option{locs} is not provided. If not provided then it will just use the median magitude of all sources (which is rarely correct, but at least it sits in the domain of the data). Sensible values are usually between mag 20-22 (for AB mag and optical bands).
}
  \item{Nstar}{
Integer scalar; optional, maximum number of fitting stars to generate in \code{profitAllStarFound2Fit}, and therefore to fit in \code{profitAllStarDoFit}.
}
  \item{autoclip}{
Logical; should ihe provided \option{image} by automatically clipped for bad pixels? If TRUE then the \option{image} median is subtracted from \option{image}, and pixels 10 times below the 0.1\% percentile and 10 times above the 99.9\% percentile are masked. if FALSE nothing is done (users should be careful to set bad pixels to NA in the \option{image}).
}
  \item{roughpedestal}{
Logical; when the initial "rough sky" is computed, should only a pedestal (based on the median of the sky map) be used for the sky? This is a good option if the image is known to contain a *very* large (many times the box size) galaxy that might otherwise be over subtracted by the initial rough sky map. This is the only \code{\link{profoundProFound}} option that is changed from the default, since our requirement is that all \option{image} inputs are already sky subracted properly.
}
  \item{plot}{
Logical; optional, should a before and after modelling plot be created? This is useful when testing, but should be set to FALSE (default) for big runs.
}
  \item{seed}{
Integer scalar; random seed to start the \code{Highlander} function with.
}
  \item{optim_iters}{
Integer scalar; number of CMA / LD loops. See \code{\link{Highlander}}. The default of 2 works pretty well in practice for stars.
}
  \item{Niters}{
Integer vector; number of iterations per CMA and LD respectively (must be length 2). See \code{\link{Highlander}}. The default of c(100,100) works pretty well in practice for stars.
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
For \code{*DoFit} functions dots are passed into the similarly titled \code{*Found2Fit} function. For \code{*Found2Fit} functions, further arguments are passed to \code{profoundProFound} for making the segmentation maps. \option{tolerance} and \option{SBdilate} are two options that are usually worth experimenting with the get ideal segmentation results for fitting purposes.
}
}
\details{
\code{profitAllStarFound2Fit} is a high-level and simplified PSF focussed version of \code{profitFound2Fit}. For many users who wish to fit stars to find the PSF this is probably the easiest entry point unless you really need to lower level control provided by \code{profitFound2Fit}.

\code{profitAllStarDoFit} is a high-level and simplified PSF focussed version of \code{profitDoFit}. For many users who wish to fit stars to find the PSF,this is probably the easiest entry point unless you really need to lower level control provided by \code{profitDoFit}.
}
\value{
The output of \code{profitAllStarFound2Fit} is:

\item{profound}{Output from the initial \code{profoundProFound} run. See \code{profoundProFound}.}
\item{Data}{The input Data structure created by \code{\link{profitSetupData}}.}
\item{star_segID}{The segID of the selected stars, which can be referred back to the \option{profound} output.}

The output of \code{profitDoFit}, \code{profitAllStarDoFit} and \code{profitMultiBandDoFit} is subset of:

\item{parm}{Best fit raw parameters being optimised (so various terms will be logged etc). See \code{Highlander}.}
\item{parm_smooth}{\code{profitMultiBandDoFit} only. Best fit smoothed parameters being optimised (so various terms will be logged etc). This is only relevant if \option{smooth.parm} smoothing functions are being supplied to create continuity between multi-band data. These are the parameters that \code{ProFit} sees internally when constructing model images.}
\item{LP}{The log posterior of the best fit. See \code{Highlander}.}
\item{diff}{Numeric scalar; LP difference between current best LP and last best LP. See \code{Highlander}.}
\item{best}{Character scalar; optimisation type of the best solution found, one of 'CMA' / 'LD_Median' / 'LD_Mean'. See \code{Highlander}.}
\item{iteration}{Integer scalar; iteration number of best solution. See \code{Highlander}.}
\item{CMA_last}{List; output of last CMA optimisation. See \code{Highlander}.}
\item{LD_last}{List; output of last MCMC optimisation. See \code{Highlander}.}
\item{call}{List; the full call made to Highlander. See \code{Highlander}.}
\item{date}{Date item; the date and time of the call. See \code{Highlander}.}
\item{time}{Numeric scalar; run time in minutes. See \code{Highlander}.}
\item{ProFit.version}{Version of \code{ProFit} used.}
\item{ProFound.version}{Version of \code{ProFound} used.}
\item{Highlander.version}{Version of \code{Highlander} used.}
\item{R.version}{Version of \code{R} used.}
\item{profound}{Output from the initial \code{profoundProFound} run. See \code{profoundProFound}.}
\item{Data}{The input Data structure created by \code{\link{profitSetupData}}.}
\item{initmodel}{The full initial model (all terms are in linear space).}
\item{finalmodel}{The full final model (all terms are in linear space).}
\item{error}{Fitting errors for all terms in \option{parm}. The errors are in the same raw units as \option{parm} (i.e. for logged units the errors are in dex).}
\item{psf}{The model Moffat image of the PSF estimated from fitting a target star with \option{Ncomp} = 0.5 (only present in this scenario) or multiple stars when using \code{profitAllStarDoFIt}. This will be created with dimensions \option{psf_dim}. The sum of the image is forced to exactly equal 1, even if the Moffat profile extends outside of the image domain. This at least ensures flux conservation during any convolution routines.}
\item{psf_modellist}{The model list for the PSF created in \option{psf}.}
\item{psf_fluxcheck}{The sum of the original PSF image flux. If far from 1 (say less than 0.95) then it implies the \option{psf_dim} is too small. If very close to 1 (say less than 0.999) then \option{psf_dim} is probably larger than it needs to be.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profuseFound2Fit}}, \code{\link{profuseDoFit}}, \code{\link{profuseMultiBandFound2Fit}}, \code{\link{profuseMultiBandDoFit}}, \code{\link{profitSetupData}}, \code{\link{profitLikeModel}}
}
\examples{
library(Rfits)
library(ProFit)

image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))$imDat

psf_modellist=list(
  moffat=list(
    xcen=75/2,
    ycen=75/2,
    mag=0,
    fwhm=3.8,
    con=2.04,
    ang=0,
    axrat=1,
    box=0
  )
)
psf = profitMakeModel(modellist=psf_modellist, dim=c(75,75))$z

\dontrun{

#Example single PSF/star fit:

star = profuseDoFit(image, loc=c(149,73), Ncomp=0.5, plot=TRUE, magzero=30, SBdilate=2)
magimage(star$psf)

#Example All PSF/star fit:

profound = profoundProFound(image, magzero=30)
sel = which(profound$segstats$mag < 23 & profound$segstats$R50 < 3.5)

#We pass in 21 possible star positions, choosing the best 4 by default:

stars_man = profuseAllStarDoFit(image, loc=profound$segstats[sel,c('xcen','ycen')],
  plot=TRUE, magzero=30, SBdilate=2)
magimage(stars_man$psf)

#We can run again, letting it figure out automatically which objects are stars:

stars_auto = profuseAllStarDoFit(image, plot=TRUE, magzero=30, SBdilate=2)
magimage(stars_auto$psf)
}
}
\concept{ fit }
