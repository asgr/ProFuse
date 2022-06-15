\name{profuseMultiBand}
\alias{profuseMultiBandFound2Fit}
\alias{profuseMultiBandDoFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProFound + ProFit + ProSpect Multi-Band Fitting
}
\description{
Functions to easily take users from multi-band image pixel data all the way to fully fitted physical profiles using \code{ProSpect}.
}
\usage{
profuseMultiBandFound2Fit(image_list, segim_list = NULL, segim_global = NULL,
  sky_list = NULL, skyRMS_list = NULL, parm_global = c("sersic.xcen1",
    "sersic.ycen1", "sersic.re1", "sersic.ang2", "sersic.axrat2"),
  Ncomp = 2, loc = NULL, cutbox = NULL, psf_list = NULL, magdiff = 2.5,
  magzero = NULL, gain = NULL, resamp = NULL, doprofit = NULL, sing_nser = 2,
  bulge_nser = 4, disk_nser = 1, sing_nser_fit = TRUE, bulge_nser_fit = FALSE,
  disk_nser_fit = FALSE, bulge_circ = TRUE, nser_upper = 5.3, star_rough = TRUE,
  fit_rough = FALSE, psf_dim = c(51, 51), star_con = 2, star_con_fit = TRUE,
  star_circ = TRUE, tightcrop = TRUE, wave = NULL, smooth.parm = NULL,
  parm_ProSpect = NULL, data_ProSpect = NULL, logged_ProSpect = NULL,
  intervals_ProSpect = NULL, autoclip = TRUE, roughpedestal = TRUE, ...)

profuseMultiBandDoFit(image_list, MF2F = NULL, parm_global = c("sersic.xcen1",
  "sersic.ycen1", "sersic.re1", "sersic.ang2", "sersic.axrat2"), Ncomp = 2,
  magzero = NULL, seed = 666, optim_iters = 5, Niters = c(200,200),
  NfinalMCMC = 1000, walltime = Inf, keepall = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_list}{
List; required for \code{profuseMultiBandFound2Fit} and for \code{profuseMultiBandDoFit} if \option{MF2F} is not provided. A list a numeric matrices (one image per target band). This list should exactly match the length of \option{segim_list}, \option{sky_list}, \option{skyRMS_list}, \option{psf_list}, \option{magzero} if these are supplied.
}
  \item{segim_list}{
List; optional, a list of integer matrix segmentation maps (one per target band). This list should exactly match the length of \option{image_list}. If not present \code{profoundProFound} will be run to estimate the stacked \option{image_list} segmentation map and use this instead.
}
  \item{segim_global}{
Integer matrix; optional, an input segmentation map to use on the stacked images to get rough initial guesses for the \code{ProFit} model. This will also be used to calculate bounds if \option{tightcrop} = TRUE. If this is not provided and \option{segim_list} is then a maximal \option{segim_global} will be computed internally using \code{\link{profuseSegimGlobal}}.
}
  \item{sky_list}{
List; optional, a list of numeric matrix sky images (one per target band). This list should exactly match the length of \option{image_list}. If not present \code{profoundProFound} will be run to estimate the per band sky.
}
  \item{skyRMS_list}{
List; optional, a list of numeric matrix sky RMS images (one per target band). This list should exactly match the length of \option{image_list}. If not present \code{profoundProFound} will be run to estimate the per band sky RMS.
}
  \item{MF2F}{
List; output of \code{profuseMultiBandFound2Fit}. If this is provided then the \code{profuseMultiBandFound2Fit} step of \code{profuseMultiBandDoFit} is skipped and this object is used directly for inference. Note you may need to run \code{\link{profuseRegenPSF_MF2F}} on the \option{MF2F} object if you are using across different session to regenerate the PSF convolution pointers.
}
  \item{Ncomp}{
Integer scalar; optional, number of components to fit to galaxy. Options for galaxies are "1" (single Sersic profile); "1.5" (disk Sersic profile with a central PSF to model the bulge); "2" (double Sersic profile, the default); and "3" (triple Sersic profile with one bulge and two [thin/thick] disks). In these cases an input \option{psf} should be supplied. To fit a star with a Moffat profile set Ncomp to "0.5". In this case no input \option{psf} should be supplied. To fit a star with a provided PSF set Ncomp to "0".
}
  \item{loc}{
Numeric vector; optional, the [X,Y] location of our target galaxy on the original \option{image} which will be targetted and/or cut out. The default is NULL, which targets the galaxy at the centre of the original \option{image}.
}
  \item{cutbox}{
Numeric vector; optional, the dimensions of the box to cut out from \option{image} centred on \option{loc}. This is ignored if \option{loc} is NULL (default). This can be useful if the target \option{image} is very large since a part of this process is to run \code{profoundProFound}, where we only really care about the solution near our object of interest.
}
  \item{parm_global}{
Vector; optional, a vector specifying the positions or names of parameters that should be considered global. These are still optimised over, but the optimised values are shared across all bands. This is usually appropriate for parameters related to xcen, ycen, axrat and ang, since these should be common between images in most cases. The default case is reasonably sensible for a two component fitting problem (the default). For convenience, if NULL then all parameters are treated as global.
}
  \item{psf_list}{
List; optional, a list of numeric matrices (one PSF per target band). This list should exactly match the length of \option{sky_list}, \option{skyRMS_list}, \option{image_list}, \option{magzero} if these are supplied. If not present \code{profitAllStarDoFit} will be run to estimate the per band PSF.
}
  \item{magdiff}{
Numeric scalar; optional, the maximum magnitude difference for modelling nearby objects (if \option{} = TRUE) that have a segment touching the object of interest. The default of 2.5 means objects have to have at least 10\% the flux of our target object to be modelled. Dropping this too low might mean a lot of time is spent modelling with ProFit for little gain in the fiddelity of the measurement.
}
  \item{magzero}{
Numeric vector; optional, the magnitude zero point for the corresponding \option{image_list}, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)). This vector should exactly match the length of \option{image_list}.
}
  \item{gain}{
Numeric vector; gain (in photo-electrons per ADU) for the corresponding \option{image_list}. This is only used to compute object shot-noise component of the flux error. If NULL then it is approximated from the image directly from Poisson statistics using \code{\link{profoundGainEst}}.
}
  \item{resamp}{
Numeric vector; optional, reampling factors for the corresponding \option{image_list}. If you have warped natively low resolution data (say, FIR) to match higher resolution optical data the resampling factor for the FIR will be larger than 1 (less than 1 is ignored) which modifies the measure sky RMS. In simple terms, by over-sampling lower resolution data and conserving flux you reduce the estimated RMS likelihood (your RMS per new pixel will be a factor \option{resamp} larger, but you have \option{resamp}^2 more pixels). It obviously does not make sense to change the likelihood just because you manipulate your data (in fact this process can only ever lose information).
}
  \item{doprofit}{
Numeric vector; optional should ProFit component structures be computed for the corresponding \option{image_list}? If Null or vector component is TRUE then this will calculate LL based on the whole image, if FALSE then the LL is only computed for flux within the relevant \option{segim_list} element. This makes sense when the source is very unresolved (typically UV/MIR/FIR) since there is very little structural information present, but we still might want to use the integrated flux properties of the source to constrain the total \code{ProSpect} SFH/ZH aspect of our model.
}
  \item{sing_nser}{
Numeric scalar; optional, start value (or fixed value if \option{sing_nser_fit} = FALSE) of the single Sersic profile index. Only relevant for \option{Ncomp} = 1 fits.
}
  \item{bulge_nser}{
Numeric scalar; optional, start value (or fixed value if \option{bulge_nser_fit} = FALSE) of the double Sersic profile bulge index. Only relevant for \option{Ncomp} = 2 fits.
}
  \item{disk_nser}{
Numeric scalar; optional, start value (or fixed value if \option{disk_nser_fit} = FALSE) of the double Sersic profile disk index. Only relevant for \option{Ncomp} = 2 or 1.5 fits.
}
  \item{sing_nser_fit}{
Logical; optional, should the single Sersic profile index be optimised during fitting? If FALSE (default) then is will be fixed at \option{sing_nser}. Only relevant for \option{Ncomp} = 1 fits.
}
  \item{bulge_nser_fit}{
Logical; optional, should the double Sersic profile bulge index be optimised during fitting? If FALSE (default) then is will be fixed at \option{bulge_nser}.  Only relevant for \option{Ncomp} = 2 fits.
}
  \item{disk_nser_fit}{
Logical; optional, should the double Sersic profile disk index be optimised during fitting? If FALSE then is will be fixed at \option{disk_nser}. Only relevant for \option{Ncomp} = 2 or 1.5 fits.
}
  \item{bulge_circ}{
Logical; optional, should the bulge be forced to be circular? If TRUE this means the \option{axrat} is fixed to be 1. Only relevant for \option{Ncomp} = 2 fits.
}
  \item{nser_upper}{
Numeric scalar; optional, upper limit of the sersic index value. Relevant for \option{Ncomp} = 1, 1.5, 2, or 3 fits.
}
  \item{star_rough}{
Logical; optional, specifies whether \code{profitAllStarDoFit} should be run with 'rough' profiles (default TRUE, since we just want to generate PSF images). This is only relevant if the \option{psf} has not been provided and we need to calculate ourself using \code{\link{profuseAllStarDoFit}}.
}
  \item{fit_rough}{
Logical; optional, specifies whether the main galaxy fitting using \code{profitDoFit} should be run with 'rough' profiles (default FALSE).
}
  \item{psf_dim}{
Integer vector; optional, the dimensions of the output PSF generated when fitting a star with \option{Ncomp} = 0.5. The PSF will be centrally located in the image. The dimensions will be forced to be odd (adding 1 if required) to ensure the PSF mode is at the centre of the middle pixel.
}
  \item{star_con}{
Numeric scalar; optional, start value (or fixed value if \option{star_con_fit} = FALSE) of the Moffat profile concentration index for star fitting. This is only relevant if the \option{psf} has not been provided and we need to calculate ourself using \code{\link{profuseAllStarDoFit}}.
}
  \item{star_con_fit}{
Logical; optional, should the Moffat profile concentration index be optimised during fitting? If FALSE then is will be fixed at \option{star_con}. This is only relevant if the \option{psf} has not been provided and we need to calculate ourself using \code{\link{profuseAllStarDoFit}}.
}
  \item{star_circ}{
Logical; optional, should the star be forced to be circular? If TRUE this means the \option{axrat} is fixed to be 1. This is only relevant if the \option{psf} has not been provided and we need to calculate ourself using \code{\link{profuseAllStarDoFit}}.
}
  \item{tightcrop}{
Logical; optional, should the image be cut down to the tightest fitting box that fully contains the segmentation map? This might be a good idea if the parent image is very large (cuts down on memory allocation etc), but might be a bad idea if you want to know the exact position of your source with the frame, e.g. \option{xcen} and \option{ycen} in the \option{modellist} are relative to the new tightly cropped image not the original cropped image.
}
  \item{wave}{
Numeric vector; optional. This is only relevant if you are providing smoothing functions via \option{smooth.parm}. If provided it should be the wavelength of each image in \option{image_list}, where \option{image_list} then needs to be in ascending (blue to red) wavelength order.
}
  \item{smooth.parm}{
List; optional. List containing named per band features (e.g. 'mag1', 'mag2', 're1', 're2') that you wish to smooth over to stop sharp discontinuities as a function of wavelength. For most purposes setting each of these to \code{\link{smooth.spline}} usually works pretty well, e.g. smooth.parm = list(mag1 = smooth.spline, re1 = smooth.spline).
}
  \item{parm_ProSpect}{
Numeric list or vector; the arguments to pass into \code{ProSpect} for fitting. The naming convention is the \code{ProSpect} argument '_N' where 'N' is the component number (1 for bulge, 2 for disk).
}
  \item{data_ProSpect}{
List; extra arguments (usually in the form of data) to pass into \code{ProSpect} for fitting. This will often be things like \option{massfunc} and \option{speclib}.
}
  \item{logged_ProSpect}{
Logical vector; optional; match length of \option{parm_ProSpect}, and elements should say which values should be raised to the power 10 when fitting. This is useful since scale parameters are better fit in log space, so \option{logged_ProSpect}=TRUE for those elements. If missing assumes all parameters are to be fit in linear space.
}
  \item{intervals_ProSpect}{
List, optional; has elements \option{lo} and \option{hi} which if provided must be the same length as \option{parm_ProSpect}. These specify the lower and upper limits in the same linear/log space as \option{parm_ProSpect}, i.e. set \option{lo} to -1 if you do not want to logged unit to go below 0.1.
}
 \item{autoclip}{
Logical; should ihe provided \option{image} by automatically clipped for bad pixels? If TRUE then the \option{image} median is subtracted from \option{image}, and pixels 10 times below the 0.1\% percentile and 10 times above the 99.9\% percentile are masked. if FALSE nothing is done (users should be careful to set bad pixels to NA in the \option{image}).
}
  \item{roughpedestal}{
Logical; when the initial "rough sky" is computed, should only a pedestal (based on the median of the sky map) be used for the sky? This is a good option if the image is known to contain a *very* large (many times the box size) galaxy that might otherwise be over subtracted by the initial rough sky map. This is the only \code{\link{profoundProFound}} option that is changed from the default, since our requirement is that all \option{image} inputs are already sky subracted properly.
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
For \code{*DoFit} functions dots are passed into the similarly titled \code{*Found2Fit} function. For \code{*Found2Fit} functions, further arguments are passed to \code{profoundProFound} for making the segmentation maps. \option{tolerance} and \option{SBdilate} are two options that are usually worth experimenting with the get ideal segmentation results for fitting purposes.
}
}
\details{
\code{profitMultiBandFound2Fit} is a high-level and simplified multi-band focussed version of \code{profitFound2Fit}. For many users who wish to fit across multiple bands of data this is probably the easiest entry point unless you really need to lower level control provided by \code{profitFound2Fit}.

\code{profitMultiBandDoFit} is a high-level and simplified multi-band focussed version of \code{profitDoFit}. For many users who wish to fit across multiple bands of data this is probably the easiest entry point unless you really need to lower level control provided by \code{profitFound2Fit}.

The default combination of \option{optim_iters}, \option{Niters} and \option{NfinalMCMC} seems to work well in practice, and was used for the initial low redshift (z<=0.06) GAMA sample.
}
\value{
The output of \code{profitFound2Fit} and \code{profitAllStarFound2Fit} is:

\item{profound}{Output from the initial \code{profoundProFound} run. See \code{profoundProFound}.}
\item{Data}{The input Data structure created by \code{\link{profitSetupData}}.}

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
\references{
None yet (paper submitted).
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profuseFound2Fit}}, \code{\link{profuseDoFit}}, \code{\link{profuseAllStarFound2Fit}}, \code{\link{profuseAllStarDoFit}}, \code{\link{profitSetupData}}, \code{\link{profitLikeModel}},
\code{\link{profuseSegimGlobal}}
}
\examples{
#See vignettes for using profuseMultiBand functions.
}
\concept{ fit }
