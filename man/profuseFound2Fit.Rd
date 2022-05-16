\name{profuseFound2Fit}
\alias{profuseFound2Fit}
\alias{profuseDoFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProFound to ProFit Profile Fitting for Single Images
}
\description{
Functions to easily take users from single image pixel data all the way to fully fitted profiles.
}
\usage{
profuseFound2Fit(image, sigma = NULL, segim = NULL, mask = NULL, Ncomp = 2,
  loc = NULL, cutbox = NULL, psf = NULL, magdiff = 2.5, magzero = 0,
  gain = NULL, resamp = NULL, loc_use = FALSE, loc_fit = TRUE,
  mag_fit = TRUE, sing_nser = 2, bulge_nser = 4, disk_nser = 1,
  sing_nser_fit = TRUE, bulge_nser_fit = FALSE, disk_nser_fit = FALSE,
  bulge_circ = TRUE, nser_upper=5.3, star_rough = TRUE, fit_rough = FALSE,
  psf_dim = c(51, 51), star_con = 2, star_con_fit = TRUE, star_circ = TRUE,
  offset = NULL, tightcrop = TRUE, deblend_extra = FALSE, fit_extra = FALSE,
  pos_delta = 10, autoclip = TRUE, roughpedestal = TRUE, ...)

profuseDoFit(image, F2F = NULL, Ncomp = 2, psf = NULL, magzero = NULL,
  psf_dim = c(51,51), plot = FALSE, seed = 666, optim_iters = 5,
  Niters = c(200,200), NfinalMCMC = 1000, walltime = Inf, keepall = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required for \code{profuseFound2Fit} and for \code{profuseDoFit} if \option{F2F} is not provided. The galaxy image we want to fit a galaxy of PSF/star model to. For convenience this can also be an Rfits_pointer object for very large on disk images. Critically, this input \option{image} needs to be well background subtracted (this is not done again internally).
}
  \item{F2F}{
List; output of \code{profuseFound2Fit}. If this is provided then the \code{profuseFound2Fit} step of \code{profuseDoFit} is skipped and this object is used directly for inference. Note you may need to run \code{\link{profuseRegenPSF_F2F}} on the \option{F2F} object if you are using across different session to regenerate the PSF convolution pointer.
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
  \item{Ncomp}{
Integer scalar; optional, number of components to fit to galaxy. Options for galaxies are "1" (single Sersic profile); "1.5" (disk Sersic profile with a central PSF to model the bulge); "2" (double Sersic profile, the default); and "3" (triple Sersic profile with one bulge and two [thin/thick] disks). In these cases an input \option{psf} should be supplied. To fit a star with a Moffat profile set Ncomp to "0.5". In this case no input \option{psf} should be supplied. To fit a star with a provided PSF set Ncomp to "0".
}
  \item{loc}{
Numeric vector; optional, the [X,Y] location of our target galaxy on the original \option{image} which will be targetted and/or cut out. The default is NULL, which targets the galaxy at the centre of the original \option{image}.
}
  \item{cutbox}{
Numeric vector; optional, the dimensions of the box to cut out from \option{image} centred on \option{loc}. This is ignored if \option{loc} is NULL (default). This can be useful if the target \option{image} is very large since a part of this process is to run \code{profoundProFound}, where we only really care about the solution near our object of interest.
}
  \item{psf}{
Image matrix; required, an empirical point spread function (PSF) image matrix that ProFit will use to convolve the image
}
  \item{magdiff}{
Numeric scalar; optional, the maximum magnitude difference for modelling nearby objects (if \option{} = TRUE) that have a segment touching the object of interest. The default of 2.5 means objects have to have at least 10\% the flux of our target object to be modelled. Dropping this too low might mean a lot of time is spent modelling with ProFit for little gain in the fiddelity of the measurement.
}
  \item{magzero}{
Numeric scalar; optional, the magnitude zero points for the \option{image}, where values become scaled by the standard scale=10^(-0.4*(mag-magzero)).
}
  \item{gain}{
Numeric scalar; optional, gain (in photo-electrons per ADU) for the \option{image}. This is only used to compute object shot-noise component of the flux error. If NULL then it is approximated from the image directly from Poisson statistics using \code{\link{profoundGainEst}}.
}
  \item{resamp}{
Numeric scalar; optional, reampling factors for the \option{image}. If you have warped natively low resolution data (say, FIR) to match higher resolution optical data the resampling factor for the FIR will be larger than 1 (less than 1 is ignored) which modifies the measure sky RMS. In simple terms, by over-sampling lower resolution data and conserving flux you reduce the estimated RMS likelihood (your RMS per new pixel will be a factor \option{resamp} larger, but you have \option{resamp}^2 fewer pixels).
}
  \item{loc_use}{
Logical; optional, should the provided \option{loc} be used to represent the starting [X,Y] coordinates of the model? If FALSE (default) then the \option{xcen} and \option{ycen} provided by \code{profoundProFound} is used instead. If users are very certain that the centre of the object should be precisely at \option{loc} then \option{loc_use} should be TRUE.
}
  \item{loc_fit}{
Logical; optional, should the [X,Y] location of optimised as part of the fitting process (TRUE) or left at the starting position (FALSE). What the starting position means will change depending on the setting of \option{loc_use}.
}
  \item{mag_fit}{
Logical; optional, should the magnitudes optimised as part of the fitting process (TRUE) or left at the starting position (FALSE). This is only really here for \code{ProSpect} SED magnitude modelling in \code{profitMultiBandFound2Fit}, so should generally not be adjusted by users manually.
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
  \item{offset}{
Numeric vector; optional. Offset argument to be passed to \code{\link{profitSetupData}} (see description there).
}
  \item{tightcrop}{
Logical; optional, should the image be cut down to the tightest fitting box that fully contains the segmentation map? This might be a good idea if the parent image is very large (cuts down on memory allocation etc), but might be a bad idea if you want to know the exact position of your source with the frame, e.g. \option{xcen} and \option{ycen} in the \option{modellist} are relative to the new tightly cropped image not the original cropped image.
}
  \item{deblend_extra}{
Logical; optional, if TRUE (and \option{fit_extra} is FALSE) then the image will be re-scaled by the weight map relevant to each segment, but only for pixels within that segment. This means that neighbouring object flux is approximately removed within each segment of a group. This obviously does not conserve flux, but it potentially makes the image within a segment more appropriate for \code{ProFit} source fitting without doing full simultaneous fitting via \option{fit_extra}. Whilst not as a good a solution as simultaneous fitting, it is often better than nothing.
}
  \item{fit_extra}{
Logical; optional, should extra objects near to the primary target also be fitted in parallel? Default is FALSE. The other important option if this is TRUE is \option{magdiff}, that controls how bright (relatively speaking) the nearby sources need to be before they are modelled by \code{ProFit}. In multi-band mode you might want to set this to FALSE because it can be hard to predict how many extra objects are found in each band, and therefore constructing the correct \option{parm} input can be very complicated. Cannot be used in conjunction with \option{deblend_extra} (i.e. if \option{fit_extra} is TRUE then \option{deblend_extra} is treated as FALSE no matter what the user sets it to).
}
  \item{pos_delta}{
Numeric scalar; allowed range of [x,y] centring compared to the initial guess. The default of 10 works well for most cases, but very large galaxies with complex structure might need this set higher to better allow the true centre to be recovered.
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
\code{profitFound2Fit} produces everything you need to start fitting (which is usually just the \option{Data} structure), and also produces the intermediate \option{profound} output to check for the quality of the automatic segmentation etc. This means a user who wants more control over the final fitting can use the outputs of \code{profitFound2Fit} and put them into their own optimisation method. This route to fitting is probably for more experienced / confident users, and for people who want to optimise the \code{profoundProFound} stage of the process. \code{profitFound2Fit} also tightly crops the image down to the pixels needed for fitting (which improves ProFit fitting performance since the memory footprint is as small as it can be).

\code{profitDoFit} also carries out the optimisation, using the default settings of \code{Highlander}. For the typical user this might be all they need to get the results they want. Additional neighbouring sources that are detected during the \code{ProFound} stage are modelled as single Sersic profiles in a separate model list, which is usually reasonable enough to get a good solution for the target source of interest. To fit a single star with a Moffat profile you can use the \option{Ncomp} = 0.5 option, but in this case the model image will not be convolved with a target PSF image (because that is what we are trying to estimate). As such, if there are any other objects close by the Moffat profile might be compromised by the imperfect light profiling of these non-convolved single Sersic sources.

If you are doing more complex fitting (e.g. potentially using a different optimiser than \code{Highlander}, or wanting to make use of parameter priors) then you should probably use \code{profitFound2Fit} to create an intermediary version of the \option{Data} structure, which you can then edit and add to as required. The default combination of \option{optim_iters}, \option{Niters} and \option{NfinalMCMC} seems to work well in practice, and was used for the initial low redshift (z<=0.06) GAMA sample.

Re \option{fit_rough} and \option{star_rough}: setting to TRUE is less accurate, and potentially much less for compact sources, but it also runs much faster. Potentially if the target object is very extended then this will be accurate enough and speed up fitting by a factor of a few usually. It is also useful for quick testing. For publishable results this should usually be left as the default of FALSE. For \code{AllStar} (so \option{star_rough}) setting this to FALSE is perhaps a good option if you do no care about the model parameters themselves (e.g. FWHM) and just want a decent PSF image.
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
\code{\link{profuseMultiBandFound2Fit}}, \code{\link{profuseMultiBandDoFit}}, \code{\link{profuseAllStarFound2Fit}}, \code{\link{profuseAllStarDoFit}}, \code{\link{profitSetupData}}, \code{\link{profitLikeModel}}
}
\examples{
library(Rfits)
library(ProFit)

\dontrun{
#Highest level mode:

image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))$imDat

lazy = profuseDoFit(image=image, magzero=30, optim_iters=2, NfinalMCMC=200)

#We can have a look at the results in a similar manner to usual:

profitLikeModel(lazy$parm, Data=lazy$Data, makeplots=TRUE, plotchisq=TRUE)

#That was easy!

#There are complex types of fits that you might need to do more manually with ProFit,
#but the above appears to cover 99\% of fitting situations in practice.

#Doing things more manually you can also specify your own PSF (if you know it):

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

#Setting SBdilate to between 2-4 is often appropriate
example_f2f = profuseFound2Fit(image=image, psf=psf, magzero=30, SBdilate=2)

plot(example_f2f$profound) #The ProFound output
profitLikeModel(example_f2f$Data$init, Data=example_f2f$Data, makeplots=TRUE,
  plotchisq=TRUE) #Starting guess for optimisation. Target source is in the top-left.

example_fit = profuseDoFit(image=image, psf=psf, magzero=30, SBdilate=2, plot=TRUE)
}
}
\concept{ fit }

