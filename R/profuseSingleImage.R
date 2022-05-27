profuseFound2Fit = function(image,
                           sigma = NULL,
                           segim = NULL,
                           mask = NULL,
                           Ncomp = 2,
                           loc = NULL,
                           cutbox = NULL,
                           psf = NULL,
                           magdiff = 2.5,
                           magzero = 0,
                           gain = NULL,
                           resamp = NULL,
                           loc_use = FALSE,
                           loc_fit = TRUE,
                           mag_fit = TRUE,
                           sing_nser = 2,
                           bulge_nser = 4,
                           disk_nser = 1,
                           sing_nser_fit = TRUE,
                           bulge_nser_fit = FALSE,
                           disk_nser_fit = FALSE,
                           bulge_circ = TRUE,
                           nser_upper=5.3,
                           star_rough = TRUE,
                           fit_rough = FALSE,
                           psf_dim = c(51, 51),
                           star_con = 2,
                           star_con_fit = TRUE,
                           star_circ = TRUE,
                           offset = NULL,
                           tightcrop = TRUE,
                           deblend_extra = FALSE,
                           fit_extra = FALSE,
                           pos_delta = 10,
                           autoclip = TRUE,
                           roughpedestal = TRUE,
                           ...) {
  if(autoclip){
    image_med = median(image, na.rm=TRUE)
    image[image - image_med < quantile(image - image_med, 0.001, na.rm=TRUE)*100] = NA
    image[image - image_med > quantile(image - image_med, 0.999, na.rm=TRUE)*100] = NA
  }

  if(Ncomp == 0.5){psf = NULL}
  if((Ncomp >= 1 & is.null(psf)) | (Ncomp == 0 & is.null(psf))){
    message('Need PSF for Ncomp >= 1 or Ncomp == 0. Running AllStarDoFit!')
    psf = profuseAllStarDoFit(image = image,
                              resamp = resamp,
                              psf_dim = psf_dim,
                              star_con = star_con,
                              star_con_fit = star_con_fit,
                              star_circ = star_circ,
                              magzero = magzero,
                              rough = star_rough,
                              autoclip = FALSE,
                              skycut = 2, #works well for stars
                              SBdilate = 2)$psf #works well for stars
  }

  if(!is.null(loc) & !is.null(cutbox)){ # this means we are wanting to cutout around the target object and centre it
    cutim = magicaxis::magcutout(image, loc=loc, box=cutbox)
    loc_cut = cutim$loc
    cutim = cutim$image

    if(!is.null(sigma)){#if we have sigma, cut it out
      cutsigma = magicaxis::magcutout(sigma, loc = loc, box = cutbox)$image
    }else{
      cutsigma = NULL
    }

    if(!is.null(segim)){#if we have segim, cut it out
      cutseg = magicaxis::magcutout(segim, loc = loc, box = cutbox)$image
    }else{
      cutseg = NULL
    }

    if(!is.null(mask)){#if we have mask, cut it out
      cutmask = magicaxis::magcutout(mask, loc = loc, box = cutbox)$image
    }else{
      cutmask = is.na(cutim)
    }
  }else{
    cutim = image
    cutsigma = sigma
    cutseg = segim
    if(!is.null(mask)){
      cutmask = mask
    }else{
      cutmask = is.na(cutim)
    }
    if(!is.null(loc)){
      loc_cut = loc
    }else{
      loc_cut = dim(cutim) / 2
    }
  }

  # if(is.null(segim)){
  #   cutseg = NULL
  # }else if(dim(segim)[1] == dim(cutim)[1] & dim(segim)[2] == dim(cutim)[2]){
  #   cutseg = segim
  # }else if (dim(segim)[1] == dim(image)[1] & dim(segim)[2] == dim(image)[2] & docutouts) {
  #   cutseg = magicaxis::magcutout(segim, loc = loc, box = cutbox)$image
  # }else{
  #   message('No input segim that matches the input image- will create one using ProFound!')
  #   cutseg = NULL
  # }

  # if(is.null(mask)){
  #   cutmask = is.na(cutim)
  # }else if(dim(mask)[1] == dim(cutim)[1] & dim(mask)[2] == dim(cutim)[2]){
  #   cutmask = mask
  # }else if (dim(mask)[1] == dim(image)[1] & dim(mask)[2] == dim(image)[2] & !is.null(loc)) {
  #   cutmask = magicaxis::magcutout(mask, loc = loc, box = cutbox)$image
  # }else{
  #   cutmask = is.na(cutim)
  # }

  message('    Running ProFound')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The ProFound package is required to run this function!')}

  if(is.null(cutseg)){
    mini_profound = ProFound::profoundProFound(
      image = cutim,
      mask = cutmask,
      sky = 0,
      redosky = FALSE,
      nearstats = TRUE,
      groupby = 'segim',
      magzero = magzero,
      verbose = FALSE,
      roughpedestal = roughpedestal,
      ...
    )
    cutseg = mini_profound$segim
  }else{
    mini_profound = ProFound::profoundProFound(
      image = cutim,
      segim = cutseg,
      mask = cutmask,
      nearstats = TRUE,
      groupby = 'segim',
      magzero = magzero,
      verbose = FALSE,
      static_photom = TRUE,
      ...
    )
  }

  if (is.null(cutsigma)) {
    if(is.null(gain)){
      gain = ProFound::profoundGainEst(cutim, objects = mini_profound$objects, sky = 0)
    }
    cutsigma = ProFound::profoundMakeSigma(
      image = cutim,
      objects = mini_profound$objects,
      gain = gain,
      sky = 0,
      skyRMS = mini_profound$skyRMS,
      plot = FALSE
    )
  }

  if(!is.null(resamp)){
    if(resamp > 1){
      cutsigma = cutsigma*resamp
    }
  }

  if(deblend_extra & fit_extra==FALSE){
    cutim = ProFound::profoundFluxDeblend(mini_profound, image_reweight=TRUE)$image
  }

  segID_tar = mini_profound$segim[ceiling(loc_cut[1]), ceiling(loc_cut[2])]
  if (segID_tar == 0) {
    message('Target appears to be sky! Consider using different ProFound parameters.')
    return(NULL)
  }
  loc_tar = which(mini_profound$segstats$segID == segID_tar)
  magID_tar = mini_profound$segstats[mini_profound$segstats$segID == segID_tar, 'mag']

  if(fit_extra){
    segID_ext = unlist(mini_profound$near$nearID[mini_profound$near$segID == segID_tar])
    if (length(segID_ext) > 0) {
      loc_ext = match(segID_ext, mini_profound$segstats$segID)
      loc_ext = loc_ext[which(mini_profound$segstats[loc_ext, "mag"] < magID_tar + magdiff)]
      segID_ext = mini_profound$segstats[loc_ext, 'segID']
      N_ext = length(loc_ext)
    } else{
      N_ext = 0
    }
  }else{
    segID_ext = {}
    N_ext = 0
  }

  region = matrix(mini_profound$segim %in% c(segID_tar, segID_ext), nrow=dim(mini_profound$segim)[1], ncol=dim(mini_profound$segim)[2])
  if(!is.null(mini_profound$mask)){
    region = region & mini_profound$mask==0L
  }
  regionlim = which(region, arr.ind=TRUE)

  if(tightcrop){
    xlo = min(regionlim[,1])
    xhi = max(regionlim[,1])
    ylo = min(regionlim[,2])
    yhi = max(regionlim[,2])

    cutim = cutim[xlo:xhi, ylo:yhi]
    region = region[xlo:xhi, ylo:yhi]
    cutsigma = cutsigma[xlo:xhi, ylo:yhi]
    cutseg = cutseg[xlo:xhi, ylo:yhi]
    cutmask = cutmask[xlo:xhi, ylo:yhi]

    if(is.null(offset)){
      if(loc_use){
        xcen = loc[1] - xlo + 1L
        ycen = loc[2] - ylo + 1L
      }else{
        xcen = mini_profound$segstats[loc_tar, 'xmax'] - xlo + 1L
        ycen = mini_profound$segstats[loc_tar, 'ymax'] - ylo + 1L
      }
      xcen_int = xcen + c(-pos_delta,pos_delta)
      ycen_int = ycen + c(-pos_delta,pos_delta)
    }else{
      offset[1] = offset[1] - xlo + 1L
      offset[2] = offset[2] - ylo + 1L
      if(loc_use){
        xcen = loc[1]
        ycen = loc[2]
      }else{
        xcen = mini_profound$segstats[loc_tar, 'xmax']
        ycen = mini_profound$segstats[loc_tar, 'ymax']
      }
      xcen_int = xcen + c(-pos_delta,pos_delta) + offset[1]
      ycen_int = ycen + c(-pos_delta,pos_delta) + offset[2]
    }
  }else{
    xlo = 1L
    ylo = 1L
    if(loc_use){
      xcen = loc[1]
      ycen = loc[2]
    }else{
      xcen = mini_profound$segstats[loc_tar, 'xmax']
      ycen = mini_profound$segstats[loc_tar, 'ymax']
    }
    xcen_int = xcen + c(-pos_delta,pos_delta)
    ycen_int = ycen + c(-pos_delta,pos_delta)
  }

  if (Ncomp == 0) {
    modellist = list(
      pointsource = list(
        xcen = mini_profound$segstats[loc_tar, 'xmax'],
        ycen = mini_profound$segstats[loc_tar, 'ymax'],
        mag = mini_profound$segstats[loc_tar, 'mag']
      )
    )
  }else if(Ncomp == 0.5) {
    if(star_circ){
      ang = 0
      axrat = 1
    }else{
      ang = mini_profound$segstats[loc_tar, 'ang']
      axrat = mini_profound$segstats[loc_tar, 'axrat']
    }

    modellist = list(
      moffat = list(
          xcen = xcen,
          ycen = ycen,
          mag = mini_profound$segstats[loc_tar, 'mag'],
          fwhm = mini_profound$segstats[loc_tar, 'R50'] * 2,
          con = star_con,
          ang = ang,
          axrat = axrat
      )
    )
  } else if (Ncomp == 1) {
    modellist = list(
      sersic = list(
        xcen = xcen,
        ycen = ycen,
        mag = mini_profound$segstats[loc_tar, 'mag'],
        re = mini_profound$segstats[loc_tar, 'R50'],
        nser = sing_nser,
        ang = mini_profound$segstats[loc_tar, 'ang'],
        axrat = mini_profound$segstats[loc_tar, 'axrat']
      )
    )
  } else if (Ncomp == 2) {
    modellist = list(
      sersic = list(
        xcen = rep(xcen, 2),
        ycen = rep(ycen, 2),
        mag = rep(mini_profound$segstats[loc_tar, 'mag'], 2) + 0.752575,
        re = mini_profound$segstats[loc_tar, 'R50'] * c(0.5, 1.5),
        nser = c(bulge_nser, disk_nser),
        ang = c(ifelse(bulge_circ, 0, mini_profound$segstats[loc_tar, 'ang']), mini_profound$segstats[loc_tar, 'ang']),
        axrat = c(ifelse(bulge_circ, 1, mini_profound$segstats[loc_tar, 'axrat']), mini_profound$segstats[loc_tar, 'axrat'])
      )
    )
  } else if (Ncomp == 1.5) {
    modellist = list(
      pointsource = list(
        xcen = xcen,
        ycen = ycen,
        mag = mini_profound$segstats[loc_tar, 'mag'] + 0.752575
      ),
      sersic = list(
        xcen = xcen,
        ycen = ycen,
        mag = mini_profound$segstats[loc_tar, 'mag'] + 0.752575,
        re = mini_profound$segstats[loc_tar, 'R50'],
        nser = disk_nser,
        ang = mini_profound$segstats[loc_tar, 'ang'],
        axrat = mini_profound$segstats[loc_tar, 'axrat']
      )
    )
  } else if (Ncomp == 3) {
    modellist = list(
      sersic = list(
        xcen = rep(xcen, 3),
        ycen = rep(ycen, 3),
        mag = rep(mini_profound$segstats[loc_tar, 'mag'], 3) + 0.4771213,
        re = mini_profound$segstats[loc_tar, 'R50'] * c(0.5, 2, 1),
        nser = c(bulge_nser, disk_nser, disk_nser),
        ang = c(ifelse(bulge_circ, 0, mini_profound$segstats[loc_tar, 'ang']), rep(mini_profound$segstats[loc_tar, 'ang'],2)),
        axrat = c(ifelse(bulge_circ, 1, mini_profound$segstats[loc_tar, 'axrat']), rep(mini_profound$segstats[loc_tar, 'axrat'],2))
      )
    )
  }

  if (Ncomp == 0) {
    tofit = list(
      pointsource = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = mag_fit
      )
    )
    constraints = NULL
  }else if (Ncomp == 0.5) {
    tofit = list(
      moffat = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = mag_fit,
        fwhm = TRUE,
        con = star_con_fit,
        ang = !star_circ,
        axrat = !star_circ
      )
    )
    constraints = NULL
  }else if (Ncomp == 1) {
    tofit = list(
      sersic = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = mag_fit,
        re = TRUE,
        nser = sing_nser_fit,
        ang = TRUE,
        axrat = TRUE
      )
    )
    constraints = NULL
  } else if (Ncomp == 2) {
    tofit = list(sersic = list(
      xcen = c(loc_fit, NA), #The NA couples the components together
      ycen = c(loc_fit, NA), #The NA couples the components together
      mag = rep(mag_fit, 2),
      re = rep(TRUE, 2),
      nser = c(bulge_nser_fit, disk_nser_fit),
      ang = c(!bulge_circ, TRUE),
      axrat = c(!bulge_circ, TRUE)
    ))
    constraints = NULL
  } else if (Ncomp == 1.5) {
    tofit = list(
      pointsource = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = mag_fit
      ),
      sersic = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = mag_fit,
        re = TRUE,
        nser = disk_nser_fit,
        ang = TRUE,
        axrat = TRUE
      )
    )
    constraints = .constraints_func1_5
  } else if (Ncomp == 3) {
    tofit = list(sersic = list(
      xcen = c(loc_fit, NA, NA), #The NA couples the components together
      ycen = c(loc_fit, NA, NA), #The NA couples the components together
      mag = rep(mag_fit, 3),
      re = rep(TRUE, 3),
      nser = c(bulge_nser_fit, disk_nser_fit, NA),
      ang = c(!bulge_circ, TRUE, NA),
      axrat = c(!bulge_circ, TRUE, NA)
    ))
    constraints = NULL
  }

  if (Ncomp == 0) {
    tolog = list(
      pointsource = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = FALSE
      )
    )
    constraints = NULL
  }else if (Ncomp == 0.5) {
    tolog = list(
      moffat = list(
        xcen = rep(FALSE, Ncomp),
        ycen = rep(FALSE, Ncomp),
        mag = rep(FALSE, Ncomp),
        fwhm = rep(TRUE, Ncomp),
        #fwhm is best fit in log space
        con = rep(TRUE, Ncomp),
        #con is best fit in log space
        ang = rep(FALSE, Ncomp),
        axrat = rep(TRUE, Ncomp) #axrat is best fit in log space
      )
    )
  } else if (Ncomp == 1 | Ncomp == 2) {
    tolog = list(
      sersic = list(
        xcen = rep(FALSE, Ncomp),
        ycen = rep(FALSE, Ncomp),
        mag = rep(FALSE, Ncomp),
        re = rep(TRUE, Ncomp),
        #re is best fit in log space
        nser = rep(TRUE, Ncomp),
        #nser is best fit in log space
        ang = rep(FALSE, Ncomp),
        axrat = rep(TRUE, Ncomp) #axrat is best fit in log space
      )
    )
  } else if (Ncomp == 1.5) {
    tolog = list(
      pointsource = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = FALSE
      ),
      sersic = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = FALSE,
        re = TRUE,
        #re is best fit in log space
        nser = TRUE,
        #nser is best fit in log space
        ang = FALSE,
        axrat = TRUE #axrat is best fit in log space
      )
    )
  } else if (Ncomp == 3) {
    tolog = list(
      sersic = list(
        xcen = rep(FALSE, 3),
        ycen = rep(FALSE, 3),
        mag = rep(FALSE, 3),
        re = rep(TRUE, 3),
        #re is best fit in log space
        nser = rep(TRUE, 3),
        #nser is best fit in log space
        ang = rep(FALSE, 3),
        axrat = rep(TRUE, 3) #axrat is best fit in log space
      )
    )
  }

  #maxsize = sqrt(dim(cutim)[1]^2 + dim(cutim)[2]^2)
  maxsize = mini_profound$segstats[loc_tar, 'R50'] * 4

  if (Ncomp == 0) {
    intervals = list(pointsource = list(
      xcen = list(xcen_int),
      ycen = list(ycen_int),
      mag = list(c(0, 40))
    ))
  } else if (Ncomp == 0.5) {
    intervals = list(moffat = list(
      xcen = list(xcen_int),
      ycen = list(ycen_int),
      mag = list(c(0, 40)),
      fwhm = list(c(0.5, maxsize)),
      con = list(c(1, 10)),
      ang = list(c(-180, 360)),
      axrat = list(c(0.5, 1))
    ))
  } else if (Ncomp == 1) {
    intervals = list(sersic = list(
      xcen = list(xcen_int),
      ycen = list(ycen_int),
      mag = list(c(0, 40)),
      re = list(c(1, maxsize)),
      nser = list(c(0.5, nser_upper)),
      ang = list(c(-180, 360)),
      axrat = list(c(0.01, 1))
    ))
  } else if (Ncomp == 2) {
    intervals = list(
      sersic = list(
        xcen = list(xcen_int, xcen_int),
        ycen = list(ycen_int, ycen_int),
        mag = list(c(0, 40), c(0, 40)),
        re = list(c(1, maxsize), c(1, maxsize)),
        nser = list(c(2, nser_upper), c(0.5, 2)),
        ang = list(c(-180, 360), c(-180, 360)),
        axrat = list(c(0.01, 1), c(0.01, 1))
      )
    )
  } else if (Ncomp == 1.5) {
    intervals = list(
      pointsource = list(
        xcen = list(xcen_int),
        ycen = list(ycen_int),
        mag = list(c(0, 40))
      ),
      sersic = list(
        xcen = list(xcen_int),
        ycen = list(ycen_int),
        mag = list(c(0, 40)),
        re = list(c(1, maxsize)),
        nser = list(c(0.5, nser_upper)),
        ang = list(c(-180, 360)),
        axrat = list(c(0.01, 1))
      )
    )
  } else if (Ncomp == 3) {
    intervals = list(
      sersic = list(
        xcen = list(xcen_int, xcen_int, xcen_int),
        ycen = list(ycen_int, ycen_int, ycen_int),
        mag = list(c(0, 40), c(0, 40), c(0, 40)),
        re = list(c(1, maxsize), c(1, maxsize), c(1, maxsize)),
        nser = list(c(2, nser_upper), c(0.5, 2), c(0.5, 2)),
        ang = list(c(-180, 360), c(-180, 360), c(-180, 360)),
        axrat = list(c(0.01, 1), c(0.01, 1), c(0.01, 1))
      )
    )
  }

  if (fit_extra & N_ext > 0) {
    modellist = c(modellist,
                  list(
                    sersic = list(
                      xcen = mini_profound$segstats[loc_ext, 'xmax'] - xlo + 1L,
                      ycen = mini_profound$segstats[loc_ext, 'ymax'] - ylo + 1L,
                      mag = mini_profound$segstats[loc_ext, 'mag'],
                      re = mini_profound$segstats[loc_ext, 'R50'],
                      nser = rep(2, N_ext),
                      ang = mini_profound$segstats[loc_ext, 'ang'],
                      axrat = mini_profound$segstats[loc_ext, 'axrat']
                    )
                  )
                )

    tofit = c(tofit,
              list(
                sersic = list(
                  xcen = rep(FALSE, N_ext),
                  ycen = rep(FALSE, N_ext),
                  mag = rep(TRUE, N_ext),
                  re = rep(TRUE, N_ext),
                  nser = rep(TRUE, N_ext),
                  ang = rep(FALSE, N_ext),
                  axrat = rep(TRUE, N_ext)
                )
              )
            )

    tolog = c(tolog,
              list(
                sersic = list(
                  xcen = rep(FALSE, N_ext),
                  ycen = rep(FALSE, N_ext),
                  mag = rep(FALSE, N_ext),
                  re = rep(TRUE, N_ext),
                  #re is best fit in log space
                  nser = rep(TRUE, N_ext),
                  #nser is best fit in log space
                  ang = rep(FALSE, N_ext),
                  axrat = rep(TRUE, N_ext) #axrat is best fit in log space
                )
              )
            )

    maxsize = max(mini_profound$segstats[loc_ext, 'R50']*4, na.rm=TRUE)

    intervals = c(intervals,
                  list(
                    sersic = list(
                      xcen = rep(list(xcen_int), N_ext),
                      ycen = rep(list(ycen_int), N_ext),
                      mag = rep(list(c(0, 40)), N_ext),
                      re = rep(list(c(1, maxsize)), N_ext),
                      nser = rep(list(c(0.5, nser_upper)), N_ext),
                      ang = rep(list(c(-180, 360)), N_ext),
                      axrat = rep(list(c(0.01, 1)), N_ext)
                    )
                  )
                )
  }

  Data = profitSetupData(
    image = cutim,
    region = region,
    sigma = cutsigma,
    segim = cutseg,
    mask = cutmask,
    psf = psf,
    modellist = modellist,
    tofit = tofit,
    tolog = tolog,
    intervals = intervals,
    constraints = constraints,
    magzero = magzero,
    algo.func = 'LD',
    verbose = FALSE,
    offset = offset,
    rough = fit_rough
  )
  Data$Nmod = Ncomp + N_ext

  F2F_single = list(profound = mini_profound, Data = Data)
  class(F2F_single) = c(class(F2F_single), 'F2F_single')

  return(invisible(F2F_single))
}

profuseDoFit = function(image,
                       F2F = NULL,
                       Ncomp = 2,
                       psf = NULL,
                       magzero = NULL,
                       psf_dim = c(51,51),
                       plot = FALSE,
                       seed = 666,
                       optim_iters = 5,
                       Niters = c(200,200),
                       NfinalMCMC = 1000,
                       walltime = Inf,
                       keepall = FALSE,
                       ...) {

  timestart = proc.time()[3] # start timer
  #call = match.call(expand.dots=TRUE)

  if(is.null(F2F)){
    message('Running Found2Fit')
    F2F = profuseFound2Fit(
      image = image,
      Ncomp = Ncomp,
      psf = psf,
      magzero = magzero,
      ...
    )
  }else{
    F2F = profuseRegenPSF_F2F(F2F)
  }

  Data = F2F$Data

  if (plot) {
    profitLikeModel(parm = Data$init,
                    Data = Data,
                    makeplots = TRUE)
    legend('topright', legend = 'Start')
  }

  lowers = unlist(Data$intervals)[c(T, F)]
  lowers[unlist(Data$tolog) == T] = log10(lowers[unlist(Data$tolog) == T])
  lowers = lowers[which(unlist(Data$tofit))]
  uppers = unlist(Data$intervals)[c(F, T)]
  uppers[unlist(Data$tolog) == T] = log10(uppers[unlist(Data$tolog) == T])
  uppers = uppers[which(unlist(Data$tofit))]

  message('Running Highlander')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The Highander package is required to run this function!')}
  highfit = Highlander::Highlander(
    parm = Data$init,
    Data = Data,
    likefunc = profitLikeModel,
    seed = seed,
    lower = lowers,
    upper = uppers,
    applyintervals = TRUE,
    applyconstraints = FALSE,
    optim_iters = optim_iters,
    Niters = Niters,
    NfinalMCMC = NfinalMCMC,
    walltime = walltime,
    parm.names = Data$parm.names,
    keepall = keepall
  )
  names(highfit$parm) = names(Data$init)

  if (plot) {
    profitLikeModel(highfit$parm, Data = Data, makeplots = TRUE)
    legend('topright', legend = 'After')
  }

  highfit$profound = F2F$profound
  highfit$Data = Data
  highfit$initmodel = profitRemakeModellist(Data$init, Data = Data)
  highfit$finalmodel = profitRemakeModellist(highfit$parm, Data = Data)

  highfit$parm[highfit$parm < lowers] = lowers[highfit$parm < lowers]
  highfit$parm[highfit$parm > uppers] = uppers[highfit$parm > uppers]

  highfit$error = apply(highfit$LD_last$Posterior1,
                         MARGIN = 2,
                         FUN = 'sd')

  if(Ncomp == 0.5){
    if(psf_dim[1] %% 2 == 0){psf_dim[1] = psf_dim[1] + 1}
    if(psf_dim[2] %% 2 == 0){psf_dim[2] = psf_dim[2] + 1}
    temp_modellist = highfit$finalmodel$modellist[[1]]
    temp_modellist$moffat$mag = 0
    temp_modellist$moffat$xcen = psf_dim[1]/2
    temp_modellist$moffat$ycen = psf_dim[2]/2

    highfit$psf = profitMakeModel(temp_modellist, dim=psf_dim)
    highfit$psf_modellist = temp_modellist

    psf_fluxcheck = sum(highfit$psf$z)

    if(psf_fluxcheck < 0.95){
      message('WARNING: psf output image contains less than 95% of the total model flux! Consider increasing the size of psf_dim.')
    }

    if(psf_fluxcheck > 0.999){
      message('WARNING: psf output image contains more than 99.9% of the total model flux! Consider decreasing the size of psf_dim.')
    }
    highfit$psf = highfit$psf$z / sum(highfit$psf$z)
    highfit$psf_fluxcheck
  }

  highfit$time = (proc.time()[3]-timestart)/60
  highfit$date = date()
  #highfit$call = call
  highfit$ProFit.version = packageVersion('ProFit')
  highfit$ProFound.version = packageVersion('ProFound')
  highfit$Highlander.version = packageVersion('Highlander')
  highfit$R.version = R.version
  highfit$LD_last$Call = NULL
  highfit$LD_last$Model = NULL

  return(highfit)
}

.constraints_func1_5 = function(modellist=NULL) {
  modellist$pointsource$xcen = modellist$sersic$xcen
  modellist$pointsource$ycen = modellist$sersic$ycen
  return(modellist)
}

# .constraints_func3 = function(modellist=NULL) { #I don't think I actually need this, can use NA
#   modellist$sersic$nser[3] = modellist$sersic$nser[2]
#   modellist$sersic$ang[3] = modellist$sersic$ang[2]
#   modellist$sersic$axrat[3] = modellist$sersic$axrat[2]
#   return(modellist)
# }
