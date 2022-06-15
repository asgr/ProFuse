profuseMultiBandFound2Fit = function(image_list,
                                     segim_list = NULL,
                                     segim_global = NULL,
                                    sky_list = NULL,
                                    skyRMS_list = NULL,
                                    parm_global = c("sersic.xcen1", "sersic.ycen1", "sersic.re1", "sersic.ang2", "sersic.axrat2"),
                                    Ncomp = 2,
                                    loc = NULL,
                                    cutbox = NULL,
                                    psf_list = NULL,
                                    magdiff = 2.5,
                                    magzero = NULL,
                                    gain = NULL,
                                    resamp = NULL,
                                    doprofit = NULL,
                                    sing_nser = 2,
                                    bulge_nser = 4,
                                    disk_nser = 1,
                                    sing_nser_fit = TRUE,
                                    bulge_nser_fit = FALSE,
                                    disk_nser_fit = FALSE,
                                    bulge_circ =  TRUE,
                                    nser_upper=5.3,
                                    star_rough = TRUE,
                                    fit_rough = FALSE,
                                    psf_dim = c(51, 51),
                                    star_con = 2,
                                    star_con_fit = TRUE,
                                    star_circ = TRUE,
                                    tightcrop = TRUE,
                                    wave = NULL,
                                    smooth.parm = NULL,
                                    parm_ProSpect = NULL,
                                    data_ProSpect = NULL, #perhaps need a way to specify extra data going to bulge/disk. Naming or list?
                                    logged_ProSpect = NULL,
                                    intervals_ProSpect = NULL,
                                    autoclip = TRUE,
                                    roughpedestal = TRUE,
                                    ...){
  Nim = length(image_list)

  if(is.null(magzero)){
    magzero = rep(0, Nim)
  }

  if(length(magzero) == 1){
    magzero = rep(magzero, Nim)
  }

  for(i in 1:Nim){
    if(autoclip){
      image_med = median(image_list[[i]], na.rm=TRUE)
      image_list[[i]][image_list[[i]] - image_med < quantile(image_list[[i]] - image_med, 0.001, na.rm=TRUE)*100] = NA
      image_list[[i]][image_list[[i]] - image_med > quantile(image_list[[i]] - image_med, 0.999, na.rm=TRUE)*100] = NA
    }

    if(is.null(sky_list[i][[1]]) | is.null(skyRMS_list[i][[1]])){ #[i][[1]] looks silly, but it will return NULL when sky_list = NULL for any i (default). [[i]] will error in this case. This does not seem to be the case as of R v4.1.0, but probably leave to be safe for now.
      message("Image ",i,": running initial ProFound")
      profound = ProFound::profoundProFound(image = image_list[[i]],
                                            segim = segim_list[[i]],
                                            sky = sky_list[i][[1]],
                                            skyRMS = skyRMS_list[i][[1]],
                                            magzero = magzero[i],
                                            roughpedestal = roughpedestal,
                                            ...)
      if(is.null(sky_list[i])){
        image_list[[i]] = image_list[[i]] - profound$sky
      }else{
        image_list[[i]] = image_list[[i]] - sky_list[i]
      }
      if(is.null(skyRMS_list[i][[1]])){
        skyRMS_list[[i]] = profound$skyRMS
      }
    }

    if(is.null(psf_list[i][[1]]) & !isFALSE(doprofit[i])){ #last bit looks weird, but means the AllStar phase is skipped when doprofit is FALSE or NULL.
      message("Image ",i,": running AllStarDoFit")
      psf_list[[i]] = profuseAllStarDoFit(image = image_list[[i]],
                                          resamp = resamp[i][[1]],
                                         psf_dim = psf_dim,
                                         star_con = star_con,
                                         star_con_fit = star_con_fit,
                                         star_circ = star_circ,
                                         magzero = magzero[i],
                                         rough = star_rough,
                                         autoclip = FALSE,
                                         roughpedestal = roughpedestal,
                                         skycut = 2, #works well for stars
                                         SBdilate = 2)$psf #works well for stars
    }
    if(is.null(psf_list[i][[1]]) & isFALSE(doprofit[i])){
      psf_list[[i]] = matrix(1,1,1)
    }
  }

  message("Making image stack")

  if(is.null(doprofit)){
    multi_stack = ProFound::profoundMakeStack(
      image_list = image_list,
      skyRMS_list = skyRMS_list,
      magzero_in = magzero,
      magzero_out = 0
    )
  }else{
    multi_stack = ProFound::profoundMakeStack(
      image_list = image_list[doprofit],
      skyRMS_list = skyRMS_list[doprofit],
      magzero_in = magzero[doprofit],
      magzero_out = 0
    )
  }

  message("Running ProFound on stack")

  if(is.null(segim_global) & !is.null(segim_list)){
    if(is.null(doprofit)){
      segim_global = profuseSegimGlobal(segim_list)
    }else{
      segim_global = profuseSegimGlobal(segim_list[doprofit])
    }
  }

  multi_stack_pro = ProFound::profoundProFound(image=multi_stack$image,
                                               segim=segim_global,
                                               sky=0,
                                               skyRMS=multi_stack$skyRMS,
                                               redosky=FALSE,
                                               static_photom=TRUE,
                                               ...)

  message("Running Found2Fit on stack")

  F2Fstack = profuseFound2Fit(image = multi_stack$image,
                             sigma = multi_stack$skyRMS, #not quite a sigma map, but doesn't matter for the stack F2F
                             segim = multi_stack_pro$segim,
                             Ncomp = Ncomp,
                             loc = loc,
                             cutbox = cutbox,
                             psf = matrix(1,1,1), #Doesn't matter what we pass in here
                             magzero = 0,
                             mag_fit = is.null(parm_ProSpect),
                             sing_nser = sing_nser,
                             bulge_nser = bulge_nser,
                             disk_nser = disk_nser,
                             sing_nser_fit = sing_nser_fit,
                             bulge_nser_fit = bulge_nser_fit,
                             disk_nser_fit = disk_nser_fit,
                             bulge_circ =  bulge_circ,
                             nser_upper = nser_upper,
                             tightcrop = FALSE,
                             fit_extra = FALSE,
                             autoclip = FALSE,
                             roughpedestal = roughpedestal
  )

  if(tightcrop){
    regionlim = which(F2Fstack$Data$region, arr.ind=TRUE)

    xlo = min(regionlim[,1])
    xhi = max(regionlim[,1])
    ylo = min(regionlim[,2])
    yhi = max(regionlim[,2])

    if(!is.null(F2Fstack$Data$modellist$sersic)){
      F2Fstack$Data$modellist$sersic$xcen = F2Fstack$Data$modellist$sersic$xcen - xlo + 1
      F2Fstack$Data$modellist$sersic$ycen = F2Fstack$Data$modellist$sersic$ycen - ylo + 1

      for(i in 1:length(F2Fstack$Data$intervals$sersic$xcen)){
        F2Fstack$Data$intervals$sersic$xcen[[i]] = F2Fstack$Data$intervals$sersic$xcen[[i]] - xlo + 1
        F2Fstack$Data$intervals$sersic$ycen[[i]] = F2Fstack$Data$intervals$sersic$ycen[[i]] - ylo + 1
      }
    }

    if(!is.null(F2Fstack$Data$modellist$moffat)){
      F2Fstack$Data$modellist$moffat$xcen = F2Fstack$Data$modellist$moffat$xcen - xlo + 1
      F2Fstack$Data$modellist$moffat$ycen = F2Fstack$Data$modellist$moffat$ycen - ylo + 1

      for(i in 1:length(F2Fstack$Data$intervals$moffat$xcen)){
        F2Fstack$Data$intervals$moffat$xcen[[i]] = F2Fstack$Data$intervals$moffat$xcen[[i]] - xlo + 1
        F2Fstack$Data$intervals$moffat$ycen[[i]] = F2Fstack$Data$intervals$moffat$ycen[[i]] - ylo + 1
      }
    }

    if(!is.null(F2Fstack$Data$modellist$pointsource)){
      F2Fstack$Data$modellist$pointsource$xcen = F2Fstack$Data$modellist$pointsource$xcen - xlo + 1
      F2Fstack$Data$modellist$pointsource$ycen = F2Fstack$Data$modellist$pointsource$ycen - ylo + 1

      for(i in 1:length(F2Fstack$Data$intervals$pointsource$xcen)){
        F2Fstack$Data$intervals$pointsource$xcen[[i]] = F2Fstack$Data$intervals$pointsource$xcen[[i]] - xlo + 1
        F2Fstack$Data$intervals$pointsource$ycen[[i]] = F2Fstack$Data$intervals$pointsource$ycen[[i]] - ylo + 1
      }
    }

    xcenloc = grep('xcen', F2Fstack$Data$parm.names)
    if(length(xcenloc) > 0){
      F2Fstack$Data$init[xcenloc] = F2Fstack$Data$init[xcenloc] - xlo + 1
    }

    ycenloc = grep('ycen', F2Fstack$Data$parm.names)
    if(length(ycenloc) > 0){
      F2Fstack$Data$init[ycenloc] = F2Fstack$Data$init[ycenloc] - ylo + 1
    }

  }else{
    xlo = 1L
    ylo = 1L
    xhi = dim(multi_stack$image)[1]
    yhi = dim(multi_stack$image)[2]
  }

  if(!is.null(parm_ProSpect)){
    F2Fstack$Data$tofit$sersic$mag[] = FALSE
  }

  mag_stack = ProFound::profoundFlux2Mag(flux=sum(multi_stack$image[F2Fstack$Data$region], na.rm=TRUE), magzero=0)

  MF2F = list()

  for(i in 1:Nim){
    if(is.null(gain[[i]])){
     gain_loc = ProFound::profoundGainEst(image_list[[i]], objects=multi_stack_pro$objects, sky = 0)
    }else{
      gain_loc = gain[[i]]
    }
    sigma = ProFound::profoundMakeSigma(
      image = image_list[[i]],
      objects = multi_stack_pro$objects,
      gain = gain_loc,
      sky = 0,
      skyRMS = skyRMS_list[[i]],
      plot = FALSE
    )

    if(!is.null(resamp[i][[1]])){
      if(resamp[i][[1]] > 1){
        sigma = sigma*resamp[i][[1]]
      }
    }

    if(is.null(segim_list[i][[1]])){
      segim_use = F2Fstack$Data$segim
    }else{
      segim_use = segim_list[i][[1]]
    }

    message("Image ",i,": running SetupData")

    region = (segim_use == which.max(tabulate(segim_use[F2Fstack$Data$region])))

    MF2F[[i]] = profitSetupData(
      image = image_list[[i]][xlo:xhi,ylo:yhi],
      region = region[xlo:xhi,ylo:yhi],
      sigma = sigma[xlo:xhi,ylo:yhi],
      segim = segim_use[xlo:xhi,ylo:yhi],
      psf = psf_list[[i]],
      modellist = F2Fstack$Data$modellist,
      tofit = F2Fstack$Data$tofit,
      tolog = F2Fstack$Data$tolog,
      intervals = F2Fstack$Data$intervals,
      constraints = F2Fstack$Data$constraints,
      magzero = magzero[i],
      algo.func = 'LD',
      verbose = FALSE,
      rough = fit_rough
    )

    #if(!is.null(doprofit)){
    MF2F[[i]]$doprofit = doprofit[i]
    MF2F[[i]]$object_flux = magAB2Jansky(profoundFlux2Mag(sum(MF2F[[i]]$image*MF2F[[i]]$region, na.rm=TRUE), magzero[i])) #only used for fitting when doprofit = FALSE, but useful outputs anyway
    MF2F[[i]]$object_fluxerr = magAB2Jansky(profoundFlux2Mag(sqrt(sum(MF2F[[i]]$sigma^2*MF2F[[i]]$region, na.rm=TRUE)), magzero[i])) #only used for fitting when doprofit = FALSE, but useful outputs anyway
    #}
  }

  names(MF2F) = names(image_list)

  #Check below for ProFuse- how this all works could be complicated... see also profitLikeModel and profitRemakeModelList
  if(is.null(parm_global)){
    parm = F2Fstack$Data$init
    for(i in 1:Nim){
      MF2F[[i]]$parmuse = 1:length(parm)
    }
  }else{
    parm_init = F2Fstack$Data$init
    if(is.character(parm_global)){
      parm_global = match(parm_global, names(parm_init))
    }
    parm = parm_init[parm_global]
    Nparm = length(parm_init)
    parm_local = 1:Nparm
    parm_local = parm_local[-parm_global] #anything not global is local - check for ProFuse
    for(i in 1:Nim){
      if(length(parm_local) > 0){
        parm_temp = F2Fstack$Data$init[parm_local] #extract local - check for ProFuse
        names(parm_temp) = paste0(names(parm_temp),'_',i) #mod names for local by adding the band number - check for ProFuse
        parm = c(parm, parm_temp) #create parent parm object (with all unique parm_local vector appended) - check for ProFuse
        parmuse = 1:Nparm
        parmuse[parm_global] = 1:length(parm_global)
        parmuse[parm_local] = length(parm_global) + 1:length(parm_local) + (i-1)*length(parm_local) #define parmuse location in parent parm object - check for ProFuse
        MF2F[[i]]$parmuse = parmuse
      }else{
        MF2F[[i]]$parmuse = 1:Nparm
      }
    }
  }

#Ideas for ProFuse. Need to calculate the parmuse positions as per they will be after all the ProSpect related parameters are removed and the relevant per band magnitudes added and named. This probably means we need a parm_ProSpect object that exhaustively identifies all of these (and need to enforce them being at the end). Probably cannot make it any more flexible just to be safe. In principle then we just need remove these ProSpect related arguments and replace that part of the parm with the mag inside profitLikeModel, which we then pass into profitRemakeModelList to make the target model images.

  #Create mag offsets based on magzero points and average mag of the stack.
  for(i in 1:Nim){
    mag_image = ProFound::profoundFlux2Mag(flux=sum(image_list[[i]][F2Fstack$Data$region], na.rm=TRUE), magzero=magzero[i])
    mag_diff = mag_stack - mag_image
    sel = grep(paste0('.*mag.*\\_',i), names(parm))
    parm[sel] = parm[sel] - mag_diff
  }

  if(is.null(data_ProSpect$LumDist_Mpc)){
    if(!is.null(data_ProSpect$z)){
      data_ProSpect$LumDist_Mpc = cosdistLumDist(data_ProSpect$z, H0 = 67.8, OmegaM = 0.308)
    }
  }

  if(is.null(data_ProSpect$magemax)){
    if(!is.null(data_ProSpect$agemax)){
      data_ProSpect$magemax = data_ProSpect$agemax/1e9
    }
  }

  if(is.null(data_ProSpect$Zagemax)){
    if(!is.null(data_ProSpect$agemax)){
      data_ProSpect$Zagemax = data_ProSpect$agemax/1e9
    }
  }

  MF2F$init = c(parm, unlist(parm_ProSpect))
  MF2F$parm.names = names(MF2F$init)
  MF2F$mon.names = F2Fstack$Data$mon.names
  MF2F$Nim = Nim #Number of images
  MF2F$Ncomp = ceiling(Ncomp) #Number of components 0.5 -> 1 and 1.5 -> 2
  MF2F$N = F2Fstack$Data$N #This is the number of fitting pixels (cannot rename)
  MF2F$wave = wave
  MF2F$smooth.parm = smooth.parm
  MF2F$parm_ProSpect = parm_ProSpect
  MF2F$data_ProSpect = data_ProSpect
  MF2F$logged_ProSpect = logged_ProSpect
  MF2F$intervals_ProSpect = intervals_ProSpect

  class(MF2F) = c(class(MF2F), 'MF2F')

  return(invisible(MF2F))
}

profuseMultiBandDoFit = function(image_list,
                                MF2F = NULL,
                                parm_global = c("sersic.xcen1", "sersic.ycen1", "sersic.re1", "sersic.ang2", "sersic.axrat2"),
                                Ncomp = 2,
                                magzero = NULL,
                                seed = 666,
                                optim_iters = 5,
                                Niters = c(200,200),
                                NfinalMCMC = 1000,
                                walltime = Inf,
                                keepall = FALSE,
                                ...) {

  timestart = proc.time()[3] # start timer
  #call = match.call(expand.dots=TRUE)

  if(is.null(MF2F)){
    message('Running MultiBandFound2Fit')
    MF2F = profuseMultiBandFound2Fit(
      image_list = image_list,
      parm_global = parm_global,
      Ncomp = Ncomp,
      magzero = magzero,
      ...
    )
  }else{
    MF2F = profuseRegenPSF_MF2F(MF2F)
  }

  if(is.null(MF2F$smooth.parm)){
    #This implies we are in ProSpect fitting mode (this is how we use it now)
    lower_profit = {}
    upper_profit = {}
    logged_profit = {}

    for(i in 1:length(MF2F[[1]]$intervals)){ #loop over profiles
      for(j in 1:length(MF2F[[1]]$intervals[[i]])){ #loop over parameters
        for(k in 1:length(MF2F[[1]]$intervals[[i]][[1]])){ #loop over components
          if(isTRUE(MF2F[[1]]$tofit[[i]][[j]][[k]])){
            lower_profit = c(lower_profit, MF2F[[1]]$intervals[[i]][[j]][[k]][1])
            upper_profit = c(upper_profit, MF2F[[1]]$intervals[[i]][[j]][[k]][2])
            logged_profit = c(logged_profit, MF2F[[1]]$tolog[[i]][[j]][[k]])
          }
        }
      }
    }

    lower_profit[logged_profit] = log10(lower_profit[logged_profit])
    upper_profit[logged_profit] = log10(upper_profit[logged_profit])

    lower = c(lower_profit, MF2F$intervals_ProSpect$lo)
    upper = c(upper_profit, MF2F$intervals_ProSpect$hi)

  }else{
    #This implies we are in smooth.spline fitting mode, i.e. not using ProSpect (we don't use this)
    #Currently I have just set these to NULL to keep things working. These should really inherit sensible limits, but we aren't really using this functionality anyway.
    lower = NULL
    upper = NULL
  }

  message('Running Highander on multi-band data')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The Highander package is required to run this function!')}
  highfit = Highlander::Highlander(
    parm = MF2F$init,
    Data = MF2F,
    likefunc = profitLikeModel,
    seed = seed,
    lower = lower,
    upper = upper,
    applyintervals = TRUE,
    applyconstraints = FALSE,
    optim_iters = optim_iters,
    Niters = Niters,
    NfinalMCMC = NfinalMCMC,
    walltime = walltime,
    parm.names = MF2F$parm.names,
    keepall = FALSE
  )

  highfit$MF2F = MF2F
  highfit$error = apply(highfit$LD_last$Posterior1,
                        MARGIN = 2,
                        FUN = 'sd')

  if(!is.null(MF2F$smooth.parm) & !is.null(MF2F$wave)){
    namevec = names(MF2F$smooth.parm)
    highfit$parm_smooth = highfit$parm
    for(i in 1:length(MF2F$smooth.parm)){
      highfit$parm_smooth = .smooth_parm(parm=highfit$parm_smooth, MF2F$parm.names, extract=namevec[i], wave=MF2F$wave, func=MF2F$smooth.parm[[i]])
    }
  }else{
    highfit$parm_smooth = NULL
  }

  highfit$time = (proc.time()[3]-timestart)/60
  highfit$date = date()
  highfit$ProFit.version = packageVersion('ProFit')
  highfit$ProFound.version = packageVersion('ProFound')
  highfit$Highlander.version = packageVersion('Highlander')
  highfit$R.version = R.version
  highfit$LD_last$Call = NULL
  highfit$LD_last$Model = NULL
  highfit$RedChi2 = highfit$RedChi2/MF2F$Nim

  class(highfit) = 'profusemulti'

  return(highfit)
}

.smooth_parm = function(parm, parm.names, extract='mag1', wave, func=smooth.spline){
  parm_loc = grep(extract,parm.names)

  parm[parm_loc] = func(log(wave),parm[parm_loc])$y
  return(parm)
}

profuseSegimGlobal = function(segim_list){
  if(requireNamespace("imager", quietly = TRUE)){
    i=NULL
    Csegim_list = foreach(i=1:length(segim_list))%do%{
      temp=segim_list[[i]]
      temp[temp==0L]=NA
      imager::as.cimg(temp)
    }
    Csegim_list = imager::as.imlist(Csegim_list)
    segim_global = imager::parmed(Csegim_list, na.rm=TRUE)
    segim_global = as.matrix(segim_global)
    segim_global[is.na(segim_global)] = 0L
    return(segim_global)
  }else{
    stop('The imager package is needed for segim merging to work. Please install from CRAN.', call. = FALSE)
  }
}
