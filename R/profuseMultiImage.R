profuseMultiImageFound2Fit = function(
    image_list,
    segim_list = NULL,
    mask_list = NULL,
    Ncomp = 2,
    loc = NULL,
    cutbox = NULL,
    psf_list = NULL,
    magzero = NULL,
    gain = NULL,
    resamp = NULL,
    tightcrop = TRUE,
    offset_list = NULL,
    ...
){
  Nim = length(image_list)

  if(is.null(magzero)){
    magzero = rep(0, Nim)
  }

  if(length(magzero) == 1){
    magzero = rep(magzero, Nim)
  }

  F2F_multi = list()

  for(i in 1:Nim){
    # if(!is.null(cutbox)){ # this means we are wanting to cutout around the target object and centre it
    #   if(is.null(loc)){
    #     loc = dim(image_list[[i]])/2
    #   }
    #   cutim = magicaxis::magcutout(image_list[[i]], loc=loc, box=cutbox)
    #   loc_cut = cutim$loc
    #   cutim = cutim$image
    #
    #   if(!is.null(segim_list[[i]])){#if we have segim, cut it out
    #     cutseg = magicaxis::magcutout(segim_list[[i]], loc = loc, box = cutbox)$image
    #   }else{
    #     cutseg = NULL
    #   }
    #
    #   if(!is.null(mask_list[[i]])){#if we have mask, cut it out
    #     cutmask = magicaxis::magcutout(mask_list[[i]], loc = loc, box = cutbox)$image
    #   }else{
    #     cutmask = is.na(cutim)
    #   }
    # }else{
    #   cutim = image_list[[i]]
    #   cutseg = segim_list[[i]]
    #   if(!is.null(mask_list[[i]])){
    #     cutmask = mask_list[[i]]
    #   }else{
    #     cutmask = is.na(cutim)
    #   }
    #   if(!is.null(loc)){
    #     loc_cut = loc
    #   }else{
    #     loc_cut = dim(cutim) / 2
    #   }
    # }

    if(is.null(loc)){
      loc = dim(image_list[[i]])/2
    }

    if(!is.null(cutbox) & !is.null(offset_list[[i]])){
      loc = c(loc[1] + offset_list[[i]][1], loc[2] + offset_list[[i]][2])
      offset_list[[i]] = loc - floor(loc)
      loc = floor(loc)
    }

    F2F = profuseFound2Fit(image = image_list[[i]],
                           segim = segim_list[[i]],
                           mask = mask_list[[i]],
                           loc = loc,
                           cutbox = cutbox,
                           Ncomp = Ncomp,
                           psf = psf_list[[i]],
                           magzero = magzero[[i]],
                           gain = gain[[i]],
                           resamp = resamp[[i]],
                           tightcrop = tightcrop,
                           offset = offset_list[[i]],
                           ...
    )$Data

    F2F_multi = c(F2F_multi, list(F2F))
  }

  F2F_multi$mon.names = F2F_multi[[1]]$mon.names
  F2F_multi$parm.names = F2F_multi[[1]]$parm.names
  F2F_multi$N = F2F_multi[[1]]$N
  F2F_multi$Nim = Nim

  names(F2F_multi)[1:Nim] = paste('image', 1:Nim, sep='')
  class(F2F_multi) = c(class(F2F_multi), 'F2F_multi')

  return(invisible(F2F_multi))
}

profuseMultiImageDoFit = function(image_list,
                        F2F = NULL,
                        Ncomp = 2,
                        psf_list = NULL,
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
    F2F = profuseMultiImageFound2Fit(
      image_list = image_list,
      Ncomp = Ncomp,
      psf_list = psf_list,
      magzero = magzero,
      ...
    )
  }else{
    F2F = profuseRegenPSF_MF2F(F2F)
  }

  Data = F2F[[1]]

  if (plot) {
    profitLikeModel(parm = Data$init,
                    Data = F2F,
                    makeplots = TRUE)
    legend('topright', legend = 'Start')
  }

  lowers = unlist(Data$intervals)[c(T, F)]
  lowers[unlist(Data$tolog) == T] = log10(lowers[unlist(Data$tolog) == T])
  lowers = as.numeric(lowers[which(unlist(Data$tofit))])
  uppers = unlist(Data$intervals)[c(F, T)]
  uppers[unlist(Data$tolog) == T] = log10(uppers[unlist(Data$tolog) == T])
  uppers = as.numeric(uppers[which(unlist(Data$tofit))])

  if(!is.null(Data$offset)){
    xcen_loc = grep('xcen',Data$parm.names)
    ycen_loc = grep('ycen',Data$parm.names)
    if(length(xcen_loc) > 0){
      lowers[xcen_loc] = lowers[xcen_loc] - Data$offset[1]
      uppers[xcen_loc] = uppers[xcen_loc] - Data$offset[1]
    }
    if(length(ycen_loc) > 0){
      lowers[ycen_loc] = lowers[ycen_loc] - Data$offset[2]
      uppers[ycen_loc] = uppers[ycen_loc] - Data$offset[2]
    }
  }

  message('Running Highlander')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The Highander package is required to run this function!')}
  highfit = Highlander::Highlander(
    parm = Data$init,
    Data = F2F,
    likefunc = profitLikeModel,
    seed = seed,
    lower = lowers,
    upper = uppers,
    applyintervals = FALSE,
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
    profitLikeModel(highfit$parm, Data=F2F, makeplots = TRUE)
    legend('topright', legend = 'After')
  }

  highfit$Data = F2F
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
