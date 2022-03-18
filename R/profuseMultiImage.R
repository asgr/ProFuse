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

  MF2F = list()

  for(i in 1:Nim){
    if(!is.null(loc) & !is.null(cutbox)){ # this means we are wanting to cutout around the target object and centre it
      cutim = magicaxis::magcutout(image_list[[i]], loc=loc, box=cutbox)
      loc_cut = cutim$loc
      cutim = cutim$image

      if(!is.null(segim_list[[i]])){#if we have segim, cut it out
        cutseg = magicaxis::magcutout(segim_list[[i]], loc = loc, box = cutbox)$image
      }

      if(!is.null(mask_list[[i]])){#if we have mask, cut it out
        cutmask = magicaxis::magcutout(mask_list[[i]], loc = loc, box = cutbox)$image
      }else{
        cutmask = is.na(cutim)
      }
    }else{
      cutim = image_list[[i]]
      cutseg = segim_list[[i]]
      if(!is.null(mask_list[[i]])){
        cutmask = mask_list[[i]]
      }else{
        cutmask = is.na(cutim)
      }
      if(!is.null(loc)){
        loc_cut = loc
      }else{
        loc_cut = dim(cutim) / 2
      }
    }

    if(!is.null(offset_list[[i]])){
      loc_cut = c(loc_cut[1] + offset_list[[i]][1], loc_cut[2] + offset_list[[i]][2])
    }

    F2F = profuseFound2Fit(image = cutim,
                           segim = cutseg,
                           mask = cutmask,
                           loc = loc_cut,
                           Ncomp = Ncomp,
                           psf = psf_list[[i]],
                           magzero = magzero[[i]],
                           gain = gain[[i]],
                           resamp = resamp[[i]],
                           tightcrop = FALSE,
                           offset = offset_list[[i]],
                           ...
    )$Data

    MF2F = c(MF2F, list(F2F))
  }

  names(MF2F) = paste('image', 1:Nim, sep='')

  return(MF2F)
}
