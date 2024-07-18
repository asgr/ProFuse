profuseRegenPSF_F2F = function(F2F, nbench=0L){
  if(nbench > 0L){
    best = profitBenchmarkConv(F2F$image, psf=F2F$psf, nbench=nbench)$best
    F2F$convopt = list(convolver = profitMakeConvolver(best, dim(F2F$image), F2F$psf), openclenv=NULL)
  }else{
    F2F$convopt = list(convolver = profitMakeConvolver('brute', dim(F2F$image), F2F$psf), openclenv=NULL)
  }
  return(F2F)
}

profuseRegenPSF_MF2F = function(MF2F, nbench=0L){
  for(i in 1:MF2F$Nim){
    if(nbench > 0L){
      best = profitBenchmarkConv(MF2F[[i]]$image, psf=MF2F[[i]]$psf, nbench=nbench)$best
      MF2F[[i]]$convopt = list(convolver = profitMakeConvolver(best, dim(MF2F[[i]]$image), MF2F[[i]]$psf), openclenv=NULL)
    }else{
      MF2F[[i]]$convopt = list(convolver = profitMakeConvolver('brute', dim(MF2F[[i]]$image), MF2F[[i]]$psf), openclenv=NULL)
    }
  }
  return(MF2F)
}
