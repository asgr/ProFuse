profuseRegenPSF_F2F = function(F2F){
  F2F$convopt = list(convolver = profitMakeConvolver('brute', dim(F2F$image), F2F$psf), openclenv=NULL)
  return(F2F)
}

profuseRegenPSF_MF2F = function(MF2F){
  for(i in 1:MF2F$Nim){
    MF2F[[i]]$convopt = list(convolver = profitMakeConvolver('brute', dim(MF2F[[i]]$image), MF2F[[i]]$psf), openclenv=NULL)
  }
  return(MF2F)
}
