test_that(".smooth_parm only modifies the matching parameters", {
  parm = c(mag1_a = 20, mag1_b = 21, mag1_c = 19, other = 5)
  wave = c(4000, 6000, 8000)
  out = ProFuse:::.smooth_parm(parm, names(parm), extract = 'mag1', wave = wave,
                               func = function(x, y) list(y = rep(mean(y), length(y))))
  expect_equal(unname(out[1:3]), rep(20, 3))
  expect_equal(out[['other']], 5)
  expect_equal(names(out), names(parm))
})

test_that(".constraints_func0_5 ties pointsource centre to first sersic", {
  ml = list(sersic = list(xcen = c(10, 20), ycen = c(11, 21)),
            pointsource = list(xcen = 0, ycen = 0))
  out = ProFuse:::.constraints_func0_5(ml)
  expect_equal(out$pointsource$xcen, 10)
  expect_equal(out$pointsource$ycen, 11)
  expect_equal(out$sersic, ml$sersic)
})

test_that("profuseSegimGlobal merges segmentation maps", {
  skip_if_not_installed("imager")
  s1 = matrix(0L, 10, 10); s1[3:5, 3:5] = 1L
  s2 = matrix(0L, 10, 10); s2[3:5, 3:5] = 1L
  out = profuseSegimGlobal(list(s1, s2))
  expect_equal(dim(out), c(10, 10))
  expect_false(anyNA(out))
  expect_true(all(out[3:5, 3:5] == 1L))
  expect_true(all(out[1:2, ] == 0L))
})

test_that("profuseRegenPSF_F2F rebuilds a convolver", {
  psf = matrix(0, 5, 5); psf[3, 3] = 1
  F2F = list(image = matrix(0, 20, 20), psf = psf)
  out = profuseRegenPSF_F2F(F2F)
  expect_named(out$convopt, c('convolver', 'openclenv'))
  expect_null(out$convopt$openclenv)
})

test_that("profuseRegenPSF_MF2F regenerates every image", {
  psf = matrix(0, 5, 5); psf[3, 3] = 1
  im = list(image = matrix(0, 20, 20), psf = psf)
  MF2F = list(im, im); MF2F$Nim = 2L
  out = profuseRegenPSF_MF2F(MF2F)
  expect_false(is.null(out[[1]]$convopt$convolver))
  expect_false(is.null(out[[2]]$convopt$convolver))
})
