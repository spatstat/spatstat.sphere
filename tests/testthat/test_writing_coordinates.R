library(spatstat.sphere)
context("Writing/exporting coordinates.")

test_that("Coordinates in x,y,z format are handled", {
  lat <- c(-45,45)
  lon <- c(0,90)
  xyz <- globe::ensure3d(data.frame(lon, lat), single = FALSE)
  colnames(xyz) <- c("x", "y", "z")
  pp <- s2pp(xyz)
  output <- data.frame(xyz, row.names = paste(1:2))

  expect_equal(as.data.frame(pp, format = "x,y,z"), output)
  expect_equal(as.data.frame(pp, format = "x,z,y"), output[,c(1,3,2)])
})

test_that("Coordinates in lat,long format are handled", {
  lat <- c(-45,45)
  lon <- c(0,90)
  xyz <- globe::ensure3d(data.frame(lon, lat), single = FALSE)
  pp <- s2pp(xyz)

  expect_equal(as.data.frame(pp), data.frame(long = lon, lat = lat))
  expect_equal(as.data.frame(pp, format = "lon,lat"), data.frame(lon = lon, lat = lat))
})
#
# test_that("Coordinates in bad format generate errors", {
#   expect_error(make_s2coords(1:3))
#   expect_error(make_s2coords(data.frame(x=1, lat=2)))
#   expect_error(make_s2coords(data.frame(x=1, y=2, lat=3))
#   expect_error(make_s2coords(as.list(1:4)))
#   expect_error(make_s2coords(cbind(1:2, 3:4)))
# })
