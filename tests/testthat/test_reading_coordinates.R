library(spatstat.sphere)
context("Reading coordinates.")

test_that("Coordinates in x,y,z format are recognised", {
  output <- matrix(1:6, ncol = 3, dimnames = list(NULL, c("x", "y", "z")))

  expect_equal(make_s2coords(list(1:2, 3:4, 5:6)), output)
  expect_equal(make_s2coords(matrix(1:6, ncol = 3)), output)
})

test_that("Coordinates in lat,long format are recognised", {
  lat <- c(-45,45)
  lon <- c(0,90)
  xyz <- globe::ensure3d(data.frame(lon, lat), single = FALSE)
  colnames(xyz) <- c("x", "y", "z")

  expect_equal(make_s2coords(cbind(lat=lat,lon=lon)), xyz)
  expect_equal(make_s2coords(data.frame(lat, lon)), xyz)
  expect_equal(make_s2coords(data.frame(long = lon, latitude = lat)), xyz)
  expect_equal(make_s2coords(data.frame(longitude = lon, latitude = lat)), xyz)
})

test_that("Coordinates in bad format generate errors", {
  expect_error(make_s2coords(1:3))
  expect_error(make_s2coords(data.frame(x=1, lat=2)))
  expect_error(make_s2coords(data.frame(x=1, y=2, lat=3)))
  expect_error(make_s2coords(as.list(1:4)))
  expect_error(make_s2coords(cbind(1:2, 3:4)))
})
