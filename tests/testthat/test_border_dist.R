library(spatstat.sphere)
context("Border distances.")

test_that("Bord distances are correct", {
  template <- data.frame(lng = c(1,-1,-1,1), lat = c(1,1,-1,-1))
  loop <- template*10
  poly <- s2polygon(loop)
  points <- data.frame(lng = c(0,5,10,30,180), lat = c(0,0,10,0,0))
  nearestborder <- data.frame(lon = c(0,10,10,10,10), lat=c(10,0,10,0,10))
  rslt <- s2dist(points, nearestborder)

  expect_equal(s2borderdist(points, poly), rslt)

  loops <- list(template*10, template*20)
  nearestborder <- data.frame(lon = c(0,10,10,20,20), lat=c(10,0,10,0,20))
  rslt <- s2dist(points, nearestborder)

  expect_equal(s2borderdist(points, s2polygon(loops)), rslt)
})
