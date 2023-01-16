#include "s2/s2cap.h"
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(s2)]]

// [[Rcpp::export]]
double cpp_s2_cap(double x, double y, double z, double height) {
  S2Point cent(x, y, z);
  S2Cap cap;
  S2Cap cap2 = cap.FromCenterHeight(cent, height);
  return cap2.height();
}
