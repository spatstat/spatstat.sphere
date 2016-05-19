# spatstatSphere: Point patterns on the sphere.

This repository is an effort to combine R code for point patterns on the sphere
and general geometrical computations on the sphere developed by:

- Adrian Baddeley
- Tom Lawrence
- Tuomas Rajala
- Ege Rubak

It depends on the R package `spatstat` (and possibly other packages) and it
could be considered a plugin to `spatstat` enabeling direct analysis of point
patterns on the sphere without having to project to a planar map.

It is currently **development mode** and backwards compatebility between versions
cannot be guaranteed.

Installation is most easily done with the `devtools` package:

```{r}
library(devtools)
install_github('spatstat/spatstatSphere')
```
