#!/bin/bash

# get the earth orientation binary spice kernels and their comments if they are not already present
# latest earth pck
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.cmt
# historical earth pck
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_720101_070426.bpc
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_720101_070426.cmt
# predicted earth pck
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_200101_990628_predict.bpc
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_200101_990628_predict.cmt
