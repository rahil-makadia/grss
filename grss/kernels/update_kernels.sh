#!/bin/bash

# get the leap seconds kernel if it is not already present
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls

# get the earth orientation binary spice kernels and their comments if they are not already present
# latest earth pck
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.cmt
# historical earth pck
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_720101_230601.bpc
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_720101_230601.cmt
# predicted earth pck
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_200101_990825_predict.bpc
wget --no-clobber https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_200101_990825_predict.cmt

python update_kernels.py
