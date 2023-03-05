\begintext
   
    Custom meta-kernel that loads de441 planetary and big16 asteroid perturber
    ephemerides and naif0012 leapseconds kernel.

\begindata

    KERNELS_TO_LOAD = ( '../kernels/naif0012.tls',
                        '../kernels/planets_big16_de441_1950_2350.bsp',
                        '../kernels/earth_200101_990628_predict.bpc',
                        '../kernels/earth_720101_070426.bpc',
                        '../kernels/earth_000101_230525_230301.bpc')

\begintext