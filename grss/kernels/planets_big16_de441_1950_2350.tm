\begintext
   
    Custom meta-kernel that loads de441 planetary and big16 asteroid perturber
    ephemerides and latest leapseconds kernel.

    The path_values variable is written by the update_kernels.py script.

\begindata

    PATH_VALUES  = ( placeholder )

    PATH_SYMBOLS = ( 'GRSS' )

    KERNELS_TO_LOAD = ( '$GRSS/latest_leapseconds.tls',
                        '$GRSS/pck00011.tpc',
                        '$GRSS/earth_200101_990825_predict.bpc',
                        '$GRSS/earth_720101_230601.bpc',
                        '$GRSS/earth_latest_high_prec.bpc' )

\begintext
