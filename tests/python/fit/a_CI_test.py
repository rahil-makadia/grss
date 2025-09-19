from astropy.time import Time

str_times = ['1976-02-28T22:23:00.096Z', '1976-02-29T22:11:00.384Z', '1996-11-26T20:10:21.792Z', '1996-11-26T20:29:28.320Z', '1996-11-26T20:48:33.120Z', '1996-11-26T20:54:49.824Z']
# save str_times to file for debugging
obs_times = Time(str_times, format='isot', scale='utc')
# obs_df['obsTimeMJD'] = obs_times.utc.mjd
# obs_df['obsTimeMJDTDB'] = obs_times.tdb.mjd
