# %% [markdown]
# # (12104) Chesley orbit determination test

# %% [markdown]
# #### Let's start by importing the necessary libraries

# %%
print('before importing grss.fit')
from grss import fit
print('imported grss.fit')
import numpy as np
np.set_printoptions(precision=40, linewidth=np.inf)

# %% [markdown]
# #### We'll then retrieve the cometary state of the asteroid (from JPL SBDB) plus any nongravitational accelerations acting on it.

# %%
body_id = '12104'
init_sol, init_cov, nongrav_info = fit.get_sbdb_info(body_id)
de_kernel = 440

# %% [markdown]
# #### Next, we'll retrieve the observations from different sources (MPC, JPL, Gaia Data Releases) and prepare them for the orbit determination process.

# %%
add_gaia_obs = True
optical_obs_file = None
t_min_tdb = None
t_max_tdb = None
debias_lowres = True
deweight = True
eliminate = False
num_obs_per_night = 4
verbose = True
accept_weights = False

print('getting obs')
from grss.fit.fit_optical import ades_column_types, _ades_mode_check, _ades_ast_cat_check, get_mpc_raw_data
if eliminate and deweight:
    raise ValueError('Cannot deweight and eliminate observations at the same time.')
if not eliminate and not deweight and verbose:
    print("WARNING: No deweighting or elimination scheme applied",
            "for observations during the same night.")
print('creating obs df')
# obs_df = create_optical_obs_df(body_id, optical_obs_file,
#                                 t_min_tdb, t_max_tdb, verbose)

import os
from astropy.time import Time
import numpy as np
import pandas as pd

if t_min_tdb is None:
    t_min_utc = -np.inf
else:
    t_min_utc = Time(t_min_tdb, format='mjd', scale='tdb').utc.mjd
if t_max_tdb is None:
    t_max_utc = np.inf
else:
    t_max_utc = Time(t_max_tdb, format='mjd', scale='tdb').utc.mjd
if optical_obs_file is None:
    obs_data = get_mpc_raw_data(body_id)
else:
    obs_data = optical_obs_file
    if not os.path.exists(obs_data):
        raise FileNotFoundError(f"File {obs_data} does not exist.")
# read in the data and add extra columns if not present
obs_df = pd.read_xml(obs_data, dtype=ades_column_types)

obs_times = Time(obs_df['obsTime'].to_list(), format='isot', scale='utc')
# obs_df['obsTimeMJD'] = obs_times.utc.mjd
# obs_df['obsTimeMJDTDB'] = obs_times.tdb.mjd
