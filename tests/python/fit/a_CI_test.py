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
print('getting obs')
obs_df = fit.get_optical_obs(body_id, optical_obs_file, t_min_tdb, t_max_tdb, debias_lowres, deweight, eliminate, num_obs_per_night, verbose)
# obs_df = fit.add_radar_obs(obs_df, t_min_tdb, t_max_tdb, verbose)
# if add_gaia_obs:
#     gaia_dr = 'gaiafpr'
#     obs_df = fit.add_gaia_obs(obs_df, t_min_tdb, t_max_tdb, gaia_dr, verbose)

# # %% [markdown]
# # #### All we need to do now is initialize the OD simulation and run the filter.

# # %%
# n_iter_max = 10
# fit_sim = fit.FitSimulation(init_sol, obs_df, init_cov, n_iter_max=n_iter_max, de_kernel=de_kernel, nongrav_info=nongrav_info)

# # %%
# fit_sim.filter_lsq()

# # %% [markdown]
# # #### Let's print some summary statistics and plot some results.

# # %%
# fit_sim.print_summary()

# # %%
# fit_sim.plot_summary(auto_close=True)

# # %%
# fit_sim.iters[-1].plot_iteration_summary(title='Postfit Residuals', auto_close=True)

# # %%
# mean_0 = np.array(list(init_sol.values())[1:])
# cov_0 = init_cov
# mean_f = np.array(list(fit_sim.x_nom.values()))
# cov_f = fit_sim.covariance

# maha_dist_f, maha_dist_0, bhattacharya, bhatt_coeff = fit.get_similarity_stats(mean_0, cov_0, mean_f, cov_f)
# print(f'Mahalonobis distance between JPL and GRSS solution: {maha_dist_f:0.2f}')
# print(f'Mahalonobis distance between GRSS and JPL solution: {maha_dist_0:0.2f}')
# print(f'Bhattacharya distance between JPL and GRSS solution: {bhattacharya:0.4f}')
# print(f'Bhattacharya coefficient between JPL and GRSS solution: {bhatt_coeff:0.4f}')

# # %% [markdown]
# # #### Finally, we'll make sure the GRSS solution is statistically consistent with the JPL SBDB solution

# # %%
# assert maha_dist_f < 5.0
# assert maha_dist_0 < 5.0
# assert bhattacharya < 0.10
# assert bhatt_coeff > 0.90

# # %%



