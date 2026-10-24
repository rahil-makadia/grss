Complete Python helper API
==========================

All Python-defined functions and classes are explicitly exported from their
implementation modules. Fitting and propagation functions are also available
through ``grss.fit`` and ``grss.prop``. Utility and download functions are
available through ``grss.utils``, ``grss.kernels`` and ``grss.debias``.
These are in-process Python APIs, not HTTP endpoints. Native C++ bindings
continue to be provided by ``grss.libgrss``; this catalog covers Python code.

Every formerly private helper has a public alias. The underscore-prefixed
spelling remains callable and exported. Both names refer to the same function;
class helpers remain methods of the relevant class, with normal instance
binding. Python constructors retain their normal ``__init__`` spelling.

Examples::

    from grss.fit import propagate_two_body, stumpff
    from grss.prop import reconstruct_one_log_file
    from grss.utils import connected_to_internet
    from grss.kernels import download_kernels
    from grss.debias import download_debiasing_data

    # Legacy private imports also work.
    from grss.fit import _propagate

    # Methods operate on an appropriately initialized instance.
    mean_weights, covariance_weights = sigma_points.get_weights()
    solution_dicts = sigma_points.get_points_dict()

Requirements and side effects
-----------------------------

* ``FitSimulation`` stages require an initialized fitter and the ordering
  described in :doc:`fit`. Setup/parsing methods rebuild caches and may append
  simulated observations; propagation methods require native GRSS objects.
  Weight construction follows observation parsing. Residual/derivative methods
  require propagated states. Outlier rejection requires a current covariance.
  These stages are not interchangeable stateless operations.
* ``IterationParams.assemble_info`` refreshes cached plotting arrays;
  ``calculate_chis`` updates statistics from those arrays. Plot helpers require
  prepared arrays and Matplotlib axes, and mutate those axes.
* ``SigmaPoints.get_weights``, ``merwe_weights`` and ``julier_weights`` return
  weights from the instance's configured dimensions/scaling. ``get_points``
  returns points/offsets; ``get_points_dict`` reads the cached ``sigma_points``.
  If parameters/covariance change, refresh the cached weights and points before
  reconstructing a transformed distribution. ``reconstruct`` requires exactly
  the corresponding ordered set of transformed sigma points.
* Propagation log reconstruction requires GRSS-format logs or parsed row
  dictionaries. Object reconstruction requires the native extension.
  ``handle_one_cloned_sim`` returns the numerical argument list consumed by
  the native parallel propagator; it does not itself propagate an orbit.
* The Gauss component functions use central-body-relative states, geometric
  unit directions and consistent units. ``propagate_two_body`` returns position
  only; ``gauss_angular_residuals`` returns residuals/ranges/maximum angle;
  ``refine_gauss_candidate`` returns a corrected state and convergence data.
  These numerical helpers do not obtain ephemerides or apply observational
  corrections. See :doc:`fit` and their docstrings for conventions.
* Network helpers perform requests when explicitly called. ``download_kernels``
  may fetch multi-GB files and refresh the latest Earth PCK. ``download_codes_file``
  replaces local code tables; already imported dictionaries are not refreshed.
  ``connected_to_internet`` tests completion of a HEAD request, not the health
  of every data source. ``initialize`` retains normal package initialization
  side effects, including data downloads and the OpenMP environment setting.
* Importing the downloader submodules does not download or extract files.
  ``download_debiasing_data`` explicitly fetches missing tables and unpacks
  existing archives, replacing same-named files. Unrelated files are kept;
  archives containing traversal paths or links are rejected. The package's
  normal top-level import still invokes its existing ``initialize`` workflow.

The catalog below lists signatures and compatibility names. Existing
NumPy-style docstrings remain the detailed input/return reference; private
method aliases are included in the generated class reference as well.

Callable catalog
----------------

grss.debias.get_debiasing_data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    grss.debias.get_debiasing_data.check_and_download(url, local_path, note)
    grss.debias.get_debiasing_data.download_debiasing_data(destination=None)
    grss.debias.get_debiasing_data.extract_debiasing_archive(archive_path, destination)

grss.fit.fit_iod
~~~~~~~~~~~~~~~~

::

    grss.fit.fit_iod.GaussIODSolution(epoch: float, state: numpy.ndarray, ranges: numpy.ndarray, root: float, residual: float, iterations: int) -> None
    grss.fit.fit_iod.GaussIODSolution.to_fit_dict(self, central_body_state)
    grss.fit.fit_iod.GaussIODSolution.__init__(self, epoch: float, state: numpy.ndarray, ranges: numpy.ndarray, root: float, residual: float, iterations: int) -> None
    grss.fit.fit_iod.GaussIODSolution.__repr__(self)
    grss.fit.fit_iod.GaussIODSolution.__eq__(self, other)
    grss.fit.fit_iod.gauss_iod(times, ra, dec, observer_positions, mu, *, max_iter=50, tol=1e-10)
    grss.fit.fit_iod._stumpff(z)
    grss.fit.fit_iod._propagate(state, dt, mu)
    grss.fit.fit_iod._residuals(state, dt, observers, los, mu)
    grss.fit.fit_iod._refine(state, dt, observers, los, mu, max_iter, tol)
    grss.fit.fit_iod.stumpff(z)
    grss.fit.fit_iod.propagate_two_body(state, dt, mu)
    grss.fit.fit_iod.gauss_angular_residuals(state, dt, observers, los, mu)
    grss.fit.fit_iod.refine_gauss_candidate(state, dt, observers, los, mu, max_iter, tol)
    grss.fit.fit_iod.universal_kepler_equation(anomaly, radius, radial_dot, alpha, dt, mu)

grss.fit.fit_optical
~~~~~~~~~~~~~~~~~~~~

::

    grss.fit.fit_optical.get_mpc_raw_data(tdes)
    grss.fit.fit_optical.validate_ades_mode(df)
    grss.fit.fit_optical.flag_unrecognized_ades_catalogs(df)
    grss.fit.fit_optical.create_optical_obs_df(body_id, optical_obs_file=None, t_min_tdb=None, t_max_tdb=None, verbose=False)
    grss.fit.fit_optical.add_psv_obs(psv_obs_file, obs_df, t_min_tdb=None, t_max_tdb=None, verbose=False)
    grss.fit.fit_optical.get_gaia_query_results(body_id, release)
    grss.fit.fit_optical.add_gaia_obs(obs_df, t_min_tdb=None, t_max_tdb=None, gaia_dr='gaiafpr', verbose=False)
    grss.fit.fit_optical.debias_obs(r_asc_vals, dec_vals, epoch_vals, catalog, biasdf, nside)
    grss.fit.fit_optical.apply_debiasing_scheme(obs_df, lowres, verbose)
    grss.fit.fit_optical.get_packed_prog_id(code)
    grss.fit.fit_optical.get_unpacked_prog_id(code)
    grss.fit.fit_optical.apply_station_weight_rules(group, obs_df, cols, verbose)
    grss.fit.fit_optical.apply_weighting_scheme(obs_df, verbose)
    grss.fit.fit_optical.deweight_obs(obs_df, eff_obs_per_night, verbose)
    grss.fit.fit_optical.eliminate_obs(obs_df, max_obs_per_night, verbose)
    grss.fit.fit_optical.get_optical_obs(body_id, optical_obs_file=None, t_min_tdb=None, t_max_tdb=None, debias_lowres=True, deweight=True, eliminate=False, num_obs_per_night=5, verbose=False, accept_weights=False)
    grss.fit.fit_optical._ades_mode_check(df)
    grss.fit.fit_optical._ades_ast_cat_check(df)
    grss.fit.fit_optical._get_gaia_query_results(body_id, release)

grss.fit.fit_radar
~~~~~~~~~~~~~~~~~~

::

    grss.fit.fit_radar.get_radar_raw_data(tdes)
    grss.fit.fit_radar.add_radar_obs(obs_df, t_min_tdb=None, t_max_tdb=None, verbose=False)

grss.fit.fit_simulation
~~~~~~~~~~~~~~~~~~~~~~~

::

    grss.fit.fit_simulation.flatten_valid_observations(arr)
    grss.fit.fit_simulation.IterationParams(iter_number, x_nom, covariance, obs, rms_u, rms_w, chi_sq)
    grss.fit.fit_simulation.IterationParams.__init__(self, iter_number, x_nom, covariance, obs, rms_u, rms_w, chi_sq)
    grss.fit.fit_simulation.IterationParams._calculate_chis(self)
    grss.fit.fit_simulation.IterationParams._assemble_info(self)
    grss.fit.fit_simulation.IterationParams._scatter_hist(self, x_val, y_val, axis, ax_histx, ax_histy, size, show_logarithmic)
    grss.fit.fit_simulation.IterationParams._plot_residuals(self, t_arr, ra_residuals, dec_residuals, delay_residuals, doppler_residuals, radar_scale, markersize, show_logarithmic, title, savefig, figname, auto_close)
    grss.fit.fit_simulation.IterationParams._plot_chi(self, t_arr, ra_chi, dec_chi, delay_chi, doppler_chi, ra_chi_squared, dec_chi_squared, delay_chi_squared, doppler_chi_squared, plot_chi_squared, sigma_limit, radar_scale, markersize, show_logarithmic, title, savefig, figname, auto_close)
    grss.fit.fit_simulation.IterationParams.plot_iteration_summary(self, show_logarithmic=False, title=None, plot_chi_squared=False, savefig=False, figname=None, auto_close=False)
    grss.fit.fit_simulation.IterationParams.calculate_chis(self)
    grss.fit.fit_simulation.IterationParams.assemble_info(self)
    grss.fit.fit_simulation.IterationParams.scatter_hist(self, x_val, y_val, axis, ax_histx, ax_histy, size, show_logarithmic)
    grss.fit.fit_simulation.IterationParams.plot_residuals(self, t_arr, ra_residuals, dec_residuals, delay_residuals, doppler_residuals, radar_scale, markersize, show_logarithmic, title, savefig, figname, auto_close)
    grss.fit.fit_simulation.IterationParams.plot_chi(self, t_arr, ra_chi, dec_chi, delay_chi, doppler_chi, ra_chi_squared, dec_chi_squared, delay_chi_squared, doppler_chi_squared, plot_chi_squared, sigma_limit, radar_scale, markersize, show_logarithmic, title, savefig, figname, auto_close)
    grss.fit.fit_simulation.FitSimulation(x_init, obs_df, cov_init=None, n_iter_max=10, de_kernel=441, nongrav_info=None, events=None, simulated_obs=None)
    grss.fit.fit_simulation.FitSimulation.__init__(self, x_init, obs_df, cov_init=None, n_iter_max=10, de_kernel=441, nongrav_info=None, events=None, simulated_obs=None)
    grss.fit.fit_simulation.FitSimulation.check_initial_solution(self, x_init, cov_init)
    grss.fit.fit_simulation.FitSimulation.prepare_priors(self)
    grss.fit.fit_simulation.FitSimulation.add_simulated_obs(self)
    grss.fit.fit_simulation.FitSimulation.parse_observation_arrays(self, obs_df)
    grss.fit.fit_simulation.FitSimulation.compute_obs_weights(self)
    grss.fit.fit_simulation.FitSimulation.get_prop_sim_past(self, name, t_eval_utc, eval_apparent_state, converged_light_time, observer_info)
    grss.fit.fit_simulation.FitSimulation.get_prop_sim_future(self, name, t_eval_utc, eval_apparent_state, converged_light_time, observer_info)
    grss.fit.fit_simulation.FitSimulation.get_prop_sims(self)
    grss.fit.fit_simulation.FitSimulation.solution_to_state(self, x_dict)
    grss.fit.fit_simulation.FitSimulation.solution_to_nongrav_params(self, x_dict)
    grss.fit.fit_simulation.FitSimulation.solution_to_events(self, x_dict)
    grss.fit.fit_simulation.FitSimulation.check_and_add_events(self, prop_sim_past, prop_sim_future, integ_body, events)
    grss.fit.fit_simulation.FitSimulation.get_perturbed_state(self, key)
    grss.fit.fit_simulation.FitSimulation.get_perturbation_info(self)
    grss.fit.fit_simulation.FitSimulation.assemble_and_propagate_bodies(self, perturbation_info)
    grss.fit.fit_simulation.FitSimulation.get_computed_obs(self, prop_sim_past, prop_sim_future, integ_body_idx)
    grss.fit.fit_simulation.FitSimulation.inflate_uncertainties(self, prop_sim_past, prop_sim_future)
    grss.fit.fit_simulation.FitSimulation.get_analytic_partials(self, prop_sim_past, prop_sim_future)
    grss.fit.fit_simulation.FitSimulation.get_numeric_partials(self, prop_sim_past, prop_sim_future, perturbation_info)
    grss.fit.fit_simulation.FitSimulation.get_partials(self, prop_sim_past, prop_sim_future, perturbation_info)
    grss.fit.fit_simulation.FitSimulation.compute_residuals_and_partials(self)
    grss.fit.fit_simulation.FitSimulation.compute_fit_statistics(self, partials, residuals, start_rejecting=False)
    grss.fit.fit_simulation.FitSimulation.add_iteration(self, iter_number, rms_u, rms_w, chi_sq)
    grss.fit.fit_simulation.FitSimulation.check_convergence(self, delta_x)
    grss.fit.fit_simulation.FitSimulation.solve_state_correction(self, partials, residuals)
    grss.fit.fit_simulation.FitSimulation.filter_lsq(self, verbose=True)
    grss.fit.fit_simulation.FitSimulation.print_summary(self, iter_idx=None, out_stream=<_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>)
    grss.fit.fit_simulation.FitSimulation.plot_summary(self, auto_close=False)
    grss.fit.fit_simulation.FitSimulation.save(self, filename)
    grss.fit.fit_simulation.FitSimulation._check_initial_solution(self, x_init, cov_init)
    grss.fit.fit_simulation.FitSimulation._add_simulated_obs(self)
    grss.fit.fit_simulation.FitSimulation._parse_observation_arrays(self, obs_df)
    grss.fit.fit_simulation.FitSimulation._compute_obs_weights(self)
    grss.fit.fit_simulation.FitSimulation._get_prop_sim_past(self, name, t_eval_utc, eval_apparent_state, converged_light_time, observer_info)
    grss.fit.fit_simulation.FitSimulation._get_prop_sim_future(self, name, t_eval_utc, eval_apparent_state, converged_light_time, observer_info)
    grss.fit.fit_simulation.FitSimulation._get_prop_sims(self)
    grss.fit.fit_simulation.FitSimulation._check_and_add_events(self, prop_sim_past, prop_sim_future, integ_body, events)
    grss.fit.fit_simulation.FitSimulation._assemble_and_propagate_bodies(self, perturbation_info)
    grss.fit.fit_simulation.FitSimulation._get_computed_obs(self, prop_sim_past, prop_sim_future, integ_body_idx)
    grss.fit.fit_simulation.FitSimulation._inflate_uncertainties(self, prop_sim_past, prop_sim_future)
    grss.fit.fit_simulation.FitSimulation._get_analytic_partials(self, prop_sim_past, prop_sim_future)
    grss.fit.fit_simulation.FitSimulation._get_numeric_partials(self, prop_sim_past, prop_sim_future, perturbation_info)
    grss.fit.fit_simulation.FitSimulation._get_partials(self, prop_sim_past, prop_sim_future, perturbation_info)
    grss.fit.fit_simulation.FitSimulation._add_iteration(self, iter_number, rms_u, rms_w, chi_sq)
    grss.fit.fit_simulation.FitSimulation._check_convergence(self, delta_x)
    grss.fit.fit_simulation.FitSimulation._check_priors(self)
    grss.fit.fit_simulation.FitSimulation._x_dict_to_state(self, x_dict)
    grss.fit.fit_simulation.FitSimulation._x_dict_to_nongrav_params(self, x_dict)
    grss.fit.fit_simulation.FitSimulation._x_dict_to_events(self, x_dict)
    grss.fit.fit_simulation.FitSimulation._get_perturbed_state(self, key)
    grss.fit.fit_simulation.FitSimulation._get_perturbation_info(self)
    grss.fit.fit_simulation.FitSimulation._get_residuals_and_partials(self)
    grss.fit.fit_simulation.FitSimulation._get_rms_and_reject_outliers(self, partials, residuals, start_rejecting=False)
    grss.fit.fit_simulation.FitSimulation._get_lsq_state_correction(self, partials, residuals)
    grss.fit.fit_simulation.FitSimulation.flatten_and_clean(self, arr)
    grss.fit.fit_simulation.FitSimulation._flatten_and_clean(self, arr)

grss.fit.fit_utils
~~~~~~~~~~~~~~~~~~

::

    grss.fit.fit_utils.mjd2et(mjd)
    grss.fit.fit_utils.et2mjd(ephem_time)
    grss.fit.fit_utils.parallax_to_lat_lon_alt(lon, rho_cos_lat, rho_sin_lat)
    grss.fit.fit_utils.get_mpc_observatory_info()
    grss.fit.fit_utils.get_radar_codes_dict()
    grss.fit.fit_utils.get_codes_dict()
    grss.fit.fit_utils.get_observer_info(obs_df)
    grss.fit.fit_utils.get_sbdb_raw_data(tdes)
    grss.fit.fit_utils.get_sbdb_elems(tdes, cov_elems=True)
    grss.fit.fit_utils.get_sbdb_info(tdes)
    grss.fit.fit_utils.get_similarity_stats(sol_1, cov_1, sol_2, cov_2)

grss.kernels.get_kernels
~~~~~~~~~~~~~~~~~~~~~~~~

::

    grss.kernels.get_kernels.check_and_download(url, local_path, note, force=False)
    grss.kernels.get_kernels.download_kernels(destination=None)

grss.prop.prop_parallel
~~~~~~~~~~~~~~~~~~~~~~~

::

    grss.prop.prop_parallel.parallel_propagate(ref_sol, ref_nongrav, ref_sim, clones, num_threads, reconstruct=False)
    grss.prop.prop_parallel.cluster_ca_or_impacts(full_list, max_duration=45, central_body=399)
    grss.prop.prop_parallel.reconstruct_all_log_files(log_dir)
    grss.prop.prop_parallel._handle_one_cloned_sim(sol, ref_nongrav)
    grss.prop.prop_parallel._reconstruct_one_log_file(log_file)
    grss.prop.prop_parallel._reconstruct_ca_and_impact(lines)
    grss.prop.prop_parallel._reconstruct_ca_params(row)
    grss.prop.prop_parallel._reconstruct_impact_params(row)
    grss.prop.prop_parallel.handle_one_cloned_sim(sol, ref_nongrav)
    grss.prop.prop_parallel.reconstruct_one_log_file(log_file)
    grss.prop.prop_parallel.reconstruct_ca_and_impact(lines)
    grss.prop.prop_parallel.reconstruct_ca_params(row)
    grss.prop.prop_parallel.reconstruct_impact_params(row)

grss.prop.prop_unscented
~~~~~~~~~~~~~~~~~~~~~~~~

::

    grss.prop.prop_unscented.SigmaPoints(x_dict, cov, sp_type, alpha=None, beta=None, kappa=None, sqrt_func=None)
    grss.prop.prop_unscented.SigmaPoints.__init__(self, x_dict, cov, sp_type, alpha=None, beta=None, kappa=None, sqrt_func=None)
    grss.prop.prop_unscented.SigmaPoints._get_weights(self)
    grss.prop.prop_unscented.SigmaPoints._merwe_weights(self)
    grss.prop.prop_unscented.SigmaPoints._julier_weights(self)
    grss.prop.prop_unscented.SigmaPoints._get_points(self)
    grss.prop.prop_unscented.SigmaPoints._get_points_dict(self)
    grss.prop.prop_unscented.SigmaPoints.reconstruct(self, transformed_sigma_points)
    grss.prop.prop_unscented.SigmaPoints.get_weights(self)
    grss.prop.prop_unscented.SigmaPoints.merwe_weights(self)
    grss.prop.prop_unscented.SigmaPoints.julier_weights(self)
    grss.prop.prop_unscented.SigmaPoints.get_points(self)
    grss.prop.prop_unscented.SigmaPoints.get_points_dict(self)

grss.prop.prop_utils
~~~~~~~~~~~~~~~~~~~~

::

    grss.prop.prop_utils.plot_solar_system(axis, xy_plane, alpha=0.25, lstyle='--')
    grss.prop.prop_utils.get_scale_factor(body_id)
    grss.prop.prop_utils.plot_ca_summary(prop_sim, flyby_body, central_body='Earth', x_min=None, x_max=None, y_min=0.1, y_max=None)
    grss.prop.prop_utils.mean_cov_to_ellipse(mean, cov, n_std, bplane_type, print_ellipse_params, units)
    grss.prop.prop_utils.partials_to_ellipse(ca, au2units, n_std, print_ellipse_params, units, analytic_info)
    grss.prop.prop_utils.data_to_ellipse(x_data, y_data, n_std, bplane_type, print_ellipse_params, units, sigma_points)
    grss.prop.prop_utils.days_to_dhms(days)
    grss.prop.prop_utils.plot_bplane(ca_list, plot_offset=False, scale_coords=False, n_std=3, units_km=False, equal_axis=True, print_ellipse_params=False, show_central_body=True, sigma_points=None, analytic_info=None)
    grss.prop.prop_utils.plot_single_bplane(axis, x_coord, y_coord, ellipse_data, bplane_type, focus_factor, show_central_body, plot_offset, scale_coords, central_body_radius, units, equal_axis)
    grss.prop.prop_utils.plot_earth_impact(impact_list, print_ellipse_params=False, sigma_points=None, zoom_size=5000.0, plot_ellipse=False, save_name=None)
    grss.prop.prop_utils.compare_earth_impact(coord1, coord2, sig1, sig2, labels, zoom_size=5000.0, save_name=None)
    grss.prop.prop_utils._add_planet_xy(axis, dist, color, alpha, lstyle)
    grss.prop.prop_utils._add_planet_xz(axis, dist, obliq, color, alpha, lstyle)
    grss.prop.prop_utils._get_bplane_labels(bplane_type, plot_offset, units)
    grss.prop.prop_utils.add_planet_xy(axis, dist, color, alpha, lstyle)
    grss.prop.prop_utils.add_planet_xz(axis, dist, obliq, color, alpha, lstyle)
    grss.prop.prop_utils.get_bplane_labels(bplane_type, plot_offset, units)

grss.utils
~~~~~~~~~~

::

    grss.utils.initialize()
    grss.utils._connected_to_internet(url='http://www.google.com/', timeout=30)
    grss.utils._download_codes_file()
    grss.utils.check_and_get_codes_file()
    grss.utils.connected_to_internet(url='http://www.google.com/', timeout=30)
    grss.utils.download_codes_file()
