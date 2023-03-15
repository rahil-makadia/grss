import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice

from . import prop
from .fit_utilities import *

__all__ = [ 'fitSimulation',
]

class iterationParams:
    def __init__(self, iter_number, x_nom, covariance, residuals, obs_array, observer_codes):
        self.iter_number = iter_number
        self.x_nom = x_nom
        self.covariance = covariance
        variance = np.sqrt(np.diag(self.covariance))
        self.variance = dict(zip([f'var_{k}' for k in self.x_nom.keys()], variance))
        self.residuals = residuals
        self.obs_array = obs_array
        self.observer_codes = observer_codes
        self.sigmas = obs_array[:, 3:5]
        self.weight_matrix = np.diag(1/self.flatten_and_clean(self.sigmas)**2)
        self.calculate_rms()
        self.calculate_chis()
        return None

    def flatten_and_clean(self, arr):
        arr = arr.flatten()
        arr = arr[~np.isnan(arr)]
        return arr

    def calculate_rms(self):
        residual_arr = self.flatten_and_clean(self.residuals)
        n_obs = len(residual_arr)
        self.unweighted_rms = float(np.sqrt(residual_arr.T @ residual_arr/n_obs))
        self.weighted_rms = float(np.sqrt(residual_arr.T @ self.weight_matrix @ residual_arr/n_obs))
        return None

    def calculate_chis(self):
        residual_arr = self.flatten_and_clean(self.residuals)
        n_obs = len(residual_arr)
        n_fit = len(self.x_nom)
        sigmas = self.flatten_and_clean(self.sigmas)
        self.chi = residual_arr/sigmas
        self.chi_squared = np.sum(self.chi**2)
        self.reduced_chi_squared = self.chi_squared/(n_obs-n_fit)
        return None
    
    # adapted from https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html
    def scatter_hist(self, x, y, ax, ax_histx, ax_histy, size, show_logarithmic):
        color = 'C0'
        fill = False
        nbins = 100
        # no labels
        ax_histx.tick_params(   top=True, labeltop=False,
                                bottom=True, labelbottom=False,
                                left=True, labelleft=True,
                                right=True, labelright=False,
                                direction='in')
        ax_histy.tick_params(   left=True, labelleft=False,
                                right=True, labelright=False,
                                top=True, labeltop=False,
                                bottom=True, labelbottom=True,
                                direction='in')
        ax.tick_params( left=True, labelleft=True,
                        right=True, labelright=False,
                        top=True, labeltop=False,
                        bottom=True, labelbottom=True,
                        direction='in')
        # the scatter plot:
        ax.scatter(x, y, s=size, c=color)
        # now determine nice limits:
        if show_logarithmic:
            hist, bins = np.histogram(x, bins=nbins)
            bins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        else:
            bins = nbins
        ax_histx.hist(x, bins=bins, orientation='vertical', color=color, edgecolor=color, linewidth=0.75, fill=fill, histtype='step')
        ax_histy.hist(y, bins=bins, orientation='horizontal', color=color, edgecolor=color, linewidth=0.75, fill=fill, histtype='step')
        return None

    def plot_residuals(self, t_arr_optical, t_arr_radar, ra_residuals, dec_residuals, ra_cosdec_residuals, delay_residuals, doppler_residuals, radar_scale, markersize, show_logarithmic):
        # sourcery skip: extract-duplicate-method
        fig = plt.figure(figsize=(21,6), dpi=150)
        iter_string = f'Iteration {self.iter_number} (prefit)' if self.iter_number == 0 else f'Iteration {self.iter_number}'
        plt.suptitle(iter_string, y=0.95)
        gs = fig.add_gridspec(1, 3, width_ratios=(1,1,1))
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(t_arr_optical, ra_residuals, '.', label='RA', markersize=markersize)
        ax1.plot(t_arr_optical, dec_residuals, '.', label='Dec', markersize=markersize)
        ax1.legend()
        ax1.set_xlabel('MJD [UTC]')
        ax1.set_ylabel('Residuals, O-C [arcsec]')
        ax1.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic: ax1.set_yscale('log')

        ax2 = gs[0,1].subgridspec(2, 2, width_ratios=(4,1), height_ratios=(1,4), wspace=0.05, hspace=0.05)
        ax2main = fig.add_subplot(ax2[1,0])
        ax2histx = fig.add_subplot(ax2[0,0], sharex=ax2main)
        ax2histy = fig.add_subplot(ax2[1,1], sharey=ax2main)
        self.scatter_hist(ra_cosdec_residuals, dec_residuals, ax2main, ax2histx, ax2histy, markersize, show_logarithmic)
        ax2main.set_xlabel('RA cos(Dec) Residuals, O-C [arcsec]')
        ax2main.set_ylabel('Dec Residuals, O-C [arcsec]')
        ax2main.grid(True, which='both', axis='both', alpha=0.2, zorder=-100)
        if show_logarithmic: ax2main.set_yscale('log'); ax2main.set_xscale('log')
        
        ax3 = fig.add_subplot(gs[0, 2])
        ax3.plot(t_arr_radar, delay_residuals, '.', mfc='C2', mec='C2', label='Delay', markersize=radar_scale*markersize)
        ax3.plot(t_arr_radar, doppler_residuals, '.', mfc='C3', mec='C3', label='Doppler', markersize=radar_scale*markersize)
        ax3.legend()
        ax3.set_xlabel('MJD [UTC]')
        ax3.set_ylabel('Residuals, O-C [$\mu$s, Hz]')
        ax3.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic: ax3.set_yscale('log')
        plt.show()

    def plot_chi(self, t_arr_optical, t_arr_radar, ra_chi, dec_chi, delay_chi, doppler_chi, ra_chi_squared, dec_chi_squared, delay_chi_squared, doppler_chi_squared, sigma_limit, radar_scale, markersize, show_logarithmic):
        # plot chi values
        plt.figure(figsize=(21,6), dpi=150)
        iter_string = f'Iteration {self.iter_number} (prefit)' if self.iter_number == 0 else f'Iteration {self.iter_number}'
        plt.suptitle(f'{iter_string}. Chi Squared: RA={np.sum(ra_chi_squared):.2f}, Dec={np.sum(dec_chi_squared):.2f}, Delay={np.sum(delay_chi_squared):.2f}, Doppler={np.sum(doppler_chi_squared):.2f}', y=0.95)
        plt.subplot(1,2,1)
        plt.plot(t_arr_optical, ra_chi, '.', markersize=markersize, label='RA')
        plt.plot(t_arr_optical, dec_chi, '.', markersize=markersize, label='Dec')
        plt.plot(t_arr_radar, delay_chi, '.', mfc='C2', mec='C2', markersize=radar_scale*markersize, label='Delay')
        plt.plot(t_arr_radar, doppler_chi, '.', mfc='C3', mec='C3', markersize=radar_scale*markersize, label='Doppler')
        plt.axhline(-sigma_limit, c='khaki', linestyle='--', alpha=1.0, label=f'$\pm{sigma_limit:.0f}\sigma$')
        plt.axhline(sigma_limit, c='khaki', linestyle='--', alpha=1.0)
        plt.axhline(-2*sigma_limit, c='red', linestyle='--', alpha=0.5, label=f'$\pm{2*sigma_limit:.0f}\sigma$')
        plt.axhline(2*sigma_limit, c='red', linestyle='--', alpha=0.5)
        plt.legend()
        plt.xlabel('MJD [UTC]')
        plt.ylabel('$\chi$, (O-C)/$\sigma$ $[\cdot]$')
        plt.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic: plt.yscale('log')

        plt.subplot(1,2,2)
        plt.plot(t_arr_optical, ra_chi_squared, '.', markersize=markersize, label='RA')
        plt.plot(t_arr_optical, dec_chi_squared, '.', markersize=markersize, label='Dec')
        plt.plot(t_arr_radar, delay_chi_squared, '.', mfc='C2', mec='C2', markersize=radar_scale*markersize, label='Delay')
        plt.plot(t_arr_radar, doppler_chi_squared, '.', mfc='C3', mec='C3', markersize=radar_scale*markersize, label='Doppler')
        plt.legend()
        plt.xlabel('MJD [UTC]')
        plt.ylabel('$\chi^2$, (O-C)$^2/\sigma^2$ $[\cdot]$')
        plt.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic: plt.yscale('log')
        plt.show()

    def plot_iteration_summary(self, show_logarithmic=False):
        markersize = 3
        sigma_limit = 3
        radar_scale = 3
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.labelsize'] = 12
        
        obs_array = self.obs_array
        residuals = self.residuals
        # opticalIdx is where neither the ra nor dec residuals are NaN
        opticalIdx = np.where(~np.isnan(residuals[:, 0]) & ~np.isnan(residuals[:, 1]))[0]
        # radarIdx is where either the ra or dec residuals are NaN
        radarIdx = np.where(np.isnan(residuals[:, 0]) | np.isnan(residuals[:, 1]))[0]
        
        t_arr_optical = obs_array[opticalIdx, 0]
        ra_obs = obs_array[opticalIdx, 1]
        dec_obs = obs_array[opticalIdx, 2]
        ra_noise = obs_array[opticalIdx, 3]
        dec_noise = obs_array[opticalIdx, 4]
        ra_residuals = residuals[opticalIdx, 0]
        dec_residuals = residuals[opticalIdx, 1]
        ra_computed = ra_obs - ra_residuals
        dec_computed = dec_obs - dec_residuals
        ra_cosdec_residuals = ra_obs*np.cos(dec_obs) - ra_computed*np.cos(dec_computed)
        ra_chi = ra_residuals/ra_noise
        dec_chi = dec_residuals/dec_noise
        ra_chi_squared = ra_chi**2
        dec_chi_squared = dec_chi**2
        ra_residuals *= 180/np.pi*3600
        dec_residuals *= 180/np.pi*3600
        ra_cosdec_residuals *= 180/np.pi*3600
        ra_residuals = np.abs(ra_residuals) if show_logarithmic else ra_residuals
        dec_residuals = np.abs(dec_residuals) if show_logarithmic else dec_residuals
        ra_cosdec_residuals = np.abs(ra_cosdec_residuals) if show_logarithmic else ra_cosdec_residuals
        ra_chi = np.abs(ra_chi) if show_logarithmic else ra_chi
        dec_chi = np.abs(dec_chi) if show_logarithmic else dec_chi
        ra_chi_squared = np.abs(ra_chi_squared) if show_logarithmic else ra_chi_squared
        dec_chi_squared = np.abs(dec_chi_squared) if show_logarithmic else dec_chi_squared
        
        t_arr_radar = obs_array[radarIdx, 0]
        delay_obs = obs_array[radarIdx, 1]
        doppler_obs = obs_array[radarIdx, 2]
        delay_noise = obs_array[radarIdx, 3]
        doppler_noise = obs_array[radarIdx, 4]
        delay_residuals = residuals[radarIdx, 0]
        doppler_residuals = residuals[radarIdx, 1]
        delay_computed = delay_obs - delay_residuals
        doppler_computed = doppler_obs - doppler_residuals
        delay_chi = delay_residuals/delay_noise
        doppler_chi = doppler_residuals/doppler_noise
        delay_chi_squared = delay_chi**2
        doppler_chi_squared = doppler_chi**2
        delay_residuals = np.abs(delay_residuals) if show_logarithmic else delay_residuals
        doppler_residuals = np.abs(doppler_residuals) if show_logarithmic else doppler_residuals
        delay_chi = np.abs(delay_chi) if show_logarithmic else delay_chi
        doppler_chi = np.abs(doppler_chi) if show_logarithmic else doppler_chi
        
        self.plot_residuals(t_arr_optical, t_arr_radar, ra_residuals, dec_residuals, ra_cosdec_residuals, delay_residuals, doppler_residuals, radar_scale, markersize, show_logarithmic)
        self.plot_chi(t_arr_optical, t_arr_radar, ra_chi, dec_chi, delay_chi, doppler_chi, ra_chi_squared, dec_chi_squared, delay_chi_squared, doppler_chi_squared, sigma_limit, radar_scale, markersize, show_logarithmic)
        return None

class fitSimulation:
    def __init__(self, x_init, cov_init=None, obs_array_optical=None, observer_codes_optical=None, obs_array_radar=None, observer_codes_radar=None, n_iter=5, DEkernel=441, DEkernelPath='', radius=0.0, nongravInfo=None, events=None):
        self.check_initial_solution(x_init, cov_init)
        self.check_input_observation_arrays(obs_array_optical, observer_codes_optical, obs_array_radar, observer_codes_radar)
        self.assemble_observation_arrays()
        self.n_iter = n_iter
        self.iters = []
        self.DEkernel = DEkernel
        self.DEkernelPath = DEkernelPath
        self.analytic_partials = False
        self.fixed_propSim_params = {'a1': 0.0, 'a2': 0.0, 'a3': 0.0, 'alpha': 1.0, 'k': 0.0, 'm': 2.0, 'n': 0.0, 'r0_au': 1.0, 'radius': radius, 'mass': 0.0}
        for key in nongravInfo:
            self.fixed_propSim_params[key] = nongravInfo[key]
        if events is not None:
            self.fixed_propSim_params['events'] = events
            self.fit_events = True
        else:
            self.fixed_propSim_params['events'] = []
            self.fit_events = False
        return None

    def check_initial_solution(self, x_init, cov_init):
        if 't' not in x_init:
            raise ValueError("Must provide a time for the initial solution.")
        if all(key in x_init for key in ("x", "y", "z", "vx", "vy", "vz")):
            self.fit_cartesian = True
            self.fit_cometary = False
        elif all(key in x_init for key in ("e", "q", "tp", "om", "w", "i")):
            self.fit_cartesian = False
            self.fit_cometary = True
        else:
            raise ValueError("Must provide at least a full cartesian or cometary state for the initial solution.")
        for key in ["a1", "a2", "a3"]:
            if key in x_init and x_init[key] != 0.0:
                self.fit_nongrav = True
        self.t = x_init['t']
        self.x_init = x_init
        self.x_nom = {key: x_init[key] for key in x_init if key != 't'} # remove t for self.x_nom
        self.n_fit = len(self.x_nom)
        if cov_init.shape != (self.n_fit, self.n_fit):
            raise ValueError("Covariance matrix must be the same size as the number of fitted parameters.")
        self.covariance_init = cov_init
        self.covariance = cov_init
        return None

    def check_input_observation_arrays(self, obs_array_optical, observer_codes_optical, obs_array_radar, observer_codes_radar):
        self.fit_optical = False
        self.fit_radar = False
        if obs_array_optical is None and obs_array_radar is None:
            raise ValueError("Must provide at least one observation array (optical or radar).")
        
        if obs_array_optical is not None and observer_codes_optical is None:
            raise ValueError("Must provide observer codes for optical observations.")
        if obs_array_optical is None and observer_codes_optical is not None:
            raise ValueError("Must provide optical observations for observer codes.")
        if obs_array_optical is not None and observer_codes_optical is not None:
            if len(obs_array_optical) != len(observer_codes_optical):
                raise ValueError("Optical observation array and observer code array must be the same length.")
            else:
                self.fit_optical = True
        if obs_array_radar is not None and observer_codes_radar is None:
            raise ValueError("Must provide observer codes for radar observations.")
        if obs_array_radar is None and observer_codes_radar is not None:
            raise ValueError("Must provide radar observations for observer codes.")
        if obs_array_radar is not None and observer_codes_radar is not None:
            if len(obs_array_radar) != len(observer_codes_radar):
                raise ValueError("Radar observation array and observer code array must be the same length.")
            else:
                self.fit_radar = True
        self.obs_array_optical = obs_array_optical
        self.observer_codes_optical = observer_codes_optical
        self.obs_array_radar = obs_array_radar
        self.observer_codes_radar = observer_codes_radar
        return None

    def flatten_and_clean(self, arr):
        arr = arr.flatten()
        arr = arr[~np.isnan(arr)]
        return arr

    def assemble_observation_arrays(self):
        if self.fit_optical and self.fit_radar:
            self.merge_observation_arrays()
        elif self.fit_optical:
            self.obs_array = self.obs_array_optical
            self.observer_codes = self.observer_codes_optical
        elif self.fit_radar:
            self.obs_array = self.obs_array_radar
            self.observer_codes = self.observer_codes_radar
        self.n_obs = np.count_nonzero(~np.isnan(self.obs_array[:, 1:3])) # number of observations is the number of non-nan values in the second and third columns of the observation array
        self.sigmas = self.obs_array[:, 3:5]
        self.weight_matrix = np.diag(1/self.flatten_and_clean(self.sigmas)**2)
        self.pastObsIdx = np.where(self.obs_array[:, 0] < self.t)[0]
        self.pastObsExist = len(self.pastObsIdx) > 0
        self.futureObsIdx = np.where(self.obs_array[:, 0] >= self.t)[0]
        self.futureObsExist = len(self.futureObsIdx) > 0
        return None

    def merge_observation_arrays(self):
        # merge optical and radar arrays
        obs_array = np.vstack((self.obs_array_optical, self.obs_array_radar))
        observer_codes = self.observer_codes_optical + self.observer_codes_radar
        # sort by time
        sort_idx = np.argsort(obs_array[:,0])
        self.obs_array = obs_array[sort_idx]
        self.observer_codes = tuple(np.array(observer_codes, dtype=tuple)[sort_idx])
        return None

    def get_propSimPast(self, name, tEvalUTC, evalApparentState, convergedLightTime, observerInfo):
        tEvalPast = self.obs_array[self.pastObsIdx, 0]
        tfPast = np.min(tEvalPast) - 1
        propSimPast = prop.propSimulation(name, self.t, self.DEkernel, self.DEkernelPath)
        propSimPast.set_integration_parameters(tfPast, tEvalPast, tEvalUTC, evalApparentState, convergedLightTime, observerInfo)
        return propSimPast

    def get_propSimFuture(self, name, tEvalUTC, evalApparentState, convergedLightTime, observerInfo):
        tEvalFuture = self.obs_array[self.futureObsIdx, 0]
        tfFuture = np.max(tEvalFuture) + 1
        propSimFuture = prop.propSimulation(name, self.t, self.DEkernel, self.DEkernelPath)
        propSimFuture.set_integration_parameters(tfFuture, tEvalFuture, tEvalUTC, evalApparentState, convergedLightTime, observerInfo)
        return propSimFuture

    def get_propSims(self, name):
        tEvalUTC = True
        evalApparentState = True
        convergedLightTime = True
        observerInfo = np.array(get_observer_info(self.observer_codes), dtype=tuple)
        observerInfoPast = tuple(observerInfo[self.pastObsIdx])
        observerInfoFuture = tuple(observerInfo[self.futureObsIdx])
        propSimPast = None
        propSimFuture = None
        if self.pastObsExist:
            propSimPast = self.get_propSimPast(f"{name}_past", tEvalUTC, evalApparentState, convergedLightTime, observerInfoPast)
        if self.futureObsExist:
            propSimFuture = self.get_propSimFuture(f"{name}_future", tEvalUTC, evalApparentState, convergedLightTime, observerInfoFuture)
        return propSimPast, propSimFuture

    def x_dict_to_state(self, x_dict):
        # sourcery skip: extract-method
        if self.fit_cartesian:
            x = x_dict['x']
            y = x_dict['y']
            z = x_dict['z']
            vx = x_dict['vx']
            vy = x_dict['vy']
            vz = x_dict['vz']
            state = [x, y, z, vx, vy, vz]
        elif self.fit_cometary:
            e = x_dict['e']
            q = x_dict['q']
            tp = x_dict['tp']
            om = x_dict['om']
            w = x_dict['w']
            i = x_dict['i']
            cometary_elements = [e, q, tp, om*np.pi/180.0, w*np.pi/180.0, i*np.pi/180.0]
            state = cometary_elements
        else:
            raise ValueError("fit_cartesian or fit_cometary must be True")
        return state

    def x_dict_to_nongrav_params(self, x_dict):
        nongravParams = prop.NongravParamaters()
        a1 = x_dict['a1'] if 'a1' in x_dict.keys() else self.fixed_propSim_params['a1']
        a2 = x_dict['a2'] if 'a2' in x_dict.keys() else self.fixed_propSim_params['a2']
        a3 = x_dict['a3'] if 'a3' in x_dict.keys() else self.fixed_propSim_params['a3']
        nongravParams.a1 = a1
        nongravParams.a2 = a2
        nongravParams.a3 = a3
        nongravParams.alpha = 1.0
        nongravParams.k = 0.0
        nongravParams.m = 2.0
        nongravParams.n = 0.0
        nongravParams.r0_au = 1.0
        return nongravParams

    def x_dict_to_events(self, x_dict):
        events = []
        for i in range(len(self.fixed_propSim_params['events'])):
            event = self.fixed_propSim_params['events'][i]
            event[1] = x_dict[f"dvx{i}"] if f"dvx{i}" in x_dict.keys() else event[1]
            event[2] = x_dict[f"dvy{i}"] if f"dvy{i}" in x_dict.keys() else event[2]
            event[3] = x_dict[f"dvz{i}"] if f"dvz{i}" in x_dict.keys() else event[3]
            event[4] = x_dict[f"multiplier{i}"] if f"multiplier{i}" in x_dict.keys() else event[4]
            events.append(tuple(event))
        return events

    def check_and_add_events(self, propSimPast, propSimFuture, integBody, events):
        for event in events:
            t_event = event[0]
            dvx = event[1]
            dvy = event[2]
            dvz = event[3]
            multiplier = event[4]
            if t_event < self.t:
                propSimPast.add_event(integBody, t_event, [dvx, dvy, dvz], multiplier)
            else:
                propSimFuture.add_event(integBody, t_event, [dvx, dvy, dvz], multiplier)
        return propSimPast, propSimFuture

    def get_perturbed_state(self, key):
        if self.fit_cartesian:
            if key in ['x', 'y', 'z']:
                fd_pert = 1e-7
            elif key in ['vx', 'vy', 'vz']:
                fd_pert = 1e-10
        elif self.fit_cometary:
            if key in ['q']:
                fd_pert = 5e-5
            elif key in ['e']:
                fd_pert = 5e-4
            elif key in ['tp']:
                fd_pert = 1e-5
            elif key in ['om']:
                fd_pert = 1e-4
            elif key in ['w']:
                fd_pert = 1e-5
            elif key in ['i']:
                fd_pert = 1e-6
        if key in ['a1', 'a2', 'a3']:
            fd_pert = 1e-3
        if key[:10] == 'multiplier':
            fd_pert = 1e-2
        if key[:3] in ['dvx', 'dvy', 'dvz']:
            fd_pert = 1e-4

        x_plus = self.x_nom.copy()
        x_minus = self.x_nom.copy()
        fd_delta = self.x_nom[key]*fd_pert # fd_pert = finite difference perturbation to nominal state for calculating derivatives
        x_plus[key] = self.x_nom[key]+fd_delta
        state_plus = self.x_dict_to_state(x_plus)
        ngParams_plus = self.x_dict_to_nongrav_params(x_plus)
        events_plus = self.x_dict_to_events(x_plus)
        x_minus[key] = self.x_nom[key]-fd_delta
        state_minus = self.x_dict_to_state(x_minus)
        ngParams_minus = self.x_dict_to_nongrav_params(x_minus)
        events_minus = self.x_dict_to_events(x_minus)
        return state_plus, ngParams_plus, events_plus, state_minus, ngParams_minus, events_minus, fd_delta

    def get_perturbation_info(self):
        perturbation_info = []
        for key in self.x_nom:
            pert_result = self.get_perturbed_state(key)
            perturbation_info.append(tuple(pert_result))
        return perturbation_info

    def assemble_and_propagate_bodies(self, perturbation_info):
        # get propagated states
        propSimPast, propSimFuture = self.get_propSims("orbit_fit_sim")
        # create nominal IntegBody object
        consts = prop.Constants()
        state_nom = self.x_dict_to_state(self.x_nom)
        ngParams_nom = self.x_dict_to_nongrav_params(self.x_nom)
        events_nom = self.x_dict_to_events(self.x_nom)
        cov_nom = self.covariance
        if self.fit_cartesian:
            integBody_nom = prop.IntegBody("integBody_nom", self.t, self.fixed_propSim_params['mass'], self.fixed_propSim_params['radius'], state_nom[:3], state_nom[3:6], cov_nom, ngParams_nom, consts)
        elif self.fit_cometary:
            integBody_nom = prop.IntegBody(self.DEkernelPath, "integBody_nom", self.t, self.fixed_propSim_params['mass'], self.fixed_propSim_params['radius'], state_nom, cov_nom, ngParams_nom, consts)
        # add the nominal IntegBody for the residuals
        if self.pastObsExist:
            propSimPast.add_integ_body(integBody_nom)
        if self.futureObsExist:
            propSimFuture.add_integ_body(integBody_nom)
        propSimPast, propSimFuture = self.check_and_add_events(propSimPast, propSimFuture, integBody_nom, events_nom)
        # add the perturbed IntegBodies for numerical derivatives
        if not self.analytic_partials and perturbation_info is not None:
            for i in range(self.n_fit):
                key = list(self.x_nom.keys())[i]
                state_plus, ngParams_plus, events_plus, state_minus, ngParams_minus, events_minus, fd_delta = perturbation_info[i]
                if self.fit_cartesian:
                    integBody_plus = prop.IntegBody(f"integBody_pert_{key}_plus", self.t, self.fixed_propSim_params['mass'], self.fixed_propSim_params['radius'], state_plus[:3], state_plus[3:6], cov_nom, ngParams_plus, consts)
                    integBody_minus = prop.IntegBody(f"integBody_pert_{key}_minus", self.t, self.fixed_propSim_params['mass'], self.fixed_propSim_params['radius'], state_minus[:3], state_minus[3:6], cov_nom, ngParams_minus, consts)
                elif self.fit_cometary:
                    integBody_plus = prop.IntegBody(self.DEkernelPath, f"integBody_pert_{key}_plus", self.t, self.fixed_propSim_params['mass'], self.fixed_propSim_params['radius'], state_plus, cov_nom, ngParams_plus, consts)
                    integBody_minus = prop.IntegBody(self.DEkernelPath, f"integBody_pert_{key}_minus", self.t, self.fixed_propSim_params['mass'], self.fixed_propSim_params['radius'], state_minus, cov_nom, ngParams_minus, consts)
                if self.pastObsExist:
                    propSimPast.add_integ_body(integBody_plus)
                    propSimPast.add_integ_body(integBody_minus)
                if self.futureObsExist:
                    propSimFuture.add_integ_body(integBody_plus)
                    propSimFuture.add_integ_body(integBody_minus)
                propSimPast, propSimFuture = self.check_and_add_events(propSimPast, propSimFuture, integBody_plus, events_plus)
                propSimPast, propSimFuture = self.check_and_add_events(propSimPast, propSimFuture, integBody_minus, events_minus)
        if self.pastObsExist:
            propSimPast.preprocess()
            propSimPast.integrate()
        if self.futureObsExist:
            propSimFuture.preprocess()
            propSimFuture.integrate()
        return propSimPast, propSimFuture

    def get_residuals(self, propSimPast, propSimFuture, integBodyIdx):
        if self.pastObsExist and self.futureObsExist:
            apparent_states_past = np.array(propSimPast.xIntegEval)
            apparent_states_future = np.array(propSimFuture.xIntegEval)
            apparent_states = np.vstack((apparent_states_past, apparent_states_future))
            radar_observations_past = np.array(propSimPast.radarObsEval)
            radar_observations_future = np.array(propSimFuture.radarObsEval)
            radar_observations = np.vstack((radar_observations_past, radar_observations_future))
        elif self.pastObsExist:
            apparent_states_past = np.array(propSimPast.xIntegEval)
            apparent_states = apparent_states_past
            radar_observations_past = np.array(propSimPast.radarObsEval)
            radar_observations = radar_observations_past
        elif self.futureObsExist:
            apparent_states_future = np.array(propSimFuture.xIntegEval)
            apparent_states = apparent_states_future
            radar_observations_future = np.array(propSimFuture.radarObsEval)
            radar_observations = radar_observations_future
        integ_body_start_col = 6*integBodyIdx
        integ_body_end_col = 6*integBodyIdx+6
        apparent_states = apparent_states[:,integ_body_start_col:integ_body_end_col]
        radar_observations = radar_observations[:,integBodyIdx]

        observed_obs = self.obs_array[:, 1:3]
        computed_obs = np.nan*np.ones_like(observed_obs)
        observerInfo = get_observer_info(self.observer_codes)
        for i in range(len(self.obs_array)):
            obs_info_len = len(observerInfo[i])
            if obs_info_len == 4:
                computed_obs[i, :] = get_radec(apparent_states[i])
            elif obs_info_len == 8: # delay measurement
                computed_obs[i, 0] = radar_observations[i]
            elif obs_info_len == 9: # dopper measurement
                computed_obs[i, 1] = radar_observations[i]
        return observed_obs - computed_obs

    def get_analytic_partials(self, propSimPast, propSimFuture):
        raise NotImplementedError("Analytic partials not yet implemented. Please use numeric partials.")

    def get_numeric_partials(self, propSimPast, propSimFuture, perturbation_info):
        partials = np.zeros((self.n_obs, self.n_fit))
        for i in range(self.n_fit):
            key = list(self.x_nom.keys())[i]
            state_plus, ngParams_plus, events_plus, state_minus, ngParams_minus, events_minus, fd_delta = perturbation_info[i]
            # get residuals for perturbed states
            residuals_plus = self.get_residuals(propSimPast, propSimFuture, integBodyIdx=2*i+1)
            residuals_minus = self.get_residuals(propSimPast, propSimFuture, integBodyIdx=2*i+2)
            residuals_plus = self.flatten_and_clean(residuals_plus)
            residuals_minus = self.flatten_and_clean(residuals_minus)
            # get partials
            partials[:, i] = (residuals_plus - residuals_minus)/(2*fd_delta)
        return partials

    def get_partials(self, propSimPast, propSimFuture, perturbation_info):
        return self.get_analytic_partials(propSimPast, propSimFuture) if self.analytic_partials else self.get_numeric_partials(propSimPast, propSimFuture, perturbation_info)

    def get_residuals_and_partials(self):
        perturbation_info = None if self.analytic_partials else self.get_perturbation_info()
        propSimPast, propSimFuture = self.assemble_and_propagate_bodies(perturbation_info)
        # get residuals
        residuals = self.get_residuals(propSimPast, propSimFuture, integBodyIdx=0)
        # get partials
        partials = self.get_partials(propSimPast, propSimFuture, perturbation_info)
        return residuals, partials

    def add_iteration(self, iter_number, residuals):
        self.iters.append(iterationParams(iter_number, self.x_nom, self.covariance, residuals, self.obs_array, self.observer_codes))
        return None

    def filter_lsq(self):
        spice.furnsh(self.DEkernelPath)
        print("Iteration\t\tUnweighted RMS\t\tWeighted RMS\t\tChi-squared\t\tReduced Chi-squared")
        for i in range(self.n_iter):
            # get residuals and partials
            residuals, a = self.get_residuals_and_partials()
            if i == 0:
                # add prefit iteration
                prefit_residuals = residuals.copy()
                self.add_iteration(0, prefit_residuals)
            b = self.flatten_and_clean(residuals)
            # get initial guess
            x0 = np.array(list(self.x_nom.values()))
            # get covariance
            cov = self.covariance
            # get state correction
            w = self.weight_matrix
            atwa = a.T @ w @ a
            P = np.linalg.inv(atwa)
            atwb = a.T @ w @ b
            dx = P @ atwb
            # get new state
            x = x0 - dx
            self.x_nom = dict(zip(self.x_nom.keys(), x))
            # get new covariance
            self.covariance = P
            # add iteration
            self.add_iteration(i+1, residuals)
            print("%d%s\t\t\t%.3f\t\t\t%.3f\t\t\t%.3f\t\t\t%.3f" % (self.iters[-1].iter_number, "", self.iters[-1].unweighted_rms, self.iters[-1].weighted_rms, self.iters[-1].chi_squared, self.iters[-1].reduced_chi_squared))
        spice.unload(self.DEkernelPath)
        return None

    def print_summary(self, iter_idx=-1):
        data = self.iters[iter_idx]
        print(f"Summary of the orbit fit calculations at iteration {data.iter_number} (of {self.n_iter}):")
        print("====================================================")
        print(f"RMS unweighted: {data.unweighted_rms}")
        print(f"RMS weighted: {data.weighted_rms}")
        print(f"chi-squared: {data.chi_squared}")
        print(f"reduced chi-squared: {data.reduced_chi_squared}")
        print("====================================================")
        print(f"t: MJD {self.t} TDB")
        print("Fitted Variable\t\tInitial Value\t\t\tUncertainty\t\t\tFitted Value\t\t\tUncertainty\t\t\tChange\t\t\t\tChange (sigma)")
        init_variance = np.sqrt(np.diag(self.covariance_init))
        final_variance = np.sqrt(np.diag(self.covariance))
        init_sol = self.iters[0].x_nom
        final_sol = data.x_nom
        with np.errstate(divide='ignore'):
            for i, key in enumerate(init_sol.keys()):
                if key[:10] == 'multiplier':
                    print(f"{key}\t\t{init_sol[key]:.11e}\t\t{init_variance[i]:.11e}\t\t{final_sol[key]:.11e}\t\t{final_variance[i]:.11e}\t\t{final_sol[key]-init_sol[key]:+.11e}\t\t{(final_sol[key]-init_sol[key])/init_variance[i]:+.3f}")
                else:
                    print(f"{key}\t\t\t{init_sol[key]:.11e}\t\t{init_variance[i]:.11e}\t\t{final_sol[key]:.11e}\t\t{final_variance[i]:.11e}\t\t{final_sol[key]-init_sol[key]:+.11e}\t\t{(final_sol[key]-init_sol[key])/init_variance[i]:+.3f}")
        return None

    def plot_summary(self):
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.labelsize'] = 12
        ticks = np.arange(1, self.n_iter+1, 1)
        start_idx = 1
        iters_for_plot = self.iters[start_idx:]
        plt.figure(figsize=(21,10), dpi=150)
        plt.subplot(2, 2, 1)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot], [iteration.unweighted_rms for iteration in iters_for_plot], label=f"Final Unweighted RMS={iters_for_plot[-1].unweighted_rms:.3e}")
        plt.xticks(ticks)
        plt.xlabel("Iteration #")
        plt.ylabel("Unweighted RMS")
        plt.legend()
        
        plt.subplot(2, 2, 2)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot], [iteration.weighted_rms for iteration in iters_for_plot], label=f"Final Weighted RMS={iters_for_plot[-1].weighted_rms:.3e}")
        plt.xticks(ticks)
        plt.xlabel("Iteration #")
        plt.ylabel("Weighted RMS")
        plt.legend()
        
        plt.subplot(2, 2, 3)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot], [iteration.chi_squared for iteration in iters_for_plot], label=f"Final $\chi^2$={iters_for_plot[-1].chi_squared:.3e}")
        plt.xticks(ticks)
        plt.xlabel("Iteration #")
        plt.ylabel("$\chi^2$")
        plt.legend()
        
        plt.subplot(2, 2, 4)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot], [iteration.reduced_chi_squared for iteration in iters_for_plot], label=f"Final Reduced $\chi^2$={iters_for_plot[-1].reduced_chi_squared:.3e}")
        plt.xticks(ticks)
        plt.xlabel("Iteration #")
        plt.ylabel("Reduced $\chi^2$")
        plt.legend()
        plt.show()
        return None
