import numpy as np
import matplotlib.pyplot as plt

from . import prop
from .fit_utilities import *

__all__ = [ 'fitSimulation',
]

class iterationParams:
    def __init__(self, iter_number, x_nom, covariance, residuals, sigmas, weight_matrix):
        self.iter_number = iter_number
        self.x_nom = x_nom
        self.covariance = covariance
        self.residuals = residuals
        self.sigmas = sigmas
        self.weight_matrix = weight_matrix
        self.calculate_rms()
        self.calculate_chis()
        return None

    def calculate_rms(self):
        residual_arr = np.array(self.residuals)
        n_obs = len(self.residuals)
        self.unweighted_rms = float(np.sqrt(residual_arr.T @ residual_arr/n_obs))
        self.weighted_rms = float(np.sqrt(residual_arr.T @ self.weight_matrix @ residual_arr/n_obs))
        return None

    def calculate_chis(self):
        residual_arr = np.array(self.residuals)
        self.chi = residual_arr/self.sigmas
        self.chi_squared = np.sum(self.chi**2)
        self.reduced_chi_squared = self.chi_squared/(len(self.residuals)-len(self.x_nom))
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

class fitSimulation:
    def __init__(self, x_init, cov_init=None, obs_array_optical=None, observer_codes_optical=None, obs_array_radar=None, observer_codes_radar=None, n_iter=5, DEkernel=441, DEkernelPath=''):
        self.check_initial_solution(x_init, cov_init)
        self.check_input_observation_arrays(obs_array_optical, observer_codes_optical, obs_array_radar, observer_codes_radar)
        self.assemble_observation_arrays()
        self.n_iter = n_iter
        self.iters = []
        self.DEkernel = DEkernel
        self.DEkernelPath = DEkernelPath
        self.analytic_partials = False
        self.fixed_nongrav_params = {'a1': 0.0, 'a2': 0.0, 'a3': 0.0}
        return None

    def check_initial_solution(self, x_init, cov_init):
        if 't' not in x_init:
            raise ValueError("Must provide a time for the initial solution.")
        if any(key not in x_init for key in ("x", "y", "z", "vx", "vy", "vz")):
            raise ValueError("Must provide at least a full cartesian state for the initial solution.")
        for key in ["a1", "a2", "a3"]:
            if key in x_init and x_init[key] != 0.0:
                self.fit_nongrav = True
        self.t = x_init['t']
        x_init = x_init
        self.x_nom = {key: x_init[key] for key in x_init if key != 't'} # remove t for self.x_nom
        self.n_fit = len(self.x_nom)
        if cov_init.shape != (self.n_fit, self.n_fit):
            raise ValueError("Covariance matrix must be the same size as the number of fitted parameters.")
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
        sigmas = self.obs_array[:, 3:5].flatten()
        sigmas = sigmas[~np.isnan(sigmas)]
        self.sigmas = sigmas
        self.weight_matrix = np.diag(1/sigmas**2)
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

    def x_dict_to_cartesian_state(self, x_dict):
        x = x_dict['x']
        y = x_dict['y']
        z = x_dict['z']
        vx = x_dict['vx']
        vy = x_dict['vy']
        vz = x_dict['vz']
        pos = [x, y, z]
        vel = [vx, vy, vz]
        return pos, vel

    def x_dict_to_nongrav_params(self, x_dict):
        nongravParams = prop.NongravParamaters()
        a1 = x_dict['a1'] if 'a1' in x_dict.keys() else self.fixed_nongrav_params['a1']
        a2 = x_dict['a2'] if 'a2' in x_dict.keys() else self.fixed_nongrav_params['a2']
        a3 = x_dict['a3'] if 'a3' in x_dict.keys() else self.fixed_nongrav_params['a3']
        nongravParams.a1 = a1
        nongravParams.a2 = a2
        nongravParams.a3 = a3
        nongravParams.alpha = 1.0
        nongravParams.k = 0.0
        nongravParams.m = 2.0
        nongravParams.n = 0.0
        nongravParams.r0_au = 1.0
        return nongravParams

    def get_perturbed_state(self, key, fd_pert):
        x_plus = self.x_nom.copy()
        x_minus = self.x_nom.copy()
        fd_delta = self.x_nom[key]*fd_pert # fd_pert = finite difference perturbation to nominal state for calculating derivatives
        x_plus[key] = self.x_nom[key]+fd_delta
        pos_plus, vel_plus = self.x_dict_to_cartesian_state(x_plus)
        ngParams_plus = self.x_dict_to_nongrav_params(x_plus)
        x_minus[key] = self.x_nom[key]-fd_delta
        pos_minus, vel_minus = self.x_dict_to_cartesian_state(x_minus)
        ngParams_minus = self.x_dict_to_nongrav_params(x_minus)
        return pos_plus, vel_plus, ngParams_plus, pos_minus, vel_minus, ngParams_minus, fd_delta

    def get_perturbation_info(self):
        perturbation_info = []
        for key in self.x_nom:
            pert_result = self.get_perturbed_state(key, fd_pert=1e-3)
            perturbation_info.append(tuple(pert_result))
        return perturbation_info

    def assemble_and_propagate_bodies(self, perturbation_info):
        # get propagated states
        propSimPast, propSimFuture = self.get_propSims("orbit_fit_sim")
        # create nominal IntegBody object
        consts = prop.Constants()
        pos_nom, vel_nom = self.x_dict_to_cartesian_state(self.x_nom)
        ngParams_nom = self.x_dict_to_nongrav_params(self.x_nom)
        cov_nom = self.covariance
        integBody_nom = prop.IntegBody("integBody_nom", self.t, 0.0, 0.0, pos_nom, vel_nom, cov_nom, ngParams_nom, consts)
        # add the nominal IntegBody for the residuals
        if self.pastObsExist:
            propSimPast.add_integ_body(integBody_nom)
        if self.futureObsExist:
            propSimFuture.add_integ_body(integBody_nom)
        # add the perturbed IntegBodies for numerical derivatives
        if not self.analytic_partials and perturbation_info is not None:
            for i in range(self.n_fit):
                key = list(self.x_nom.keys())[i]
                pos_plus, vel_plus, ngParams_plus, pos_minus, vel_minus, ngParams_minus, fd_delta = perturbation_info[i]
                integBody_plus = prop.IntegBody(f"integBody_pert_{key}_plus", self.t, 0.0, 0.0, pos_plus, vel_plus, cov_nom, ngParams_plus, consts)
                integBody_minus = prop.IntegBody(f"integBody_pert_{key}_minus", self.t, 0.0, 0.0, pos_minus, vel_minus, cov_nom, ngParams_minus, consts)
                if self.pastObsExist:
                    propSimPast.add_integ_body(integBody_plus)
                    propSimPast.add_integ_body(integBody_minus)
                if self.futureObsExist:
                    propSimFuture.add_integ_body(integBody_plus)
                    propSimFuture.add_integ_body(integBody_minus)
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
        
        residuals = observed_obs - computed_obs
        residuals = residuals.flatten()
        residuals = residuals[~np.isnan(residuals)]
        return residuals

    def get_analytic_partials(self, propSimPast, propSimFuture):
        raise NotImplementedError("Analytic partials not yet implemented. Please use numeric partials.")

    def get_numeric_partials(self, propSimPast, propSimFuture, perturbation_info):
        partials = np.zeros((self.n_obs, self.n_fit))
        for i in range(self.n_fit):
            key = list(self.x_nom.keys())[i]
            pos_plus, vel_plus, ngParams_plus, pos_minus, vel_minus, ngParams_minus, fd_delta = perturbation_info[i]
            # get residuals for perturbed states
            residuals_plus = self.get_residuals(propSimPast, propSimFuture, integBodyIdx=2*i+1)
            residuals_minus = self.get_residuals(propSimPast, propSimFuture, integBodyIdx=2*i+2)
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
        # destroy propSim objects
        del propSimPast
        del propSimFuture
        return residuals, partials

    def add_iteration(self, iter_number, residuals):
        self.iters.append(iterationParams(iter_number, self.x_nom, self.covariance, residuals, self.sigmas, self.weight_matrix))
        return None
    
    def filter_lsq(self):
        print("Iteration\tUnweighted RMS\tWeighted RMS\tChi Squared\tReduced Chi Squared")
        for i in range(self.n_iter):
            # get residuals and partials
            b, a = self.get_residuals_and_partials()
            # add iteration
            self.add_iteration(i, b)
            print("%d\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f" % (i, self.iters[-1].unweighted_rms, self.iters[-1].weighted_rms, self.iters[-1].chi_squared, self.iters[-1].reduced_chi_squared))
            # get initial guess
            x0 = np.array(list(self.x_nom.values()))
            # get covariance
            cov = self.covariance
            # get state correction
            w = self.weight_matrix
            P = np.linalg.inv(a.T @ w @ a)
            dx = P @ a.T @ w @ b
            print("dx: ", dx)
            # get new state
            x = x0 + dx
            self.x_nom = dict(zip(self.x_nom.keys(), x))
            # get new covariance
            self.covariance = P
        return None

    def print_summary(self):
        print("Iteration\tUnweighted RMS\tWeighted RMS\tChi-squared\tReduced Chi-squared")
        for iteration in self.iters:
            print("%d\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f" % (iteration.iter_number, iteration.unweighted_rms, iteration.weighted_rms, iteration.chi_squared, iteration.reduced_chi_squared))
        
        print("Summary of the orbit fit calculations:")
        print("====================================================")
        print(f"Number of iterations: {len(self.iters)}")
        print(f"Final RMS unweighted: {self.iters[-1].unweighted_rms}")
        print(f"Final RMS weighted: {self.iters[-1].weighted_rms}")
        print(f"Final chi-squared: {self.iters[-1].chi_squared}")
        print(f"Final reduced chi-squared: {self.iters[-1].reduced_chi_squared}")
        print("====================================================")
        return None

    def plot_summary(self):
        plt.figure(figsize=(25,10), dpi=100)
        plt.subplot(2, 2, 1)
        plt.semilogy([iteration.iter_number for iteration in self.iters], [iteration.unweighted_rms for iteration in self.iters], label=f"Unweighted RMS={self.iters[-1].unweighted_rms:.3e}")
        plt.xticks(np.arange(0, len(self.iters), 1))
        plt.xlabel("Iteration #")
        plt.ylabel("Unweighted RMS")
        plt.legend()

        plt.subplot(2, 2, 2)
        plt.semilogy([iteration.iter_number for iteration in self.iters], [iteration.weighted_rms for iteration in self.iters], label=f"Weighted RMS={self.iters[-1].weighted_rms:.3e}")
        plt.xticks(np.arange(0, len(self.iters), 1))
        plt.xlabel("Iteration #")
        plt.ylabel("Weighted RMS")
        plt.legend()

        plt.subplot(2, 2, 3)
        plt.semilogy([iteration.iter_number for iteration in self.iters], [iteration.chi_squared for iteration in self.iters], label=f"$\chi^2$={self.iters[-1].chi_squared:.3e}")
        plt.xticks(np.arange(0, len(self.iters), 1))
        plt.xlabel("Iteration #")
        plt.ylabel("$\chi^2$")
        plt.legend()

        plt.subplot(2, 2, 4)
        plt.semilogy([iteration.iter_number for iteration in self.iters], [iteration.reduced_chi_squared for iteration in self.iters], label=f"Reduced $\chi^2$={self.iters[-1].reduced_chi_squared:.3e}")
        plt.xticks(np.arange(0, len(self.iters), 1))
        plt.xlabel("Iteration #")
        plt.ylabel("Reduced $\chi^2$")
        plt.legend()
        plt.show()
        return None
