"""Simulation classes for the Python GRSS orbit determination code"""
# pylint: disable=no-name-in-module, no-member, too-many-lines, useless-return
import sys
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time

from .. import libgrss
from ..utils import default_kernel_path, grss_path
from .fit_utils import get_observer_info, get_similarity_stats
from .fit_ades import special_codes

__all__ = [ 'FitSimulation',
]

class IterationParams:
    """
    Class for storing the iteration parameters for a single iteration of the
    orbit determination process. It is also used for plotting the residuals
    and chi-squared values for each iteration.
    """
    def __init__(self, iter_number, x_nom, covariance, obs, rms_u, rms_w, chi_sq):
        """
        Constructor for the IterationParams class

        Parameters
        ----------
        iter_number : int
            Iteration number
        x_nom : dict
            Dictionary of nominal state vector values at the current iteration
        covariance : array
            Covariance matrix at the current iteration
        obs : pandas DataFrame
            Observation data for the orbit fit
        rms_u : float
            Unweighted RMS of the residuals.
        rms_w : float
            Weighted RMS of the residuals.
        chi_sq : float
            Chi-squared value for the residuals
        """
        self.iter_number = iter_number
        self.x_nom = x_nom
        self.covariance = covariance
        variance = np.sqrt(np.diag(self.covariance))
        self.variance = dict(zip([f'var_{k}' for k in self.x_nom.keys()], variance))
        self.obs = obs.copy()
        force_del_idx = self.obs[self.obs['selAst'] == 'd'].index
        auto_del_idx = self.obs[self.obs['selAst'] == 'D'].index
        self.rejected_idx = np.concatenate((force_del_idx, auto_del_idx))
        self.accepted_idx = np.setdiff1d(np.arange(len(self.obs)), self.rejected_idx)
        self.unweighted_rms = rms_u
        self.weighted_rms = rms_w
        self._assemble_info()
        accepted_obs = self.obs.iloc[self.accepted_idx]
        all_obs = (
            accepted_obs.ra.tolist() +
            accepted_obs.dec.tolist() +
            accepted_obs.delay.tolist() +
            accepted_obs.doppler.tolist()
        )
        n_obs = np.sum(~np.isnan(all_obs))
        n_fit = covariance.shape[0]
        self.chi_squared = chi_sq
        self.reduced_chi_squared = self.chi_squared/(n_obs-n_fit)
        return None

    def _calculate_chis(self):
        """
        Calculates the chi and chi-squared values for the residuals

        Returns
        -------
        None : NoneType
            None
        """
        chi = np.zeros((len(self.obs), 4))
        chi[:, 0] = self.all_info['ra_chi']
        chi[:, 1] = self.all_info['dec_chi']
        chi[:, 2] = self.all_info['delay_chi']
        chi[:, 3] = self.all_info['doppler_chi']
        # n_obs is number of non-nan entries in self.chi
        n_obs = np.sum(~np.isnan(chi[self.accepted_idx, :]))
        n_fit = self.covariance.shape[0]
        self.chi_squared = np.nansum(chi[self.accepted_idx, :]**2)
        self.reduced_chi_squared = self.chi_squared/(n_obs-n_fit)
        return None

    def _assemble_info(self):
        """
        Assembles the information for the iteration into a dictionary

        Returns
        -------
        None : NoneType
            None
        """
        # sourcery skip: extract-duplicate-method
        self.all_info = {
            'ra_res': self.obs['resRA'].values.copy(),
            'dec_res': self.obs['resDec'].values.copy(),
            'ra_obs': self.obs['ra'].values.copy(),
            'dec_obs': self.obs['dec'].values.copy(),
            'ra_noise': self.obs['sigRA'].values.copy(),
            'dec_noise': self.obs['sigDec'].values.copy(),
            'delay_res': self.obs['resDelay'].values.copy(),
            'doppler_res': self.obs['resDoppler'].values.copy(),
            'delay_obs': self.obs['delay'].values.copy(),
            'doppler_obs': self.obs['doppler'].values.copy(),
            'delay_noise': self.obs['sigDelay'].values.copy(),
            'doppler_noise': self.obs['sigDoppler'].values.copy(),
        }
        self.all_info['ra_chi'] = self.all_info['ra_res']/self.all_info['ra_noise']
        self.all_info['dec_chi'] = self.all_info['dec_res']/self.all_info['dec_noise']
        self.all_info['ra_chi_squared'] = self.all_info['ra_chi']**2
        self.all_info['dec_chi_squared'] = self.all_info['dec_chi']**2
        self.all_info['delay_chi'] = self.all_info['delay_res']/self.all_info['delay_noise']
        self.all_info['doppler_chi'] = self.all_info['doppler_res']/self.all_info['doppler_noise']
        self.all_info['delay_chi_squared'] = self.all_info['delay_chi']**2
        self.all_info['doppler_chi_squared'] = self.all_info['doppler_chi']**2
        return None

    # adapted from https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html
    def _scatter_hist(self, x_val, y_val, axis, ax_histx, ax_histy, size, show_logarithmic):
        """
        Create a scatter plot with histograms on the right and top 4x4 grid

        Parameters
        ----------
        x_val : vector
            x values to plot
        y_val : vector
            y values to plot
        axis : matplotlib axis
            axis to plot the scatter points on
        ax_histx : matplotlib axis
            axis to plot the x histogram on
        ax_histy : matplotlib axis
            axis to plot the y histogram on
        size : float
            size of the scatter points
        show_logarithmic : bool
            Flag to show the logarithmic scale on the histograms

        Returns
        -------
        None : NoneType
            None
        """
        color = 'C0'
        fill = False
        nbins = 100
        # no labels
        ax_histx.tick_params(top=True, labeltop=False,
                                bottom=True, labelbottom=False,
                                left=True, labelleft=True,
                                right=True, labelright=False,
                                direction='in')
        ax_histy.tick_params(left=True, labelleft=False,
                                right=True, labelright=False,
                                top=True, labeltop=False,
                                bottom=True, labelbottom=True,
                                direction='in')
        axis.tick_params(left=True, labelleft=True,
                            right=True, labelright=False,
                            top=True, labeltop=False,
                            bottom=True, labelbottom=True,
                            direction='in')
        # the scatter plot:
        axis.scatter(x_val, y_val, s=size, c=color)
        # now determine nice limits:
        if show_logarithmic:
            _, bins = np.histogram(x_val, bins=nbins)
            bins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        else:
            bins = nbins
        ax_histx.hist(x_val, bins=bins, orientation='vertical', color=color, edgecolor=color,
                        linewidth=0.75, fill=fill, histtype='step')
        ax_histy.hist(y_val, bins=bins, orientation='horizontal', color=color, edgecolor=color,
                        linewidth=0.75, fill=fill, histtype='step')
        return None

    def _plot_residuals(self, t_arr, ra_residuals, dec_residuals,
                        delay_residuals, doppler_residuals, radar_scale, markersize,
                        show_logarithmic, title, savefig, figname, auto_close):
        """
        Plot the residuals

        Parameters
        ----------
        t_arr : vector
            time array
        ra_residuals : vector
            right ascension residuals
        dec_residuals : vector
            declination residuals
        delay_residuals : vector
            delay residuals
        doppler_residuals : vector
            doppler residuals
        radar_scale : float
            scale factor for the radar scatter points
        markersize : float
            size of the scatter points
        show_logarithmic : bool
            Flag to show the logarithmic scale on the plot
        title : str
            title of the plot
        savefig : bool
            Flag to save the figure
        figname : str
            Name of the figure to save
        auto_close : bool, optional
            Flag to automatically close the figure

        Returns
        -------
        None : NoneType
            None
        """
        # sourcery skip: extract-duplicate-method
        rejected_idx = self.rejected_idx
        fig = plt.figure(figsize=(15,6), dpi=150)
        if self.iter_number == 0:
            iter_string = f'Iteration {self.iter_number} (prefit)'
        else:
            iter_string = f'Iteration {self.iter_number}'
        iter_string = title if title is not None else iter_string
        plt.suptitle(iter_string, y=0.95)
        grid_spec = fig.add_gridspec(1, 3, width_ratios=(1,1,1))
        ax1 = fig.add_subplot(grid_spec[0, 0])
        ax1.plot(t_arr, ra_residuals, '.', label='RA cos(Dec)', markersize=markersize,
                    color='C1', alpha=0.5)
        ax1.plot(t_arr, dec_residuals, '.', label='Dec', markersize=markersize,
                    color='C0', alpha=0.5)
        ax1.plot(t_arr[rejected_idx], ra_residuals[rejected_idx], 'ro',
                    markersize=2*markersize, markerfacecolor='none', alpha=0.25)
        ax1.plot(t_arr[rejected_idx], dec_residuals[rejected_idx], 'ro',
                    markersize=2*markersize, markerfacecolor='none', alpha=0.25)
        ax1.legend()
        ax1.set_xlabel('Time [UTC]')
        ax1.set_ylabel('Residuals, O-C [arcsec]')
        ax1.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic:
            ax1.set_yscale('log')
        ax2 = grid_spec[0,1].subgridspec(2, 2, width_ratios=(4,1), height_ratios=(1,4),
                                    wspace=0.05, hspace=0.05)
        ax2main = fig.add_subplot(ax2[1,0])
        ax2histx = fig.add_subplot(ax2[0,0], sharex=ax2main)
        ax2histy = fig.add_subplot(ax2[1,1], sharey=ax2main)
        self._scatter_hist(ra_residuals, dec_residuals, ax2main, ax2histx, ax2histy,
                            markersize, show_logarithmic)
        ax2main.plot(ra_residuals[rejected_idx], dec_residuals[rejected_idx], 'ro',
                        markersize=2*markersize, markerfacecolor='none', alpha=0.25)
        ax2main.set_xlabel('RA cos(Dec) Residuals, O-C [arcsec]')
        ax2main.set_ylabel('Dec Residuals, O-C [arcsec]')
        ax2main.grid(True, which='both', axis='both', alpha=0.2, zorder=-100)
        if show_logarithmic:
            ax2main.set_yscale('log')
            ax2main.set_xscale('log')
        ax3 = fig.add_subplot(grid_spec[0, 2])
        ax3.plot(t_arr, doppler_residuals, '.', mfc='C3', mec='C3', label='Doppler',
                    markersize=radar_scale*markersize)
        ax3.plot(t_arr, delay_residuals, '.', mfc='C2', mec='C2', label='Delay',
                    markersize=radar_scale*markersize)
        ax3.legend()
        ax3.set_xlabel('Time [UTC]')
        ax3.set_ylabel(r'Residuals, O-C [$\mu$s, Hz]')
        ax3.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic:
            ax3.set_yscale('log')
        fig.autofmt_xdate()
        plt.tight_layout()
        if savefig:
            figname = f'residuals_iter_{self.iter_number}' if figname is None else figname
            plt.savefig(f'{figname}_residuals.pdf', bbox_inches='tight')
        block = not auto_close
        plt.show(block=block)
        if auto_close:
            plt.close(fig)
        return None

    def _plot_chi(self, t_arr, ra_chi, dec_chi, delay_chi, doppler_chi,
                    ra_chi_squared, dec_chi_squared, delay_chi_squared, doppler_chi_squared,
                    plot_chi_squared, sigma_limit, radar_scale, markersize,
                    show_logarithmic, title, savefig, figname, auto_close):
        # sourcery skip: low-code-quality
        """
        Plot the chi and chi-squared values for the residuals

        Parameters
        ----------
        t_arr : vector
            time array
        ra_chi : vector
            right ascension chi values
        dec_chi : vector
            declination chi values
        delay_chi : vector
            delay chi values
        doppler_chi : vector
            doppler chi values
        ra_chi_squared : vector
            right ascension chi-squared values
        dec_chi_squared : vector
            declination chi-squared values
        delay_chi_squared : vector
            delay chi-squared values
        doppler_chi_squared : vector
            doppler chi-squared values
        plot_chi_squared : bool
            Flag to plot the chi-squared values
        sigma_limit : float
            reference sigma limit for plotting the chi values
        radar_scale : float
            scale factor for the radar scatter points
        markersize : float
            size of the scatter points
        show_logarithmic : bool
            Flag to show the chi values on a logarithmic scale
        title : str
            title for the plot
        savefig : bool
            Flag to save the figure
        figname : str
            name of the figure to save
        auto_close : bool, optional
            Flag to automatically close the figure

        Returns
        -------
        None : NoneType
            None
        """
        rejected_idx = self.rejected_idx
        accepted_idx = self.accepted_idx
        # plot chi values
        factor = 2.5 if plot_chi_squared else 1
        plt.figure(figsize=(factor*7,6), dpi=150)
        if self.iter_number == 0:
            iter_string = f'Iteration {self.iter_number} (prefit)'
        else:
            iter_string = f'Iteration {self.iter_number}'
        iter_string = title if title is not None else iter_string
        msg = iter_string
        if plot_chi_squared:
            msg += ('. Chi Squared: '
                'RA='
                f'{np.nansum(self.all_info["ra_chi_squared"][accepted_idx]):.2f}, '
                'Dec='
                f'{np.nansum(self.all_info["dec_chi_squared"][accepted_idx]):.2f}, '
                'Delay='
                f'{np.nansum(self.all_info["delay_chi_squared"][accepted_idx]):.2f}, '
                'Doppler='
                f'{np.nansum(self.all_info["doppler_chi_squared"][accepted_idx]):.2f}')
        plt.suptitle(msg, y=0.95)
        if plot_chi_squared:
            plt.subplot(1,2,1)
        if not np.all(np.isnan(ra_chi)) and not np.all(np.isnan(dec_chi)):
            plt.plot(t_arr, ra_chi, '.', markersize=markersize,
                        label='RA cos(Dec)', color='C1', alpha=0.5)
            plt.plot(t_arr, dec_chi, '.', markersize=markersize,
                        label='Dec', color='C0', alpha=0.5)
            plt.plot(t_arr[rejected_idx], ra_chi[rejected_idx], 'ro',
                        markersize=2*markersize, markerfacecolor='none', alpha=0.25, label='Rejected')
            plt.plot(t_arr[rejected_idx], dec_chi[rejected_idx], 'ro',
                        markersize=2*markersize, markerfacecolor='none', alpha=0.25)
        if not np.all(np.isnan(doppler_chi)):
            plt.plot(t_arr, doppler_chi, '.', mfc='C3', mec='C3',
                        markersize=radar_scale*markersize, label='Doppler')
        if not np.all(np.isnan(delay_chi)):
            plt.plot(t_arr, delay_chi, '.', mfc='C2', mec='C2',
                        markersize=radar_scale*markersize, label='Delay')
        plt.axhline(-sigma_limit, c='red', linestyle='--', alpha=0.5,
                        label=fr'$\pm{sigma_limit:.0f}\sigma$')
        plt.axhline(sigma_limit, c='red', linestyle='--', alpha=0.5)
        plt.legend(ncol=2)
        plt.xlabel('Time [UTC]')
        plt.ylabel(r'Weighted Residuals, (O-C) $[\sigma]$')
        plt.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic:
            plt.yscale('log')
        else:
            lim = np.max(np.abs(plt.ylim()))
            plt.ylim(-lim, lim)
        plt.gcf().autofmt_xdate()
        if plot_chi_squared:
            plt.subplot(1,2,2)
            plt.plot(t_arr, ra_chi_squared, '.', markersize=markersize,
                        label='RA cos(Dec)', color='C1', alpha=0.5)
            plt.plot(t_arr, dec_chi_squared, '.', markersize=markersize,
                        label='Dec', color='C0', alpha=0.5)
            plt.plot(t_arr[rejected_idx], ra_chi_squared[rejected_idx], 'ro',
                        markersize=2*markersize, markerfacecolor='none', alpha=0.25)
            plt.plot(t_arr[rejected_idx], dec_chi_squared[rejected_idx], 'ro',
                        markersize=2*markersize, markerfacecolor='none', alpha=0.25)
            plt.plot(t_arr, doppler_chi_squared, '.', mfc='C3', mec='C3',
                        markersize=radar_scale*markersize, label='Doppler')
            plt.plot(t_arr, delay_chi_squared, '.', mfc='C2', mec='C2',
                        markersize=radar_scale*markersize, label='Delay')
            plt.legend(ncol=2)
            plt.xlabel('Time [UTC]')
            plt.ylabel(r'$\chi^2$')
            plt.grid(True, which='both', axis='both', alpha=0.2)
            # if show_logarithmic: plt.yscale('log')
            plt.yscale('log')
            plt.tight_layout()
            plt.gcf().autofmt_xdate()
        if savefig:
            figname = f'chi_iter_{self.iter_number}' if figname is None else figname
            plt.savefig(f'{figname}_chi.pdf', bbox_inches='tight')
        block = not auto_close
        plt.show(block=block)
        if auto_close:
            plt.close()
        return None

    def plot_iteration_summary(self, show_logarithmic=False, title=None, plot_chi_squared=False,
                                savefig=False, figname=None, auto_close=False):
        """
        Plot the summary of the iteration, including residuals, chi values,
        and chi squared values.

        Parameters
        ----------
        show_logarithmic : bool, optional
            Flag to show the plotted values on a logarithmic scale, by default False
        title : str, optional
            Title of the plot, by default None
        plot_chi_squared : bool, optional
            Flag to plot the chi squared values, by default False
        savefig : bool, optional
            Flag to save the figure, by default False
        figname : str, optional
            Name of the figure, by default None
        auto_close : bool, optional
            Flag to automatically close the figure, by default False

        Returns
        -------
        None : NoneType
            None
        """
        markersize = 3
        sigma_limit = 3
        radar_scale = 3
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.labelsize'] = 12
        t_arr = self.obs['obsTimeMJD'].values
        t_arr = Time(t_arr, format='mjd', scale='utc').utc.datetime
        if show_logarithmic:
            ra_residuals = np.abs(self.all_info['ra_res'])
            dec_residuals = np.abs(self.all_info['dec_res'])
            ra_chi = np.abs(self.all_info['ra_chi'])
            dec_chi = np.abs(self.all_info['dec_chi'])
        else:
            ra_residuals = self.all_info['ra_res']
            dec_residuals = self.all_info['dec_res']
            ra_chi = self.all_info['ra_chi']
            dec_chi = self.all_info['dec_chi']
        ra_chi_squared = self.all_info['ra_chi_squared']
        dec_chi_squared = self.all_info['dec_chi_squared']
        if show_logarithmic:
            delay_residuals = np.abs(self.all_info['delay_res'])
            doppler_residuals = np.abs(self.all_info['doppler_res'])
            delay_chi = np.abs(self.all_info['delay_chi'])
            doppler_chi = np.abs(self.all_info['doppler_chi'])
        else:
            delay_residuals = self.all_info['delay_res']
            doppler_residuals = self.all_info['doppler_res']
            delay_chi = self.all_info['delay_chi']
            doppler_chi = self.all_info['doppler_chi']
        delay_chi_squared = self.all_info['delay_chi_squared']
        doppler_chi_squared = self.all_info['doppler_chi_squared']
        self._plot_residuals(t_arr, ra_residuals, dec_residuals, delay_residuals,
                                doppler_residuals, radar_scale, markersize,
                                show_logarithmic, title, savefig, figname, auto_close)
        self._plot_chi(t_arr, ra_chi, dec_chi, delay_chi, doppler_chi,
                        ra_chi_squared, dec_chi_squared, delay_chi_squared, doppler_chi_squared,
                        plot_chi_squared, sigma_limit, radar_scale, markersize,
                        show_logarithmic, title, savefig, figname, auto_close)
        return None

class FitSimulation:
    """
    Class to perform an orbit fit simulation.
    """
    def __init__(self, x_init, obs_df, cov_init=None, n_iter_max=10,
                    de_kernel=441, nongrav_info=None, events=None,
                    simulated_obs=None):
        """
        Constructor of the FitSimulation class.

        Parameters
        ----------
        x_init : dict
            Initial guess of the state vector.
        obs_df : pandas DataFrame
            Observation data for the orbit fit
        cov_init : array, optional
            Initial covariance matrix, by default None
        n_iter_max : int, optional
            Number of maximum iterations to correct the orbit estimate, by default 10
        de_kernel : int, optional
            SPICE kernel version, by default 441
        nongrav_info : dict, optional
            Dictionary containing the non-gravitational parameters, by default None
        events : list, optional
            List of lists containing impulsive maneuvers, by default None
        simulated_obs : dict, optional
            Dictionary containing simulated observation data, by default None

        Returns
        -------
        None : NoneType
            None
        """
        perm_id = obs_df.iloc[-1]['permID']
        prov_id = obs_df.iloc[-1]['provID']
        body_id = perm_id if isinstance(perm_id, str) else prov_id
        self.name = body_id
        self.t_sol = None
        self.x_init = None
        self.x_nom = None
        self.covariance_init = cov_init
        self.covariance = None
        self.fit_cartesian = False
        self.fit_cometary = False
        self.n_fit = None
        self._check_initial_solution(x_init, cov_init)
        self.constraint_dir = None
        self.obs = None
        self.observer_info = None
        self.optical_idx = None
        self.delay_idx = None
        self.doppler_idx = None
        self.radar_idx = None
        self.n_obs = None
        self.past_obs_idx = None
        self.past_obs_exist = None
        self.future_obs_idx = None
        self.future_obs_exist = None
        self.simulated_obs = simulated_obs
        self._parse_observation_arrays(obs_df)
        self.obs_cov = None
        self.obs_weight = None
        self._compute_obs_weights()
        self.n_iter = 0
        self.n_iter_max = n_iter_max
        self.iters = [[]]
        self.de_kernel = de_kernel
        self.de_kernel_path = default_kernel_path
        self.analytic_partials = False
        self.prop_sims = [None, None]
        self.fixed_propsim_params = {'a1': 0.0, 'a2': 0.0, 'a3': 0.0,
                                        'alpha': 1.0, 'k': 0.0, 'm': 2.0, 'n': 0.0,
                                        'r0_au': 1.0, 'radius': 0.0, 'mass': 0.0}
        if nongrav_info is not None:
            for key in nongrav_info:
                self.fixed_propsim_params[key] = nongrav_info[key]
        self.fixed_propsim_params['events'] = events if events is not None else []
        self.prior_est = None
        self.prior_sig = None
        self._priors_given = False
        self._xbar0 = None
        self._info0 = None
        self._prior_constant = None
        self.reject_outliers = True
        self.reject_criteria = [3.0, 2.8]
        # number of rejected obs is the count of ['D', 'd'] in selAst
        sel_ast = self.obs['selAst'].values
        num_auto_rejected = np.sum(sel_ast == 'D')
        num_force_rejected = np.sum(sel_ast == 'd')
        self.num_rejected = num_auto_rejected + num_force_rejected
        self.info_mats = []
        self.converged = False
        return None

    def _check_initial_solution(self, x_init, cov_init):
        """
        Check the initial solution provided by the user and make 
        sure it is valid.

        Parameters
        ----------
        x_init : dict
            Initial guess of the state vector.
        cov_init : array
            Initial covariance matrix.

        Returns
        -------
        None : NoneType
            None

        Raises
        ------
        ValueError
            If the initial solution does not contain the time.
        ValueError
            If the initial solution does not contain a full cartesian
            or cometary state.
        ValueError
            If the size of the covariance matrix does not match the
            number of fitted parameters.
        """
        if 't' not in x_init:
            raise ValueError("Must provide a time for the initial solution.")
        if all(key in x_init for key in ("x", "y", "z", "vx", "vy", "vz")):
            self.fit_cartesian = True
            self.fit_cometary = False
        elif all(key in x_init for key in ("e", "q", "tp", "om", "w", "i")):
            self.fit_cartesian = False
            self.fit_cometary = True
        else:
            msg = ("Must provide at least a full cartesian",
                    "or cometary state for the initial solution.")
            raise ValueError(msg)
        self.t_sol = x_init['t']
        self.x_init = {key: x_init[key] for key in x_init if key != 't'}
        self.x_nom = self.x_init.copy()
        self.n_fit = len(self.x_nom)
        if cov_init.shape != (self.n_fit, self.n_fit):
            msg = ("Covariance matrix must be the same size "
                    "as the number of fitted parameters.")
            raise ValueError(msg)
        return None

    def _check_priors(self):
        """
        Check the prior estimates and sigmas provided by the user
        and make sure they are valid.

        Returns
        -------
        None : NoneType
            None

        Raises
        ------
        ValueError
            If the prior estimates and sigmas do not have the same length.
        """
        # check that the keys in self.x_init are the same as in self.prior_est and self.prior_sig
        if self.prior_est is not None and self.prior_sig is not None:
            self._priors_given = True
            if len(self.prior_est) != len(self.prior_sig):
                raise ValueError("Prior estimates and sigmas must have the same length.")
            for key in self.prior_est.keys():
                if key not in self.x_init:
                    raise ValueError(f"Key {key} not found in initial state vector.")
            self._info0 = np.zeros_like(self.covariance_init)
            self._xbar0 = np.zeros(self.n_fit)
            for key in self.prior_est.keys():
                idx = list(self.x_init.keys()).index(key)
                self._info0[idx, idx] = 1/self.prior_sig[key]**2
                self._xbar0[idx] = self.prior_est[key] - self.x_init[key]
            self._prior_constant = list(self.x_nom.values())+self._xbar0
        return None

    def _flatten_and_clean(self, arr):
        """
        Flatten an array and remove any NaN values.

        Parameters
        ----------
        arr : array
            Array to flatten and clean.

        Returns
        -------
        arr : array
            Flattened and cleaned array.
        """
        arr = arr.flatten()
        arr = arr[~np.isnan(arr)]
        return arr

    def _add_simulated_obs(self):
        """
        Add the simulated observation data to the observation data.

        Returns
        -------
        None : NoneType
            None
        """
        times = self.simulated_obs['times']
        modes = self.simulated_obs['modes']
        weights = self.simulated_obs['weights']
        if not len(times) == len(modes) == len(weights):
            raise ValueError("Simulated observation data must have the same length.")
        perm_id = self.obs['permID'][0]
        prov_id = self.obs['provID'][0]
        new_sim_obs = []
        for i, time in enumerate(times):
            this_obs = {}
            mode = modes[i]
            if mode not in {'SIM_CCD', 'SIM_OCC', 'SIM_RAD_DEL', 'SIM_RAD_DOP'}:
                raise ValueError(f"Unknown simulated observation mode {mode}.")
            weight = weights[i]
            if isinstance(perm_id, str):
                this_obs['permID'] = f'SIM_{perm_id}'
            if isinstance(prov_id, str):
                this_obs['provID'] = f'SIM_{prov_id}'
            this_obs['obsTime'] = f'{time.utc.isot}Z'
            this_obs['obsTimeMJD'] = time.utc.mjd
            this_obs['obsTimeMJDTDB'] = time.tdb.mjd
            this_obs['mode'] = mode
            if mode in {'SIM_CCD', 'SIM_OCC'}:
                this_obs['stn'] = 'SIM'
                this_obs['ra'] = np.nan
                this_obs['rmsRA'] = weight[0]
                this_obs['sigRA'] = weight[0]
                this_obs['dec'] = np.nan
                this_obs['rmsDec'] = weight[1]
                this_obs['sigDec'] = weight[1]
                this_obs['rmsCorr'] = weight[2]
                this_obs['sigCorr'] = weight[2]
                this_obs['cosDec'] = 1.0
            elif mode in {'SIM_RAD_DEL', 'SIM_RAD_DOP'}:
                this_obs['trx'] = 'SIM'
                this_obs['rcv'] = 'SIM'
                this_obs['com'] = 1
                if mode == 'SIM_RAD_DEL':
                    this_obs['delay'] = np.nan
                    this_obs['rmsDelay'] = weight[0]
                    this_obs['sigDelay'] = weight[0]
                    this_obs['doppler'] = np.nan
                    this_obs['rmsDoppler'] = np.nan
                    this_obs['sigDoppler'] = np.nan
                else:
                    this_obs['delay'] = np.nan
                    this_obs['rmsDelay'] = np.nan
                    this_obs['sigDelay'] = np.nan
                    this_obs['doppler'] = np.nan
                    this_obs['rmsDoppler'] = weight[0]
                    this_obs['sigDoppler'] = weight[0]
                    this_obs['frq'] = 8560.0
            this_obs['ctr'] = 399
            this_obs['sys'] = 'ITRF'
            this_obs['pos1'] = 0.0
            this_obs['pos2'] = 0.0
            this_obs['pos3'] = 0.0
            this_obs['selAst'] = 'a'
            new_sim_obs.append(this_obs)
        sim_obs_df = pd.DataFrame(new_sim_obs)
        self.obs = pd.concat([self.obs, sim_obs_df], ignore_index=True)
        return None

    def _parse_observation_arrays(self, obs_df):
        """
        Parse the observation data for the orbit fit.

        Returns
        -------
        None : NoneType
            None
        """
        self.obs = obs_df
        sel_ast_map = {'': 'A', ' ': 'A', np.nan: 'A', 'nan': 'A'}
        self.obs['selAst'] = self.obs['selAst'].replace(sel_ast_map)
        if self.simulated_obs is not None:
            self._add_simulated_obs()
        self.obs.sort_values(by='obsTimeMJD', inplace=True, ignore_index=True)
        self.observer_info = get_observer_info(self.obs)
        self.optical_idx = []
        self.delay_idx = []
        self.doppler_idx = []
        self.observer_info_lengths = [len(info) for info in self.observer_info]
        for i, length in enumerate(self.observer_info_lengths):
            if length in {4,7}: # optical
                self.optical_idx.append(i)
            elif length == 9: # delay
                self.delay_idx.append(i)
            elif length == 10: # doppler
                self.doppler_idx.append(i)
            else:
                raise ValueError('Unknown observer type')
        self.radar_idx = self.delay_idx + self.doppler_idx
        self.n_obs = 2*len(self.optical_idx)+len(self.radar_idx)

        self.past_obs_idx = self.obs.query('obsTimeMJDTDB < @self.t_sol').index
        self.past_obs_exist = len(self.past_obs_idx) > 0
        self.future_obs_idx = self.obs.query('obsTimeMJDTDB >= @self.t_sol').index
        self.future_obs_exist = len(self.future_obs_idx) > 0
        return None

    def _compute_obs_weights(self):
        """
        Assembles the weight matrix for the orbit fit based
        on observation uncertainties and correlations.

        Returns
        -------
        None : NoneType
            None
        """
        self.obs_cov = []
        self.obs_weight = []
        fields = ['mode', 'sigRA', 'sigDec', 'sigCorr', 'sigDelay', 'sigDoppler']
        for info in zip(*[self.obs[field] for field in fields]):
            mode, sig_ra, sig_dec, sig_corr, sig_delay, sig_doppler = info
            if mode in {'RAD', 'SIM_RAD_DEL', 'SIM_RAD_DOP'}:
                sig = sig_doppler if np.isfinite(sig_doppler) else sig_delay
                if sig == 0.0 or np.isnan(sig):
                    raise ValueError("Radar uncertainty is 0 or NaN.")
                self.obs_cov.append([[sig**2]])
                self.obs_weight.append([[1.0/sig**2]])
            else:
                off_diag = 0.0 if np.isnan(sig_corr) else sig_corr*sig_ra*sig_dec
                cov = np.array([[sig_ra**2, off_diag],
                                    [off_diag, sig_dec**2]])
                det = cov[0, 0]*cov[1, 1] - cov[0, 1]*cov[1, 0]
                if det == 0.0:
                    print(cov)
                    raise ValueError("Optical covariance matrix is singular.")
                inv = np.array([[cov[1, 1], -cov[0, 1]],
                                [-cov[1, 0], cov[0, 0]]])/det
                self.obs_cov.append(cov)
                self.obs_weight.append(inv)
        return None

    def _get_prop_sim_past(self, name, t_eval_utc, eval_apparent_state,
                                converged_light_time, observer_info):
        """
        Get a propSim object for the past observations.

        Parameters
        ----------
        name : str
            Name of the propSim object.
        t_eval_utc : bool
            Flag for whether the observation times are in UTC.
        eval_apparent_state : bool
            Flag for whether to evaluate the apparent state.
        converged_light_time : bool
            Flag for whether to use converged Newtonian light time correction.
        observer_info : tuple
            Observer information for each past observation.

        Returns
        -------
        prop_sim_past : libgrss.PropSimulation object
            propSim object for the past observations.
        """
        # pylint: disable=no-member
        t_eval_past = self.obs.obsTimeMJD.values[self.past_obs_idx]
        tf_past = np.min(self.obs.obsTimeMJDTDB.values[self.past_obs_idx])
        prop_sim_past = libgrss.PropSimulation(name, self.t_sol,
                                                self.de_kernel, self.de_kernel_path)
        prop_sim_past.tEvalMargin = 1.0
        # flip t_eval_past and observer_info to go in reverse time order
        t_eval_past = t_eval_past[::-1]
        observer_info = observer_info[::-1]
        prop_sim_past.set_integration_parameters(tf_past, t_eval_past, t_eval_utc,
                                                    eval_apparent_state, converged_light_time,
                                                    observer_info)
        prop_sim_past.evalMeasurements = True
        # set prop_sim_past.obsType = 3 where gaia_flag is True
        gaia_flag = (self.obs.loc[self.past_obs_idx]['stn'] == '258').values[::-1]
        obs_types = prop_sim_past.obsType
        for idx in np.where(gaia_flag)[0]:
            obs_types[idx] = 3
        prop_sim_past.obsType = obs_types
        return prop_sim_past

    def _get_prop_sim_future(self, name, t_eval_utc, eval_apparent_state,
                                converged_light_time, observer_info):
        """
        Get a propSim object for the future observations.

        Parameters
        ----------
        name : str
            Name of the propSim object.
        t_eval_utc : bool
            Flag for whether the observation times are in UTC.
        eval_apparent_state : bool
            Flag for whether to evaluate the apparent state.
        converged_light_time : bool
            Flag for whether to use converged Newtonian light time correction.
        observer_info : tuple
            Observer information for each future observation.

        Returns
        -------
        prop_sim_future : libgrss.PropSimulation object
            propSim object for the future observations.
        """
        # pylint: disable=no-member
        t_eval_future = self.obs.obsTimeMJD.values[self.future_obs_idx]
        tf_future = np.max(self.obs.obsTimeMJDTDB.values[self.future_obs_idx])
        prop_sim_future = libgrss.PropSimulation(name, self.t_sol,
                                                    self.de_kernel, self.de_kernel_path)
        prop_sim_future.tEvalMargin = 1.0
        prop_sim_future.set_integration_parameters(tf_future, t_eval_future, t_eval_utc,
                                                    eval_apparent_state, converged_light_time,
                                                    observer_info)
        prop_sim_future.evalMeasurements = True
        # set prop_sim_future.obsType = 3 where gaia_flag is True
        gaia_flag = (self.obs.loc[self.future_obs_idx]['stn'] == '258').values
        obs_types = prop_sim_future.obsType
        for idx in np.where(gaia_flag)[0]:
            obs_types[idx] = 3
        prop_sim_future.obsType = obs_types
        return prop_sim_future

    def _get_prop_sims(self):
        """
        Get propSim objects for the past and future observations.

        Parameters
        ----------
        name : str
            Name of the propSim objects.

        Returns
        -------
        prop_sim_past : libgrss.PropSimulation object
            propSim object for the past observations.
        prop_sim_future : libgrss.PropSimulation object
            propSim object for the future observations.
        """
        t_eval_utc = True
        eval_apparent_state = True
        converged_light_time = True
        observer_info = np.array(self.observer_info, dtype=tuple)
        observer_info_past = observer_info[self.past_obs_idx]
        observer_info_future = observer_info[self.future_obs_idx]
        prop_sim_past = None
        prop_sim_future = None
        if self.past_obs_exist:
            prop_sim_past = self._get_prop_sim_past(f"{self.name}_past", t_eval_utc,
                                                    eval_apparent_state, converged_light_time,
                                                    observer_info_past)
        if self.future_obs_exist:
            prop_sim_future = self._get_prop_sim_future(f"{self.name}_future", t_eval_utc,
                                                        eval_apparent_state, converged_light_time,
                                                        observer_info_future)
        return prop_sim_past, prop_sim_future

    def _x_dict_to_state(self, x_dict):
        """
        Convert a dictionary of nominal state to a state vector.

        Parameters
        ----------
        x_dict : dict
            Dictionary of the nominal state.

        Returns
        -------
        state : list
            State vector.

        Raises
        ------
        ValueError
            If the state type is not well-defined.
        """
        # sourcery skip: extract-method
        if self.fit_cartesian:
            pos_x = x_dict['x']
            pos_y = x_dict['y']
            pos_z = x_dict['z']
            vel_x = x_dict['vx']
            vel_y = x_dict['vy']
            vel_z = x_dict['vz']
            state = [pos_x, pos_y, pos_z, vel_x, vel_y, vel_z]
        elif self.fit_cometary:
            ecc = x_dict['e']
            peri_dist = x_dict['q']
            time_peri = x_dict['tp']
            omega = x_dict['om']
            arg_peri = x_dict['w']
            inc = x_dict['i']
            cometary_elements = [ecc, peri_dist, time_peri,
                                    omega, arg_peri, inc]
            state = cometary_elements
        else:
            raise ValueError("fit_cartesian or fit_cometary must be True")
        return state

    def _x_dict_to_nongrav_params(self, x_dict):
        """
        Convert a dictionary of nominal state to non-gravitational parameters.

        Parameters
        ----------
        x_dict : dict
            Dictionary of the nominal state.

        Returns
        -------
        nongrav_params : libgrss.NongravParameters object
            Non-gravitational parameters for the fitted body.
        """
        # pylint: disable=no-member
        nongrav_params = libgrss.NongravParameters()
        a1_est = 'a1' in x_dict.keys()
        a2_est = 'a2' in x_dict.keys()
        a3_est = 'a3' in x_dict.keys()
        nongrav_params.a1 = x_dict['a1'] if a1_est else self.fixed_propsim_params['a1']
        nongrav_params.a2 = x_dict['a2'] if a2_est else self.fixed_propsim_params['a2']
        nongrav_params.a3 = x_dict['a3'] if a3_est else self.fixed_propsim_params['a3']
        nongrav_params.a1Est = a1_est
        nongrav_params.a2Est = a2_est
        nongrav_params.a3Est = a3_est
        nongrav_params.alpha = self.fixed_propsim_params['alpha']
        nongrav_params.k = self.fixed_propsim_params['k']
        nongrav_params.m = self.fixed_propsim_params['m']
        nongrav_params.n = self.fixed_propsim_params['n']
        nongrav_params.r0_au = self.fixed_propsim_params['r0_au']
        return nongrav_params

    def _x_dict_to_events(self, x_dict):
        """
        Convert a dictionary of nominal state to events.

        Parameters
        ----------
        x_dict : dict
            Dictionary of the nominal state.

        Returns
        -------
        events : list
            List of events for the fitted body.
        """
        events = []
        for i in range(len(self.fixed_propsim_params['events'])):
            fixed_event = tuple(self.fixed_propsim_params['events'][i])
            fixed_event_impulsive = len(fixed_event) == 5
            fixed_event_continuous = len(fixed_event) == 6
            if fixed_event_continuous and fixed_event[4] == 0.0:
                raise ValueError("Continuous event time constant cannot be 0.")
            assert fixed_event_impulsive or fixed_event_continuous
            event = [None]*5
            event[0] = fixed_event[0] # time of delta-v
            if fixed_event_impulsive:
                event[1] = x_dict[f"dvx{i}"] if f"dvx{i}" in x_dict.keys() else fixed_event[1]
                event[2] = x_dict[f"dvy{i}"] if f"dvy{i}" in x_dict.keys() else fixed_event[2]
                event[3] = x_dict[f"dvz{i}"] if f"dvz{i}" in x_dict.keys() else fixed_event[3]
                event[4] = x_dict[f"mult{i}"] if f"mult{i}" in x_dict.keys() else fixed_event[4]
            elif fixed_event_continuous:
                event[1] = x_dict[f"ax{i}"] if f"ax{i}" in x_dict.keys() else fixed_event[1]
                event[2] = x_dict[f"ay{i}"] if f"ay{i}" in x_dict.keys() else fixed_event[2]
                event[3] = x_dict[f"az{i}"] if f"az{i}" in x_dict.keys() else fixed_event[3]
                event[4] = x_dict[f"dt{i}"] if f"dt{i}" in x_dict.keys() else fixed_event[4]
                event.append(fixed_event[5])
            events.append(tuple(event))
        return events

    def _check_and_add_events(self, prop_sim_past, prop_sim_future, integ_body, events):
        """
        Check if events are in the past or future and add them to the appropriate prop_sim.

        Parameters
        ----------
        prop_sim_past : libgrss.PropSimulation object
            PropSimulation object for the past.
        prop_sim_future : libgrss.PropSimulation object
            PropSimulation object for the future.
        integ_body : libgrss.IntegBody object
            IntegBody object for the body being fitted.
        events : list
            List of events for the fitted body.

        Returns
        -------
        prop_sim_past : libgrss.PropSimulation object
            PropSimulation object for the past.
        prop_sim_future : libgrss.PropSimulation object
            PropSimulation object for the future.
        """
        est_keys = self.x_nom.keys()
        for i, event_list in enumerate(events):
            t_event = event_list[0]
            dvx_or_ax = event_list[1]
            dvy_or_ay = event_list[2]
            dvz_or_az = event_list[3]
            multiplier_or_dt = event_list[4]
            continuous_flag = event_list[5] if len(event_list) == 6 else False
            event = libgrss.Event()
            event.t = t_event
            event.bodyName = integ_body.name
            event.isContinuous = continuous_flag
            if continuous_flag:
                event.expAccel0 = [dvx_or_ax, dvy_or_ay, dvz_or_az]
                event.tau = multiplier_or_dt
                event.expAccel0Est = (
                    f"ax{i}" in est_keys and
                    f"ay{i}" in est_keys and
                    f"az{i}" in est_keys
                )
                event.tauEst = f"dt{i}" in est_keys
                event.eventEst = event.expAccel0Est or event.tauEst
            else:
                event.deltaV = [dvx_or_ax, dvy_or_ay, dvz_or_az]
                event.multiplier = multiplier_or_dt
                event.deltaVEst = (
                    f"dvx{i}" in est_keys and
                    f"dvy{i}" in est_keys and 
                    f"dvz{i}" in est_keys
                )
                event.multiplierEst = f"mult{i}" in est_keys
                event.eventEst = event.deltaVEst or event.multiplierEst
            if t_event < self.t_sol:
                prop_sim_past.add_event(event)
            else:
                prop_sim_future.add_event(event)
        return prop_sim_past, prop_sim_future

    def _get_perturbed_state(self, key):
        """
        Get the perturbed state for a given nominal state parameter.

        Parameters
        ----------
        key : str
            Key for the nominal state parameter.

        Returns
        -------
        state_plus : list
            Perturbed state in the positive direction.
        ng_params_plus : libgrss.NongravParameters object
            Perturbed non-gravitational parameters in the positive direction.
        events_plus : list
            Perturbed events in the positive direction.
        state_minus : list
            Perturbed state in the negative direction.
        ng_params_minus : libgrss.NongravParameters object
            Perturbed non-gravitational parameters in the negative direction.
        events_minus : list
            Perturbed events in the negative direction.
        fd_delta : float
            Finite difference perturbation.
        """
        fd_pert = np.inf
        if self.fit_cartesian:
            if key in {'x', 'y', 'z'}:
                fd_pert = 1e-8
            elif key in {'vx', 'vy', 'vz'}:
                fd_pert = 1e-10
        elif self.fit_cometary:
            fd_pert = 1e-9
            if key == 'e':
                fd_pert = 7.5e-9
            elif key == 'q':
                fd_pert = 1e-9
            elif key == 'tp':
                fd_pert = 1e-7
        if key in {'a1', 'a2', 'a3'}:
            fd_pert = 1e-12
        if key[:4] == 'mult':
            fd_pert = 1e0
        if key[:3] in {'dvx', 'dvy', 'dvz'}:
            fd_pert = 1e-11
        if key[:2] in {'ax', 'ay', 'az'}:
            fd_pert = 1e-11
        if key[:2] == 'dt':
            # 10 percent of the time constant
            fd_pert = 0.1*self.x_nom[key]
        if fd_pert == np.inf:
            raise ValueError("Finite difference perturbation not defined"
                                " for the given state parameter.")
        x_plus = self.x_nom.copy()
        x_minus = self.x_nom.copy()
        # fd_pert = finite difference perturbation to nominal state for calculating derivatives
        fd_delta = fd_pert
        x_plus[key] = self.x_nom[key]+fd_delta
        state_plus = self._x_dict_to_state(x_plus)
        ng_params_plus = self._x_dict_to_nongrav_params(x_plus)
        events_plus = self._x_dict_to_events(x_plus)
        x_minus[key] = self.x_nom[key]-fd_delta
        state_minus = self._x_dict_to_state(x_minus)
        ng_params_minus = self._x_dict_to_nongrav_params(x_minus)
        events_minus = self._x_dict_to_events(x_minus)
        return (state_plus, ng_params_plus, events_plus, state_minus,
                    ng_params_minus, events_minus, fd_delta)

    def _get_perturbation_info(self):
        """
        Get the perturbation information for all nominal state parameters.

        Returns
        -------
        perturbation_info : list
            List of tuples containing the perturbation information
            for each nominal state parameter.
        """
        perturbation_info = []
        for key in self.x_nom:
            pert_result = self._get_perturbed_state(key)
            perturbation_info.append(tuple(pert_result))
        return perturbation_info

    def _assemble_and_propagate_bodies(self, perturbation_info):
        """
        Assemble and propagate the body for the orbit fit.

        Parameters
        ----------
        perturbation_info : list
            List of tuples containing the perturbation information
            for each nominal state parameter.

        Returns
        -------
        prop_sim_past : libgrss.PropSimulation object
            propagated PropSimulation object for the past.
        prop_sim_future : libgrss.PropSimulation object
            propagated PropSimulation object for the future.
        """
        # sourcery skip: low-code-quality
        # pylint: disable=no-member
        # get propagated states
        prop_sim_past, prop_sim_future = self._get_prop_sims()
        # create nominal integ_body object
        state_nom = self._x_dict_to_state(self.x_nom)
        ng_params_nom = self._x_dict_to_nongrav_params(self.x_nom)
        events_nom = self._x_dict_to_events(self.x_nom)
        if self.fit_cartesian:
            integ_body_nom = libgrss.IntegBody("integ_body_nom", self.t_sol,
                                            self.fixed_propsim_params['mass'],
                                            self.fixed_propsim_params['radius'],
                                            state_nom[:3], state_nom[3:6], ng_params_nom)
        elif self.fit_cometary:
            integ_body_nom = libgrss.IntegBody("integ_body_nom",
                                            self.t_sol, self.fixed_propsim_params['mass'],
                                            self.fixed_propsim_params['radius'], state_nom,
                                            ng_params_nom)
        integ_body_nom.caTol = 0.0 # turn off close approach detection
        # add the nominal integ_body for the residuals
        if self.analytic_partials:
            integ_body_nom.prepare_stm()
        if self.past_obs_exist:
            prop_sim_past.add_integ_body(integ_body_nom)
        if self.future_obs_exist:
            prop_sim_future.add_integ_body(integ_body_nom)
        prop_sim_past, prop_sim_future = self._check_and_add_events(prop_sim_past, prop_sim_future,
                                                                    integ_body_nom, events_nom)
        # add the perturbed IntegBodies for numerical derivatives
        if not self.analytic_partials and perturbation_info is not None:
            for i in range(self.n_fit):
                key = list(self.x_nom.keys())[i]
                (state_plus, ng_params_plus, events_plus, state_minus,
                        ng_params_minus, events_minus, _) = perturbation_info[i]
                if self.fit_cartesian:
                    integ_body_plus = libgrss.IntegBody(f"integBody_pert_{key}_plus", self.t_sol,
                                                        self.fixed_propsim_params['mass'],
                                                        self.fixed_propsim_params['radius'],
                                                        state_plus[:3], state_plus[3:6],
                                                        ng_params_plus)
                    integ_body_minus = libgrss.IntegBody(f"integBody_pert_{key}_minus", self.t_sol,
                                                        self.fixed_propsim_params['mass'],
                                                        self.fixed_propsim_params['radius'],
                                                        state_minus[:3], state_minus[3:6],
                                                        ng_params_minus)
                elif self.fit_cometary:
                    integ_body_plus = libgrss.IntegBody(f"integBody_pert_{key}_plus", self.t_sol,
                                                        self.fixed_propsim_params['mass'],
                                                        self.fixed_propsim_params['radius'],
                                                        state_plus, ng_params_plus)
                    integ_body_minus = libgrss.IntegBody(f"integBody_pert_{key}_minus", self.t_sol,
                                                        self.fixed_propsim_params['mass'],
                                                        self.fixed_propsim_params['radius'],
                                                        state_minus, ng_params_minus)
                integ_body_plus.caTol = 0.0
                integ_body_minus.caTol = 0.0
                if self.past_obs_exist:
                    prop_sim_past.add_integ_body(integ_body_plus)
                    prop_sim_past.add_integ_body(integ_body_minus)
                if self.future_obs_exist:
                    prop_sim_future.add_integ_body(integ_body_plus)
                    prop_sim_future.add_integ_body(integ_body_minus)
                prop_sim_past, prop_sim_future = self._check_and_add_events(prop_sim_past,
                                                                            prop_sim_future,
                                                                            integ_body_plus,
                                                                            events_plus)
                prop_sim_past, prop_sim_future = self._check_and_add_events(prop_sim_past,
                                                                            prop_sim_future,
                                                                            integ_body_minus,
                                                                            events_minus)
        if self.past_obs_exist:
            prop_sim_past.integrate()
        if self.future_obs_exist:
            prop_sim_future.integrate()
        return prop_sim_past, prop_sim_future

    def _get_computed_obs(self, prop_sim_past, prop_sim_future, integ_body_idx):
        """
        Computes the optical and radar observations from the propagated states.

        Parameters
        ----------
        prop_sim_past : libgrss.PropSimulation object
            The propagated PropSimulation object for the past.
        prop_sim_future : libgrss.PropSimulation object
            The propagated PropSimulation object for the future.
        integ_body_idx : int
            The index of the integ_body to use for calculating the observations.

        Returns
        -------
        computed_obs : array
            The computed observations from the propagated states.

        Raises
        ------
        ValueError
            If the observer information is not well-defined.
        """
        if self.past_obs_exist and self.future_obs_exist:
            optical_obs = prop_sim_past.opticalObs + prop_sim_future.opticalObs
            optical_obs_corr = prop_sim_past.opticalObsCorr + prop_sim_future.opticalObsCorr
            radar_obs = prop_sim_past.radarObs + prop_sim_future.radarObs
        elif self.past_obs_exist:
            optical_obs = prop_sim_past.opticalObs
            optical_obs_corr = prop_sim_past.opticalObsCorr
            radar_obs = prop_sim_past.radarObs
        elif self.future_obs_exist:
            optical_obs = prop_sim_future.opticalObs
            optical_obs_corr = prop_sim_future.opticalObsCorr
            radar_obs = prop_sim_future.radarObs
        computed_obs = np.nan*np.ones((len(self.obs), 2))
        cos_dec = self.obs.cosDec.values
        for i, obs_info_len in enumerate(self.observer_info_lengths):
            if obs_info_len in {4, 7}:
                computed_obs[i, :] = optical_obs[i][2*integ_body_idx:2*integ_body_idx+2]
                computed_obs[i, 0] *= cos_dec[i]
                correction = optical_obs_corr[i][2*integ_body_idx:2*integ_body_idx+2]
                computed_obs[i, :] += correction
            elif obs_info_len == 9: # delay measurement
                computed_obs[i, 0] = radar_obs[i][integ_body_idx]
            elif obs_info_len == 10: # dopper measurement
                computed_obs[i, 1] = radar_obs[i][integ_body_idx]
            else:
                raise ValueError("Observer info length not recognized.")
        nom_body = integ_body_idx == 0
        if nom_body:
            self._inflate_uncertainties(prop_sim_past, prop_sim_future)
        return computed_obs

    def _inflate_uncertainties(self, prop_sim_past, prop_sim_future):
        """
        Apply time uncertainties to the optical observations weights.

        Parameters
        ----------
        computed_obs_dot : array
            Computed optical observations dot.

        Returns
        -------
        None : NoneType
            None
        """
        stations = self.obs.stn.values
        sig_times = self.obs.sigTime.values
        sig_ra_vals = self.obs.sigRA.values
        sig_dec_vals = self.obs.sigDec.values
        sig_corr_vals = self.obs.sigCorr.values
        sys_vals = self.obs.sys.values
        p11_vals = self.obs.posCov11.values
        p12_vals = self.obs.posCov12.values
        p13_vals = self.obs.posCov13.values
        p22_vals = self.obs.posCov22.values
        p23_vals = self.obs.posCov23.values
        p33_vals = self.obs.posCov33.values
        if self.past_obs_exist and self.future_obs_exist:
            optical_obs_dot = prop_sim_past.opticalObsDot + prop_sim_future.opticalObsDot
            optical_obs_partials = prop_sim_past.opticalPartials + prop_sim_future.opticalPartials
            state_eval = prop_sim_past.xIntegEval + prop_sim_future.xIntegEval
        elif self.past_obs_exist:
            optical_obs_dot = prop_sim_past.opticalObsDot
            optical_obs_partials = prop_sim_past.opticalPartials
            state_eval = prop_sim_past.xIntegEval
        elif self.future_obs_exist:
            optical_obs_dot = prop_sim_future.opticalObsDot
            optical_obs_partials = prop_sim_future.opticalPartials
            state_eval = prop_sim_future.xIntegEval
        computed_obs_dot = np.array(optical_obs_dot)[:,0:2]
        computed_obs_dot[:, 0] *= self.obs.cosDec.values
        radius_nonzero = self.fixed_propsim_params['radius'] != 0.0
        fac = 0.0
        if radius_nonzero:
            rel_dists = np.linalg.norm([state[:3] for state in state_eval], axis=1)
            lmbda = 0.3 # from fuentes-munoz et al. 2024
            au2m = 1.495978707e11
            fac = (lmbda*self.fixed_propsim_params['radius']/au2m/rel_dists*180/np.pi*3600)**2
        modes = self.obs['mode'].values
        no_time_uncert = (
            special_codes["gaia"]
            | special_codes["occultation"]
            | special_codes["spacecraft"]
        )
        for i in range(len(self.observer_info_lengths)):
            if modes[i] not in {'RAD', 'SIM_RAD_DEL', 'SIM_RAD_DOP', 'SIM_CCD', 'SIM_OCC'}:
                sig_ra = sig_ra_vals[i]
                sig_dec = sig_dec_vals[i]
                sig_corr = sig_corr_vals[i]
                off_diag = 0.0 if np.isnan(sig_corr) else sig_corr*sig_ra*sig_dec
                cov = np.array([[sig_ra**2, off_diag],
                                [off_diag, sig_dec**2]])
                time_uncert = sig_times[i]
                # apply time uncertainties to the optical observations
                if time_uncert != 0.0 and stations[i] not in no_time_uncert:
                    ra_dot_cos_dec = computed_obs_dot[i, 0]
                    dec_dot = computed_obs_dot[i, 1]
                    off_diag_time = ra_dot_cos_dec*dec_dot
                    cov_time = np.array([[ra_dot_cos_dec**2, off_diag_time],
                                        [off_diag_time, dec_dot**2]])*(time_uncert/86400)**2
                    cov += cov_time
                    sig_times[i] = time_uncert
                # apply Gaia astrometric handling
                if radius_nonzero and stations[i] == '258':
                    cov_fac = np.array([[fac[i], 0.0],
                                        [0.0, fac[i]]])
                    cov += cov_fac
                # apply station position uncertainties
                if np.isfinite(p11_vals[i]):
                    cov_stn_pos = np.array([
                        [p11_vals[i], p12_vals[i], p13_vals[i]],
                        [p12_vals[i], p22_vals[i], p23_vals[i]],
                        [p13_vals[i], p23_vals[i], p33_vals[i]]])
                    # make sure all cov values are finite
                    if not np.isfinite(cov_stn_pos).all():
                        raise ValueError("Station position covariance matrix is not finite.")
                    if sys_vals[i] not in {'ICRF_AU', 'ICRF_KM'}:
                        raise ValueError("Station position covariance matrix system not defined/recognized.")
                    conv = 1 if sys_vals[i] == 'ICRF_AU' else 1/1.495978707e8
                    cov_stn_pos *= conv**2
                    partials = np.array(optical_obs_partials[i]).reshape(2, 6)[:, :3]
                    partials[0, :] *= self.obs.cosDec.values[i]
                    cov_stn_unc = partials @ cov_stn_pos @ partials.T
                    cov += cov_stn_unc
                det = cov[0, 0]*cov[1, 1] - cov[0, 1]*cov[1, 0]
                inv = np.array([[cov[1, 1], -cov[0, 1]],
                                [-cov[1, 0], cov[0, 0]]])/det
                self.obs_cov[i] = cov
                self.obs_weight[i] = inv
        self.obs.sigTime = sig_times
        return None

    def _get_analytic_partials(self, prop_sim_past, prop_sim_future):
        """
        Computes the analytic partials of the observations with respect to the
        initial nominal state.

        Parameters
        ----------
        prop_sim_past : libgrss.PropSimulation object
            The propagated PropSimulation object for the past.
        prop_sim_future : libgrss.PropSimulation object
            The propagated PropSimulation object for the future.
        """
        partials = np.zeros((self.n_obs, self.n_fit))
        len_past_idx = len(self.past_obs_idx) if self.past_obs_exist else 0
        partials_idx = 0
        if self.past_obs_exist:
            past_optical_partials = np.array(prop_sim_past.opticalPartials)
            past_radar_partials = np.array(prop_sim_past.radarPartials)
            past_states = np.array(prop_sim_past.xIntegEval)
        if self.future_obs_exist:
            future_optical_partials = np.array(prop_sim_future.opticalPartials)
            future_radar_partials = np.array(prop_sim_future.radarPartials)
            future_states = np.array(prop_sim_future.xIntegEval)
        cos_dec = self.obs.cosDec.values
        for i, obs_info_len in enumerate(self.observer_info_lengths):
            if obs_info_len in {4, 7}:
                is_optical = True
                size = 2
            elif obs_info_len in {9, 10}:
                is_optical = False
                size = 1
            else:
                raise ValueError("Observer info length not recognized.")
            # part is partial of observable with respect to state at observation time
            part = np.zeros((size, 6))
            if self.past_obs_exist and i < len_past_idx:
                sim_idx = i
                optical_partials = past_optical_partials
                radar_partials = past_radar_partials
                state = past_states
            else:
                sim_idx = i-len_past_idx
                optical_partials = future_optical_partials
                radar_partials = future_radar_partials
                state = future_states
            # stm is partial of state at observation time with respect to nominal state
            stm = libgrss.reconstruct_stm(state[sim_idx][6:])[:6]
            if is_optical:
                part[0, :6] = optical_partials[sim_idx, :6]*cos_dec[i]
                part[1, :6] = optical_partials[sim_idx, 6:12]
            else:
                part[0, :6] = radar_partials[sim_idx, :6]
            # partial is partial of observable with respect to nominal state
            partial = part @ stm
            partials[partials_idx:partials_idx+size, :partial.shape[1]] = partial
            partials_idx += size
        return partials

    def _get_numeric_partials(self, prop_sim_past, prop_sim_future, perturbation_info):
        """
        Computes the numeric partials of the observations with respect to the
        initial nominal state.

        Parameters
        ----------
        prop_sim_past : libgrss.PropSimulation object
            The propagated PropSimulation object for the past.
        prop_sim_future : libgrss.PropSimulation object
            The propagated PropSimulation object for the future.
        perturbation_info : list
            A list of tuples containing the perturbation information
            for each nominal state parameter.

        Returns
        -------
        partials : array
            The numeric partials of the observations with respect to the
            initial nominal state.
        """
        partials = np.zeros((self.n_obs, self.n_fit))
        for i in range(self.n_fit):
            _ = list(self.x_nom.keys())[i]
            _, _, _, _, _, _, fd_delta = perturbation_info[i]
            # get computed_obs for perturbed states
            computed_obs_plus = self._get_computed_obs(prop_sim_past, prop_sim_future,
                                                        integ_body_idx=2*i+1)
            computed_obs_minus = self._get_computed_obs(prop_sim_past, prop_sim_future,
                                                        integ_body_idx=2*i+2)
            computed_obs_plus = self._flatten_and_clean(computed_obs_plus)
            computed_obs_minus = self._flatten_and_clean(computed_obs_minus)
            # get partials
            partials[:, i] = (computed_obs_plus - computed_obs_minus)/(2*fd_delta)
        return partials

    def _get_partials(self, prop_sim_past, prop_sim_future, perturbation_info):
        """
        Computes the partials of the observations with respect to the
        initial nominal state.

        Parameters
        ----------
        prop_sim_past : libgrss.PropSimulation object
            The propagated PropSimulation object for the past.
        prop_sim_future : libgrss.PropSimulation object
            The propagated PropSimulation object for the future.
        perturbation_info : list
            A list of tuples containing the perturbation information
            for each nominal state parameter.

        Returns
        -------
        partials : array
            The partials of the observations with respect to the
            initial nominal state.
        """
        if self.analytic_partials:
            return self._get_analytic_partials(prop_sim_past, prop_sim_future)
        return self._get_numeric_partials(prop_sim_past, prop_sim_future, perturbation_info)

    def _get_residuals_and_partials(self):
        """
        Computes the residuals and partials of the observations with respect to the
        initial nominal state.

        Returns
        -------
        residuals : array
            The residuals of the observations
        partials : array
            The partials of the observations with respect to the
            initial nominal state.
        """
        perturbation_info = None if self.analytic_partials else self._get_perturbation_info()
        prop_sim_past, prop_sim_future = self._assemble_and_propagate_bodies(perturbation_info)
        self.prop_sims = (prop_sim_past, prop_sim_future)
        # get partials
        partials = self._get_partials(prop_sim_past, prop_sim_future, perturbation_info)
        # get residuals
        computed_obs = self._get_computed_obs(prop_sim_past, prop_sim_future, integ_body_idx=0)
        residuals = [None]*len(self.obs)
        ra_res = self.obs['resRA'].values
        dec_res = self.obs['resDec'].values
        delay_res = self.obs['resDelay'].values
        doppler_res = self.obs['resDoppler'].values
        fields = ['mode', 'ra', 'dec', 'cosDec', 'biasRA', 'biasDec', 'delay', 'doppler']
        for info in self.obs[fields].itertuples():
            i, mode, ra, dec, cos_dec, bias_ra, bias_dec, delay, doppler = info
            if mode.startswith('SIM'):
                if mode in {'SIM_RAD_DEL', 'SIM_RAD_DOP'}:
                    residuals[i] = np.array([0.0])
                else:
                    residuals[i] = np.array([0.0, 0.0])
            elif mode == 'RAD':
                if i in self.delay_idx:
                    delay_res[i] = delay*1e6 - computed_obs[i, 0]
                    residuals[i] = np.array([delay_res[i]])
                else:
                    doppler_res[i] = doppler - computed_obs[i, 1]
                    residuals[i] = np.array([doppler_res[i]])
            else:
                ra_res[i] = ra*3600*cos_dec - bias_ra - computed_obs[i, 0]
                dec_res[i] = dec*3600 - bias_dec - computed_obs[i, 1]
                residuals[i] = np.array([ra_res[i], dec_res[i]])
        self.obs['resRA'] = ra_res
        self.obs['resDec'] = dec_res
        self.obs['resDelay'] = delay_res
        self.obs['resDoppler'] = doppler_res
        return residuals, partials

    def _get_rms_and_reject_outliers(self, partials, residuals, start_rejecting):
        # sourcery skip: low-code-quality
        """
        Outlier rejection algorithm for the residuals.

        Parameters
        ----------
        partials : array
            The partials of the observations with respect to the
            initial nominal state.
        residuals : array
            The residuals of the observations
        start_rejecting : bool
            Flag for whether to start rejecting outliers.

        Returns
        -------
        rms_u : float
            Unweighted RMS of the residuals.
        rms_w : float
            Weighted RMS of the residuals.
        chi_sq : float
            Chi-squared value of the residuals.

        Raises
        ------
        ValueError
            If chi_recover is greater than chi_reject.
        ValueError
            If the observer information is not well-defined.
        """
        chi_reject = self.reject_criteria[0]
        chi_recover = self.reject_criteria[1]
        if chi_recover > chi_reject:
            raise ValueError("chi_recover must be less than chi_reject. "
                                "Use default values if unsure. "
                                "Default values are chi_reject=3.0 and chi_recover=2.8 "
                                "(Implemented as FitSimulation.reject_criteria=[3.0, 2.8])")
        full_cov = self.covariance
        rms_u = 0
        chi_sq = 0
        j = 0
        sel_ast = self.obs['selAst'].values
        res_chisq_vals = np.nan*np.ones(len(self.observer_info_lengths))
        for i, obs_info_len in enumerate(self.observer_info_lengths):
            if obs_info_len in {4, 7}:
                size = 2
            elif obs_info_len in {9, 10}:
                size = 1
            else:
                raise ValueError("Observer info length not recognized.")
            resid = residuals[i]
            # calculate chi-squared for each residual if after the first iteration
            if start_rejecting and size == 2:
                obs_cov = self.obs_cov[i]
                obs_partials = partials[j:j+size, :]
                if sel_ast[i] in {'D', 'd'}:
                    resid_cov = obs_cov + obs_partials @ full_cov @ obs_partials.T
                else:
                    resid_cov = obs_cov - obs_partials @ full_cov @ obs_partials.T
                resid_cov_det = resid_cov[0,0]*resid_cov[1,1] - resid_cov[0,1]*resid_cov[1,0]
                resid_cov_inv = np.array([[resid_cov[1,1], -resid_cov[0,1]],
                                            [-resid_cov[1,0], resid_cov[0,0]]])/resid_cov_det
                outlier_chisq = resid @ resid_cov_inv @ resid.T
                # outlier rejection, only reject RA/Dec measurements
                if abs(outlier_chisq) > chi_reject**2 and sel_ast[i] not in {'a', 'd'}:
                    if sel_ast[i] == 'A':
                        self.num_rejected += 1
                    sel_ast[i] = 'D'
                elif abs(outlier_chisq) < chi_recover**2 and sel_ast[i] == 'D':
                    sel_ast[i] = 'A'
                    self.num_rejected -= 1
            res_chisq_vals[i] = resid @ self.obs_weight[i] @ resid.T
            if sel_ast[i] not in {'D', 'd'}:
                rms_u += resid @ resid.T
                chi_sq += res_chisq_vals[i]
            j += size
        # # write res_chisq_vals to file if any values are negative
        # if np.any(res_chisq_vals < 0):
        #     print("WARNING: Negative outlier rejection chi-squared values detected. "
        #             "Writing to file 'warn_neg_chisq_vals.txt'.")
        #     np.savetxt("warn_neg_chisq_vals.txt", res_chisq_vals, fmt="%.16f")
        self.obs['selAst'] = sel_ast
        self.obs['resChi'] = res_chisq_vals**0.5
        rejected_fraction = self.num_rejected/len(self.obs)
        if rejected_fraction > 0.25:
            print("WARNING: More than 25% of observations rejected. Consider changing the",
                    "rejection criteria, or turning off outlier rejection altogether.")
        rms_u = np.sqrt(rms_u/self.n_obs)
        rms_w = np.sqrt(chi_sq/self.n_obs)
        return rms_u, rms_w, chi_sq

    def _add_iteration(self, iter_number, rms_u, rms_w, chi_sq):
        """
        Adds an iteration to the list of iterations in the FitSimulation object.

        Parameters
        ----------
        iter_number : int
            Iteration number.
        rms_u : float
            Unweighted RMS of the residuals.
        rms_w : float
            Weighted RMS of the residuals.
        chi_sq : float
            Chi-squared value of the residuals.

        Returns
        -------
        None : NoneType
            None
        """
        self.iters.append(IterationParams(iter_number, self.x_nom, self.covariance,
                                            self.obs, rms_u, rms_w, chi_sq))
        return None

    def _check_convergence(self, delta_x):
        """
        Checks if the orbit fit has converged.

        Parameters
        ----------
        delta_x : array
            The state correction.

        Returns
        -------
        None : NoneType
            None
        """
        # check for convergence based on weighted rms
        if self.n_iter > 1:
            del_rms_convergence = 1e-4
            curr_rms = self.iters[-1].weighted_rms
            prev_rms = self.iters[-2].weighted_rms
            del_rms = abs(prev_rms - curr_rms)#/prev_rms
            if del_rms < del_rms_convergence:
                # print("Converged based on weighted RMS.")
                self.converged = True
        # check for convergence based on magnitude of corrections
        sigmas = np.sqrt(np.diag(self.covariance))
        corrections = np.abs(delta_x/sigmas)
        max_correction = np.max(corrections)
        if max_correction < 1e-2:
            # print("Converged based on magnitude of corrections.")
            self.converged = True
        return None

    def _get_lsq_state_correction(self, partials, residuals):
        """
        Get the state correction using least-squares.

        Parameters
        ----------
        partials : array
            The partials of the observations with respect to the
            initial nominal state.
        residuals : array
            The residuals of the observations

        Returns
        -------
        delta_x : array
            The state correction.
        """
        if self._priors_given:
            atwa = self._info0.copy()
            atwb = (self._info0 @ self._xbar0).copy()
        else:
            atwa = np.zeros((self.n_fit, self.n_fit))
            atwb = np.zeros(self.n_fit)
        j = 0
        sel_ast = self.obs['selAst'].values
        self.info_mats = []
        for i, obs_info_len in enumerate(self.observer_info_lengths):
            if obs_info_len in {4, 7}:
                size = 2
            elif obs_info_len in {9, 10}:
                size = 1
            else:
                raise ValueError("Observer info length not recognized.")
            if sel_ast[i] in {'D', 'd'}:
                j += size
                self.info_mats.append(atwa.copy())
                continue
            atwa += partials[j:j+size, :].T @ self.obs_weight[i] @ partials[j:j+size, :]
            atwb += partials[j:j+size, :].T @ self.obs_weight[i] @ residuals[i]
            j += size
            self.info_mats.append(atwa.copy())
        # use pseudo-inverse if the data arc is less than 7 days
        data_arc = self.obs.obsTimeMJD.max() - self.obs.obsTimeMJD.min()
        if data_arc < 1.0:
            self.covariance = np.linalg.pinv(atwa, rcond=1e-10, hermitian=True)
        elif data_arc < 7.0:
            self.covariance = np.linalg.pinv(atwa, rcond=1e-20, hermitian=True)
        else:
            self.covariance = np.array(libgrss.matrix_inverse(atwa))
        delta_x = self.covariance @ atwb
        return delta_x.ravel()

    def filter_lsq(self, verbose=True):
        """
        Performs a least-squares fit on the observations.

        Parameters
        ----------
        verbose : bool, optional
            Flag for printing the iteration information while fitting, by default True.

        Returns
        -------
        None : NoneType
            None
        """
        self._check_priors()
        start_rejecting = False
        if verbose:
            print("Iteration\t\tUnweighted RMS\t\tWeighted RMS",
                    "\t\tChi-squared\t\tReduced Chi-squared")
        delta_x = np.zeros(self.n_fit)
        for i in range(self.n_iter_max):
            if self._priors_given:
                self._xbar0 -= delta_x
                a_priori_constant = list(self.x_nom.values())+self._xbar0
                a_priori_constant_violation = np.linalg.norm(a_priori_constant-self._prior_constant)
                if a_priori_constant_violation >= 1e-11:
                    print(f"a priori constraint constant violation before iteration {i+1}: "
                            f"{a_priori_constant_violation}. Solution may not be trustworthy. "
                            "Violation for each estimated parameter:")
                    print(dict(zip(self.x_nom.keys(), a_priori_constant-self._prior_constant)))
            self.n_iter = i+1
            # get residuals and partials
            residuals, partials = self._get_residuals_and_partials()
            # calculate rms and reject outliers here if desired
            rms_u, rms_w, chi_sq = self._get_rms_and_reject_outliers(partials, residuals,
                                                                            start_rejecting)
            # get current state
            curr_state = list(self.x_nom.values())
            # get state correction
            delta_x = self._get_lsq_state_correction(partials, residuals)
            if self.constraint_dir is not None:
                constr_hat = self.constraint_dir/np.linalg.norm(self.constraint_dir)
                delta_x -= np.dot(delta_x, constr_hat)*constr_hat
            # make sure eccentricity is non-negative and
            # any keys starting with dt are 0.1 at minimum
            for key in self.x_nom:
                if self.fit_cometary and key == 'e':
                    idx = list(self.x_nom.keys()).index(key)
                    if curr_state[idx] + delta_x[idx] < 0.0:
                        delta_x[idx] = -curr_state[idx]
                        print("WARNING: Eccentricity is negative per least squares "
                                "state correction. Setting to 0.0. This solution "
                                "may not be trustworthy.")
                if key[:2] == 'dt':
                    idx = list(self.x_nom.keys()).index(key)
                    if curr_state[idx] + delta_x[idx] < 0.1:
                        delta_x[idx] = 0.1 - curr_state[idx]
            # get new nominal state
            next_state = curr_state + delta_x
            self.x_nom = dict(zip(self.x_nom.keys(), next_state))
            # add iteration
            self._add_iteration(i+1, rms_u, rms_w, chi_sq)
            if verbose:
                print(f"{self.iters[-1].iter_number}\t\t\t",
                        f"{self.iters[-1].unweighted_rms:.3f}\t\t\t",
                        f"{self.iters[-1].weighted_rms:.3f}\t\t\t",
                        f"{self.iters[-1].chi_squared:.3f}\t\t\t",
                        f"{self.iters[-1].reduced_chi_squared:.3f}")
            self._check_convergence(delta_x)
            if self.converged:
                if self.reject_outliers and start_rejecting:
                    if verbose:
                        print(f"Converged after rejecting outliers. Rejected {self.num_rejected}",
                                f"out of {len(self.optical_idx)} optical observations.")
                    break
                msg = "Converged without rejecting outliers."
                if self.reject_outliers:
                    msg += " Starting outlier rejection now..."
                    start_rejecting = True
                    self.converged = False
                if verbose:
                    print(msg)
                if not self.reject_outliers:
                    break
        # add postfit iteration if converged
        if self.converged:
            residuals, partials = self._get_residuals_and_partials()
            rms_u, rms_w, chi_sq = self._get_rms_and_reject_outliers(partials, residuals,
                                                                        start_rejecting)
            self._add_iteration(self.n_iter+1, rms_u, rms_w, chi_sq)
        if self.n_iter == self.n_iter_max and not self.converged:
            print("WARNING: Maximum number of iterations reached without converging.")
        return None

    def print_summary(self, iter_idx=None, out_stream=sys.stdout):
        """
        Prints a summary of the orbit fit calculations at a given iteration.

        Parameters
        ----------
        iter_idx : int, optional
            The iteration number to print the summary for, by default -1.

        Returns
        -------
        None : NoneType
            None
        """
        iter_idx = len(self.iters)-1 if iter_idx is None else iter_idx
        data = self.iters[iter_idx]
        arc = self.obs.obsTimeMJD.max() - self.obs.obsTimeMJD.min()
        str_to_print = "Summary of the orbit fit calculations"
        if data.iter_number <= self.n_iter:
            str_to_print += f" after iteration {data.iter_number} (of {self.n_iter}):\n"
        else:
            str_to_print += " after postfit pass:\n"
        str_to_print += "==============================================================\n"
        str_to_print += f"RMS unweighted: {data.unweighted_rms}\n"
        str_to_print += f"RMS weighted: {data.weighted_rms}\n"
        str_to_print += f"chi-squared: {data.chi_squared}\n"
        str_to_print += f"reduced chi-squared: {data.reduced_chi_squared}\n"
        str_to_print += f"square root of reduced chi-squared: {np.sqrt(data.reduced_chi_squared)}\n"
        str_to_print += "--------------------------------------------------------------\n"
        str_to_print += f"Solution Time: MJD {self.t_sol:0.3f} TDB = "
        str_to_print += f"{Time(self.t_sol, format='mjd', scale='tdb').iso} TDB\n"
        str_to_print += f"Solution Observation Arc: {arc:0.2f} days ({arc/365.25:0.2f} years)\n"
        str_to_print += "--------------------------------------------------------------\n"
        str_to_print += "Fitted Variable\t\tInitial Value\t\t\tUncertainty\t\t\tFitted Value"
        str_to_print += "\t\t\tUncertainty\t\t\tChange\t\t\t\tChange (sigma)\n"
        init_variance = np.sqrt(np.diag(self.covariance_init))
        final_variance = np.sqrt(np.diag(self.covariance))
        init_sol = self.x_init
        final_sol = data.x_nom
        with np.errstate(divide='ignore'):
            for i, key in enumerate(init_sol.keys()):
                init_val = init_sol[key]
                init_unc = init_variance[i]
                final_val = final_sol[key]
                final_unc = final_variance[i]
                if self.fit_cometary and key in ['om', 'w', 'i']:
                    init_val *= 180/np.pi
                    init_unc *= 180/np.pi
                    final_val *= 180/np.pi
                    final_unc *= 180/np.pi
                str_to_print += f"{key}\t\t\t{init_val:.11e}\t\t{init_unc:.11e}"
                str_to_print += f"\t\t{final_val:.11e}\t\t{final_unc:.11e}"
                str_to_print += f"\t\t{final_val-init_val:+.11e}"
                str_to_print += f"\t\t{(final_val-init_val)/final_unc:+.3f}\n"
        print(str_to_print, file=out_stream)
        return None

    def plot_summary(self, auto_close=False):
        """
        Plots a summary of the orbit fit calculations.

        Parameters
        ----------
        auto_close : bool, optional
            Flag to automatically close the figure, by default False

        Returns
        -------
        None : NoneType
            None
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.labelsize'] = 12
        ticks = np.arange(1, self.n_iter+2, 1)
        labels = list(ticks)
        labels[-1] = "Postfit"
        start_idx = 1
        iters_for_plot = self.iters[start_idx:]
        final_iter = iters_for_plot[-1]
        plt.figure(figsize=(21,10), dpi=150)
        plt.subplot(2, 2, 1)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot],
                        [iteration.unweighted_rms for iteration in iters_for_plot],
                        label=f"Final Unweighted RMS={final_iter.unweighted_rms:.3e}")
        plt.xticks(ticks)
        plt.gca().set_xticklabels(labels)
        plt.xlabel("Iteration")
        plt.ylabel("Unweighted RMS")
        plt.legend()
        plt.subplot(2, 2, 2)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot],
                        [iteration.weighted_rms for iteration in iters_for_plot],
                        label=f"Final Weighted RMS={final_iter.weighted_rms:.3e}")
        plt.xticks(ticks)
        plt.gca().set_xticklabels(labels)
        plt.xlabel("Iteration")
        plt.ylabel("Weighted RMS")
        plt.legend()
        plt.subplot(2, 2, 3)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot],
                        [iteration.chi_squared for iteration in iters_for_plot],
                        label=fr"Final $\chi^2$={final_iter.chi_squared:.3e}")
        plt.xticks(ticks)
        plt.gca().set_xticklabels(labels)
        plt.xlabel("Iteration")
        plt.ylabel(r"$\chi^2$")
        plt.legend()
        plt.subplot(2, 2, 4)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot],
                        [iteration.reduced_chi_squared for iteration in iters_for_plot],
                        label=fr"Final Reduced $\chi^2$={final_iter.reduced_chi_squared:.3e}")
        plt.xticks(ticks)
        plt.gca().set_xticklabels(labels)
        plt.xlabel("Iteration")
        plt.ylabel(r"Reduced $\chi^2$")
        plt.legend()
        plt.tight_layout()
        block = not auto_close
        plt.show(block=block)
        if auto_close:
            plt.close()
        return None

    def save(self, filename):
        """
        Saves the FitSimulation log to a file.

        Parameters
        ----------
        filename : str
            The filename to save the log to.

        Returns
        -------
        None : NoneType
            None
        """
        max_chars = 135
        half_tab = "  "
        subsection_full = "-"*max_chars
        section_full = "="*max_chars
        header_section_half = "="*((max_chars-13)//2)
        try:
            with open(f'{grss_path}/version.txt', 'r', encoding='utf-8') as f:
                version = f.read().strip()
        except FileNotFoundError:
            version = "INFTY"
        units_dict = {
            'e': '', 'q': 'AU', 'tp': 'MJD [TDB]', 'om': 'deg',
            'w': 'deg', 'i': 'deg', 'x': 'AU', 'y': 'AU', 'z': 'AU',
            'vx': 'AU/day', 'vy': 'AU/day', 'vz': 'AU/day',
            'a1': 'AU/day^2', 'a2': 'AU/day^2', 'a3': 'AU/day^2',
            'dv': 'AU/day', 'ax': 'AU/day^2', 'ay': 'AU/day^2', 'az': 'AU/day^2',
            'dt': 'day'
        }
        with open(filename, 'x', encoding='utf-8') as f:
            f.write(section_full + '\n')
            f.write(f"{header_section_half} GRSS v{version} {header_section_half}" + '\n')
            f.write(section_full + '\n')

            dt = datetime.datetime.now(datetime.timezone.utc)
            time_str = dt.strftime("%a %b %d %Y %H:%M:%S UTC")
            time_buff = " "*int((max_chars-len(time_str))/2)
            f.write(f'{time_buff}{time_str}{time_buff}\n\n')
            name_buff = "-"*int((max_chars-len(self.name))/2)
            f.write(f'{name_buff}{self.name}{name_buff}\n')
            f.write(subsection_full + '\n')
            obs_min = self.obs.obsTimeMJD.argmin()
            obs_max = self.obs.obsTimeMJD.argmax()
            f.write(f"Observation data arc from {self.obs.loc[obs_min, 'obsTime']} to ")
            f.write(f"{self.obs.loc[obs_max, 'obsTime']} [UTC] ")
            range_mjd = self.obs.obsTimeMJD.max() - self.obs.obsTimeMJD.min()
            f.write(f"({range_mjd:.2f} days, {range_mjd/365.25:.2f} years)\n")
            f.write(f"{len(self.obs)} observations ")
            f.write(f"({len(self.optical_idx)} optical, {len(self.delay_idx)} delay, ")
            f.write(f"{len(self.doppler_idx)} doppler)\n")
            f.write(f"DE kernel version: {self.de_kernel}\n")
            f.write("\n")
            init_state_str = "Initial Nominal State"
            init_state_buff = "-"*int((max_chars-len(init_state_str))/2)
            f.write(f'{init_state_buff}{init_state_str}{init_state_buff}\n')
            f.write(subsection_full + '\n')
            f.write(f"T: {Time(self.t_sol, format='mjd', scale='tdb').iso} ")
            f.write(f"({self.t_sol:.8f} MJD) [TDB]\n")
            keys = list(self.x_init.keys())
            init_variance = np.sqrt(np.diag(self.covariance_init))
            for i, key in enumerate(keys):
                val = self.x_init[key]
                sig = init_variance[i]
                if key in ['om', 'w', 'i']:
                    val *= 180/np.pi
                    sig *= 180/np.pi
                if key[:2] in {'a1', 'a2', 'a3', 'dv', 'ax', 'ay', 'az'}:
                    f.write(f"{key.upper()}: {val:.11e} \u00B1 {sig:.11e}")
                else:
                    f.write(f"{key.upper()}: {val:.11f} \u00B1 {sig:.11f}")
                if key[:2] in units_dict:
                    f.write(f" {units_dict[key[:2]]}")
                if self._priors_given and key in self.prior_est:
                    p_val = self.prior_est[key]
                    p_sig = self.prior_sig[key]
                    if key in ['om', 'w', 'i']:
                        p_val *= 180/np.pi
                        p_sig *= 180/np.pi
                    if key[:2] in {'a1', 'a2', 'a3', 'dv', 'ax', 'ay', 'az'}:
                        f.write(f" (a priori: {p_val:.11e} \u00B1 {p_sig:.11e}")
                    else:
                        f.write(f" (a priori: {p_val:.11f} \u00B1 {p_sig:.11f}")
                    if key[:2] in units_dict:
                        f.write(f" {units_dict[key[:2]]}")
                    f.write(")")
                f.write("\n")
            f.write("\n")
            f.write("Initial Covariance Matrix:\n")
            for row in self.covariance_init:
                f.write(half_tab.join([f"{val:18.11e}" for val in row]) + '\n')
            f.write("\n")

            final_state_str = "Final Nominal State"
            final_state_buff = "-"*int((max_chars-len(final_state_str))/2)
            f.write(f'{final_state_buff}{final_state_str}{final_state_buff}\n')
            f.write(subsection_full + '\n')
            f.write(f"T: {Time(self.t_sol, format='mjd', scale='tdb').iso} ")
            f.write(f"({self.t_sol:.8f} MJD) [TDB]\n")
            keys = list(self.x_nom.keys())
            final_variance = np.sqrt(np.diag(self.covariance))
            for i, key in enumerate(keys):
                val = self.x_nom[key]
                sig = final_variance[i]
                if key in ['om', 'w', 'i']:
                    val *= 180/np.pi
                    sig *= 180/np.pi
                if key[:2] in {'a1', 'a2', 'a3', 'dv', 'ax', 'ay', 'az'}:
                    f.write(f"{key.upper()}: {val:.11e} \u00B1 {sig:.11e}")
                else:
                    f.write(f"{key.upper()}: {val:.11f} \u00B1 {sig:.11f}")
                if key[:2] in units_dict:
                    f.write(f" {units_dict[key[:2]]}")
                if self._priors_given and key in self.prior_est:
                    p_val = self.prior_est[key]
                    p_sig = self.prior_sig[key]
                    if key in ['om', 'w', 'i']:
                        p_val *= 180/np.pi
                        p_sig *= 180/np.pi
                    if key[:2] in {'a1', 'a2', 'a3', 'dv', 'ax', 'ay', 'az'}:
                        f.write(f" (a priori: {p_val:.11e} \u00B1 {p_sig:.11e}")
                    else:
                        f.write(f" (a priori: {p_val:.11f} \u00B1 {p_sig:.11f}")
                    if key[:2] in units_dict:
                        f.write(f" {units_dict[key[:2]]}")
                    f.write(")")
                f.write("\n")
            f.write("\n")
            f.write("Final Covariance Matrix:\n")
            for row in self.covariance:
                f.write(half_tab.join([f"{val:18.11e}" for val in row]) + '\n')
            f.write("\n")

            postfit_summary_str = "Postfit Summary"
            postfit_summary_buff = "-"*int((max_chars-len(postfit_summary_str))/2)
            f.write(f'{postfit_summary_buff}{postfit_summary_str}{postfit_summary_buff}\n')
            f.write(subsection_full + '\n')
            widths = {
                'obsTime': 29, 'stn': 7, 'ra': 16, 'dec': 18,
                'sigRA': 16, 'sigDec': 17, 'resRA': 16, 'resDec': 15,
            }
            f.write(f'|{"Observations".center(69)}|')
            f.write(f'{"Sigmas".center(32)}|')
            f.write(f'{"Residuals".center(29)}|')
            f.write("\n")
            f.write(f'{"Time [UTC]".center(widths["obsTime"])}')
            f.write(f'{"Obs".center(widths["stn"])}')
            f.write(f'{"RA[deg]/Del[s]".center(widths["ra"])}')
            f.write(f'{"Dec[deg]/Dop[Hz]".center(widths["dec"])}')
            f.write(f'{f"RA[as]/Del[{chr(0x00b5)}s]".center(widths["sigRA"])}')
            f.write(f'{"Dec[as]/Dop[Hz]".center(widths["sigDec"])}')
            f.write(f'{f"RA[as]/Del[{chr(0x00b5)}s]".center(widths["resRA"])}')
            f.write(f'{"Dec[as]/Dop[Hz]".center(widths["resDec"])}')
            f.write("\n")
            for i, obs in self.obs[::-1].iterrows():
                f.write(f'{obs["selAst"]+obs["obsTime"]:<{widths["obsTime"]}}')
                if obs['mode'] in {'RAD', 'SIM_RAD_DEL', 'SIM_RAD_DOP'}:
                    f.write(f'{obs["trx"]+"/"+obs["rcv"]:^{widths["stn"]}}')
                    f.write(f'{obs["delay"]:^{widths["ra"]}.13g}')
                    f.write(f'{obs["doppler"]:^+{widths["dec"]}.12g}')
                    f.write(f'{obs["sigDelay"]:^{widths["sigRA"]}.4f}')
                    f.write(f'{obs["sigDoppler"]:^{widths["sigDec"]}.4f}')
                    f.write(f'{obs["resDelay"]:^+{widths["resRA"]}.7f}')
                    f.write(f'{obs["resDoppler"]:^+{widths["resDec"]}.7f}')
                else:
                    f.write(f'{obs["stn"]:^{widths["stn"]}}')
                    f.write(f'{obs["ra"]:^{widths["ra"]}.13g}')
                    f.write(f'{obs["dec"]:^+{widths["dec"]}.12g}')
                    f.write(f'{obs["sigRA"]:^{widths["sigRA"]}.4f}')
                    f.write(f'{obs["sigDec"]:^{widths["sigDec"]}.4f}')
                    f.write(f'{obs["resRA"]:^+{widths["resRA"]}.7f}')
                    f.write(f'{obs["resDec"]:^+{widths["resDec"]}.7f}')
                f.write("\n")
            f.write("\n")

            iter_summary_str = "Iteration Summary"
            iter_summary_buff = "-"*int((max_chars-len(iter_summary_str))/2)
            f.write(f'{iter_summary_buff}{iter_summary_str}{iter_summary_buff}\n')
            f.write(subsection_full + '\n')
            for itrn in self.iters[1:]:
                f.write("Summary of the orbit fit calculations")
                if itrn.iter_number <= self.n_iter:
                    f.write(f" after iteration {itrn.iter_number} (of {self.n_iter}):\n")
                else:
                    f.write(" after postfit pass:\n")
                f.write(f"RMS unweighted: {itrn.unweighted_rms}\n")
                f.write(f"RMS weighted: {itrn.weighted_rms}\n")
                f.write(f"chi-squared: {itrn.chi_squared}\n")
                f.write(f"reduced chi-squared: {itrn.reduced_chi_squared}\n")
                f.write(f"square root of reduced chi-squared: {itrn.reduced_chi_squared**0.5}\n\n")
                f.write(f'{"Variable":<8}{"Initial Value":>20}')
                f.write(f'{"Uncertainty":>20}{"Fitted Value":>20}')
                f.write(f'{"Uncertainty":>20}{"Change":>20}{"Change":>8}'+' (\u03C3)\n')
                init_sol = self.x_init
                final_sol = itrn.x_nom
                with np.errstate(divide='ignore'):
                    for i, key in enumerate(init_sol.keys()):
                        init_val = init_sol[key]
                        init_unc = init_variance[i]
                        final_val = final_sol[key]
                        final_unc = final_variance[i]
                        if key in ['om', 'w', 'i']:
                            init_val *= 180/np.pi
                            init_unc *= 180/np.pi
                            final_val *= 180/np.pi
                            final_unc *= 180/np.pi
                        f.write(f'{key.upper():<8}{half_tab}')
                        f.write(f'{init_val:18.11e}{half_tab}')
                        f.write(f'{init_unc:18.11e}{half_tab}')
                        f.write(f'{final_val:18.11e}{half_tab}')
                        f.write(f'{final_unc:18.11e}{half_tab}')
                        f.write(f'{final_val-init_val:+18.11e}{half_tab}')
                        f.write(f'{(final_val-init_val)/final_unc:+10.3f}\n')
                f.write("\n")
                if itrn.iter_number <= self.n_iter:
                    f.write(subsection_full + '\n')

            mean_0 = np.array(list(self.x_init.values()))
            cov_0 = self.covariance_init
            mean_f = np.array(list(self.x_nom.values()))
            cov_f = self.covariance
            md_f, md_0, b_dist, b_coeff = get_similarity_stats(mean_0, cov_0, mean_f, cov_f)
            f.write(f'Mahalonobis distance between initial and final solution: {md_f:0.2f}')
            f.write("\n")
            f.write(f'Mahalonobis distance between final and initial solution: {md_0:0.2f}')
            f.write("\n")
            f.write(f'Bhattacharya distance between initial and final solution: {b_dist:0.4f}')
            f.write("\n")
            f.write(f'Bhattacharya coefficient between initial and final solution: {b_coeff:0.4f}')
            f.write("\n")
            f.write(subsection_full + '\n')

            f.write(section_full + '\n')
            f.write(f"{header_section_half} END OF FILE {header_section_half}" + '\n')
            f.write(section_full + '\n')
        return None
