"""Simulation classes for the Python GRSS orbit determination code"""
#pylint: disable=no-member, too-many-lines
import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice

from . import prop
from .fit_utils import get_observer_info, get_radec

__all__ = [ 'FitSimulation',
            'create_simulated_obs_arrays',
]

class IterationParams:
    """
    Class for storing the iteration parameters for a single iteration of the
    orbit determination process. It is also used for plotting the residuals
    and chi-squared values for each iteration.
    """
    def __init__(self, iter_number, x_nom, covariance, residuals,
                    obs_array, obs_weight_rej, observer_codes, rejection_flags):
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
        residuals : array
            Residuals at the current iteration
        obs_array : array
            Observation array for the orbit fit
        obs_weight_rej : array
            Observation weight matrix after outlier rejection
        observer_codes : tuple
            Observer locations for each observation in obs_array
        rejection_flags : list
            List of rejection flags for each observation in obs_array
        """
        self.iter_number = iter_number
        self.x_nom = x_nom
        self.covariance = covariance
        variance = np.sqrt(np.diag(self.covariance))
        self.variance = dict(zip([f'var_{k}' for k in self.x_nom.keys()], variance))
        self.residuals = residuals
        self.obs_array = obs_array
        self.observer_codes = observer_codes
        self.is_accepted = np.where(np.logical_not(rejection_flags))[0]
        self.is_rejected = np.where(rejection_flags)[0]
        self.sigmas = obs_array[:, 3:5]
        self.obs_weight_rej = obs_weight_rej
        self.calculate_rms()
        self.calculate_chis()
        self.assemble_info()
        return None

    def flatten_and_clean(self, arr):
        """
        Flatten an array and remove NaNs

        Parameters
        ----------
        arr : array
            Array to be flattened and cleaned

        Returns
        -------
        arr : array
            Flattened and cleaned array
        """
        arr = arr.flatten()
        arr = arr[~np.isnan(arr)]
        return arr

    def calculate_rms(self):
        """
        Calculate the weighted and unweighted RMS values for the residuals

        Returns
        -------
        None : NoneType
            None
        """
        resid_arr = self.flatten_and_clean(self.residuals[self.is_accepted])
        n_obs = len(resid_arr)
        self.unweighted_rms = float(np.sqrt(resid_arr.T @ resid_arr/n_obs))
        self.weighted_rms = float(np.sqrt(resid_arr.T @ self.obs_weight_rej @ resid_arr/n_obs))
        return None

    def calculate_chis(self):
        """
        Calculates the chi and chi-squared values for the residuals

        Returns
        -------
        None : NoneType
            None
        """
        resid_arr = self.flatten_and_clean(self.residuals[self.is_accepted])
        n_obs = len(resid_arr)
        n_fit = len(self.x_nom)
        sigmas = self.flatten_and_clean(self.sigmas[self.is_accepted])
        self.chi = resid_arr/sigmas
        self.chi_squared = np.sum(self.chi**2)
        self.reduced_chi_squared = self.chi_squared/(n_obs-n_fit)
        return None

    def assemble_info(self):
        """
        Assembles the information for the iteration into a dictionary

        Returns
        -------
        None : NoneType
            None
        """
        # sourcery skip: extract-duplicate-method
        # optical_idx is where neither the RA nor dec residuals are NaN
        optical_idx = np.where(~np.isnan(self.residuals[:, 0]) &
                                ~np.isnan(self.residuals[:, 1]))[0]
        # radar_idx is where either the RA or dec residuals are NaN
        radar_idx = np.where(np.isnan(self.residuals[:, 0]) |
                                np.isnan(self.residuals[:, 1]))[0]
        optical_info = self.obs_array.copy()
        optical_info[radar_idx, :] = np.nan
        ra_obs = optical_info[:, 1]
        dec_obs = optical_info[:, 2]
        ra_noise = optical_info[:, 3]
        dec_noise = optical_info[:, 4]
        optical_residuals = self.residuals.copy()
        optical_residuals[radar_idx, :] = np.nan
        radar_info = self.obs_array.copy()
        radar_info[optical_idx, :] = np.nan
        delay_obs = radar_info[:, 1]
        doppler_obs = radar_info[:, 2]
        delay_noise = radar_info[:, 3]
        doppler_noise = radar_info[:, 4]
        radar_residuals = self.residuals.copy()
        radar_residuals[optical_idx, :] = np.nan
        self.all_info = {
            'ra_res': optical_residuals[:, 0],
            'dec_res': optical_residuals[:, 1],
            'ra_obs': ra_obs,
            'dec_obs': dec_obs,
            'ra_noise': ra_noise,
            'dec_noise': dec_noise,
            'ra_comp': ra_obs - optical_residuals[:, 0],
            'dec_comp': dec_obs - optical_residuals[:, 1],
            'delay_res': radar_residuals[:, 0],
            'doppler_res': radar_residuals[:, 1],
            'delay_obs': delay_obs,
            'doppler_obs': doppler_obs,
            'delay_noise': delay_noise,
            'doppler_noise': doppler_noise,
            'delay_comp': delay_obs - radar_residuals[:, 0],
            'doppler_comp': doppler_obs - radar_residuals[:, 1],
        }
        self.all_info['ra_cosdec_res'] = self.all_info['ra_res']*np.cos(dec_obs/3600*np.pi/180)
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
    def scatter_hist(self, x_val, y_val, axis, ax_histx, ax_histy, size, show_logarithmic):
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

    def plot_residuals(self, t_arr, ra_residuals, dec_residuals, ra_cosdec_residuals,
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
        ra_cosdec_residuals : vector
            right ascension residuals multiplied by the cosine of the declination
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
        is_rejected = self.is_rejected
        fig = plt.figure(figsize=(21,6), dpi=150)
        if self.iter_number == 0:
            iter_string = f'Iteration {self.iter_number} (prefit)'
        else:
            iter_string = f'Iteration {self.iter_number}'
        iter_string = title if title is not None else iter_string
        plt.suptitle(iter_string, y=0.95)
        grid_spec = fig.add_gridspec(1, 3, width_ratios=(1,1,1))
        ax1 = fig.add_subplot(grid_spec[0, 0])
        ax1.plot(t_arr, ra_residuals, '.', label='RA', markersize=markersize)
        ax1.plot(t_arr, dec_residuals, '.', label='Dec', markersize=markersize)
        ax1.plot(t_arr[is_rejected], ra_residuals[is_rejected], 'ro',
                    markersize=2*markersize, markerfacecolor='none')
        ax1.plot(t_arr[is_rejected], dec_residuals[is_rejected], 'ro',
                    markersize=2*markersize, markerfacecolor='none')
        ax1.legend()
        ax1.set_xlabel('MJD [UTC]')
        ax1.set_ylabel('Residuals, O-C [arcsec]')
        ax1.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic:
            ax1.set_yscale('log')
        ax2 = grid_spec[0,1].subgridspec(2, 2, width_ratios=(4,1), height_ratios=(1,4),
                                    wspace=0.05, hspace=0.05)
        ax2main = fig.add_subplot(ax2[1,0])
        ax2histx = fig.add_subplot(ax2[0,0], sharex=ax2main)
        ax2histy = fig.add_subplot(ax2[1,1], sharey=ax2main)
        self.scatter_hist(ra_cosdec_residuals, dec_residuals, ax2main, ax2histx, ax2histy,
                            markersize, show_logarithmic)
        ax2main.plot(ra_cosdec_residuals[is_rejected], dec_residuals[is_rejected], 'ro',
                        markersize=2*markersize, markerfacecolor='none')
        ax2main.set_xlabel('RA cos(Dec) Residuals, O-C [arcsec]')
        ax2main.set_ylabel('Dec Residuals, O-C [arcsec]')
        ax2main.grid(True, which='both', axis='both', alpha=0.2, zorder=-100)
        if show_logarithmic:
            ax2main.set_yscale('log')
            ax2main.set_xscale('log')
        ax3 = fig.add_subplot(grid_spec[0, 2])
        ax3.plot(t_arr, delay_residuals, '.', mfc='C2', mec='C2', label='Delay',
                    markersize=radar_scale*markersize)
        ax3.plot(t_arr, doppler_residuals, '.', mfc='C3', mec='C3', label='Doppler',
                    markersize=radar_scale*markersize)
        ax3.legend()
        ax3.set_xlabel('MJD [UTC]')
        ax3.set_ylabel(r'Residuals, O-C [$\mu$s, Hz]')
        ax3.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic:
            ax3.set_yscale('log')
        plt.tight_layout()
        if savefig:
            figname = f'residuals_iter_{self.iter_number}' if figname is None else figname
            plt.savefig(f'{figname}_residuals.pdf', bbox_inches='tight')
        block = not auto_close
        plt.show(block=block)
        if auto_close:
            plt.close(fig)
        return None

    def plot_chi(self, t_arr, ra_chi, dec_chi, delay_chi, doppler_chi,
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
        is_rejected = self.is_rejected
        is_accepted = self.is_accepted
        # optical_idx is where neither the RA nor dec residuals are NaN
        optical_idx = np.where(~np.isnan(self.residuals[:, 0]) &
                                ~np.isnan(self.residuals[:, 1]))[0]
        # radar_idx is where either the RA or dec residuals are NaN
        radar_idx = np.where(np.isnan(self.residuals[:, 0]) |
                                np.isnan(self.residuals[:, 1]))[0]
        # plot chi values
        factor = 3 if plot_chi_squared else 1
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
                f'{np.sum(ra_chi_squared[np.intersect1d(is_accepted, optical_idx)]):.2f}, '
                'Dec='
                f'{np.sum(dec_chi_squared[np.intersect1d(is_accepted, optical_idx)]):.2f}, '
                'Delay='
                f'{np.nansum(delay_chi_squared[np.intersect1d(is_accepted, radar_idx)]):.2f}, '
                'Doppler='
                f'{np.nansum(doppler_chi_squared[np.intersect1d(is_accepted, radar_idx)]):.2f}')
        plt.suptitle(msg, y=0.95)
        if plot_chi_squared:
            plt.subplot(1,2,1)
        plt.plot(t_arr, ra_chi, '.', markersize=markersize, label='RA')
        plt.plot(t_arr, dec_chi, '.', markersize=markersize, label='Dec')
        plt.plot(t_arr[is_rejected], ra_chi[is_rejected], 'ro',
                    markersize=2*markersize, markerfacecolor='none')
        plt.plot(t_arr[is_rejected], dec_chi[is_rejected], 'ro',
                    markersize=2*markersize, markerfacecolor='none')
        plt.plot(t_arr, delay_chi, '.', mfc='C2', mec='C2',
                    markersize=radar_scale*markersize, label='Delay')
        plt.plot(t_arr, doppler_chi, '.', mfc='C3', mec='C3',
                    markersize=radar_scale*markersize, label='Doppler')
        plt.axhline(-sigma_limit, c='khaki', linestyle='--', alpha=1.0,
                        label=fr'$\pm{sigma_limit:.0f}\sigma$')
        plt.axhline(sigma_limit, c='khaki', linestyle='--', alpha=1.0)
        plt.axhline(-2*sigma_limit, c='red', linestyle='--', alpha=0.5,
                        label=fr'$\pm{2*sigma_limit:.0f}\sigma$')
        plt.axhline(2*sigma_limit, c='red', linestyle='--', alpha=0.5)
        plt.legend(ncol=3)
        plt.xlabel('MJD [UTC]')
        plt.ylabel(r'$\chi$, (O-C)/$\sigma$ $[\cdot]$')
        plt.grid(True, which='both', axis='both', alpha=0.2)
        if show_logarithmic:
            plt.yscale('log')
        if plot_chi_squared:
            plt.subplot(1,2,2)
            plt.plot(t_arr, ra_chi_squared, '.', markersize=markersize, label='RA')
            plt.plot(t_arr, dec_chi_squared, '.', markersize=markersize, label='Dec')
            plt.plot(t_arr[is_rejected], ra_chi_squared[is_rejected], 'ro',
                        markersize=2*markersize, markerfacecolor='none')
            plt.plot(t_arr[is_rejected], dec_chi_squared[is_rejected], 'ro',
                        markersize=2*markersize, markerfacecolor='none')
            plt.plot(t_arr, delay_chi_squared, '.', mfc='C2', mec='C2',
                        markersize=radar_scale*markersize, label='Delay')
            plt.plot(t_arr, doppler_chi_squared, '.', mfc='C3', mec='C3',
                        markersize=radar_scale*markersize, label='Doppler')
            plt.legend(ncol=2)
            plt.xlabel('MJD [UTC]')
            plt.ylabel(r'$\chi^2$, (O-C)$^2/\sigma^2$ $[\cdot]$')
            plt.grid(True, which='both', axis='both', alpha=0.2)
            # if show_logarithmic: plt.yscale('log')
            plt.yscale('log')
            plt.tight_layout()
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
        t_arr = self.obs_array[:, 0]
        if show_logarithmic:
            ra_residuals = np.abs(self.all_info['ra_res'])
            dec_residuals = np.abs(self.all_info['dec_res'])
            ra_cosdec_residuals = np.abs(self.all_info['ra_cosdec_res'])
            ra_chi = np.abs(self.all_info['ra_chi'])
            dec_chi = np.abs(self.all_info['dec_chi'])
        else:
            ra_residuals = self.all_info['ra_res']
            dec_residuals = self.all_info['dec_res']
            ra_cosdec_residuals = self.all_info['ra_cosdec_res']
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
        self.plot_residuals(t_arr, ra_residuals, dec_residuals, ra_cosdec_residuals,
                            delay_residuals, doppler_residuals, radar_scale, markersize,
                            show_logarithmic, title, savefig, figname, auto_close)
        self.plot_chi(t_arr, ra_chi, dec_chi, delay_chi, doppler_chi,
                        ra_chi_squared, dec_chi_squared, delay_chi_squared, doppler_chi_squared,
                        plot_chi_squared, sigma_limit, radar_scale, markersize,
                        show_logarithmic, title, savefig, figname, auto_close)
        return None

class FitSimulation:
    """
    Class to perform an orbit fit simulation.
    """
    def __init__(self, x_init, cov_init=None, obs_array_optical=None, observer_codes_optical=None,
                    obs_array_radar=None, observer_codes_radar=None, n_iter_max=10,
                    de_kernel=441, de_kernel_path='', radius=0.0, nongrav_info=None, events=None):
        """
        Constructor of the FitSimulation class.

        Parameters
        ----------
        x_init : dict
            Initial guess of the state vector.
        cov_init : array, optional
            Initial covariance matrix, by default None
        obs_array_optical : array, optional
            Optical observation data for the orbit fit, by default None
        observer_codes_optical : tuple, optional
            Observer locations for each observation in obs_array_optical, by default None
        obs_array_radar : array, optional
            Radar observation data for the orbit fit, by default None
        observer_codes_radar : tuple, optional
            Observer locations for each observation in obs_array_radar, by default None
        n_iter_max : int, optional
            Number of maximum iterations to correct the orbit estimate, by default 10
        de_kernel : int, optional
            SPICE kernel version, by default 441
        de_kernel_path : str, optional
            Path to the SPICE kernel, by default ''
        radius : float, optional
            Radius of the body, by default 0.0
        nongrav_info : dict, optional
            Dictionary containing the non-gravitational parameters, by default None
        events : list, optional
            List of lists containing impulsive maneuvers, by default None

        Returns
        -------
        None : NoneType
            None
        """
        self.x_nom = None
        self.covariance = None
        self.obs_array = None
        self.observer_codes = None
        self.check_initial_solution(x_init, cov_init)
        self.check_input_observation_arrays(obs_array_optical, observer_codes_optical,
                                            obs_array_radar, observer_codes_radar)
        self.sigmas = None
        self.obs_cov = None
        self.obs_weight = None
        self.assemble_observation_arrays()
        self.n_iter = 0
        self.n_iter_max = n_iter_max
        self.iters = []
        self.de_kernel = de_kernel
        self.de_kernel_path = de_kernel_path
        self.analytic_partials = False
        self.fixed_propsim_params = {'a1': 0.0, 'a2': 0.0, 'a3': 0.0,
                                        'alpha': 1.0, 'k': 0.0, 'm': 2.0, 'n': 0.0,
                                        'r0_au': 1.0, 'radius': radius, 'mass': 0.0}
        for key in nongrav_info:
            self.fixed_propsim_params[key] = nongrav_info[key]
        if events is not None:
            self.fixed_propsim_params['events'] = events
            self.fit_events = True
        else:
            self.fixed_propsim_params['events'] = []
            self.fit_events = False
        self.residual_chi_squared = None
        self.converged = False
        return None

    def check_initial_solution(self, x_init, cov_init):
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
        for key in ["a1", "a2", "a3"]:
            if key in x_init and x_init[key] != 0.0:
                self.fit_nongrav = True
        self.t_sol = x_init['t']
        self.x_init = x_init
        self.x_nom = {key: x_init[key] for key in x_init if key != 't'} # remove t for self.x_nom
        self.n_fit = len(self.x_nom)
        if cov_init.shape != (self.n_fit, self.n_fit):
            msg = ("Covariance matrix must be the same size "
                    "as the number of fitted parameters.")
            raise ValueError(msg)
        self.covariance_init = cov_init
        self.covariance = cov_init
        return None

    def check_input_observation_arrays(self, obs_array_optical, observer_codes_optical,
                                        obs_array_radar, observer_codes_radar):
        """
        Check the input observation arrays provided by the user and
        make sure they are valid.

        Parameters
        ----------
        obs_array_optical : array
            Optical observation data for the orbit fit.
        observer_codes_optical : tuple
            Observer locations for each observation in obs_array_optical.
        obs_array_radar : array
            Radar observation data for the orbit fit.
        observer_codes_radar : tuple
            Observer locations for each observation in obs_array_radar.

        Returns
        -------
        None : NoneType
            None

        Raises
        ------
        ValueError
            If no observation arrays are provided.
        ValueError
            If optical observations are provided but no observer codes.
        ValueError
            If observer codes are provided but no optical observations.
        ValueError
            If the length of the optical observation array and observer
            code array are not the same.
        ValueError
            If radar observations are provided but no observer codes.
        ValueError
            If observer codes are provided but no radar observations.
        ValueError
            If the length of the radar observation array and observer
            code array are not the same.
        """
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
                msg = ("Optical observation array and observer "
                        "code array must be the same length.")
                raise ValueError(msg)
            else:
                self.fit_optical = True
        if obs_array_radar is not None and observer_codes_radar is None:
            raise ValueError("Must provide observer codes for radar observations.")
        if obs_array_radar is None and observer_codes_radar is not None:
            raise ValueError("Must provide radar observations for observer codes.")
        if obs_array_radar is not None and observer_codes_radar is not None:
            if len(obs_array_radar) != len(observer_codes_radar):
                msg = ("Radar observation array and observer "
                        "code array must be the same length.")
                raise ValueError(msg)
            else:
                self.fit_radar = True
        self.obs_array_optical = obs_array_optical
        self.observer_codes_optical = observer_codes_optical
        self.obs_array_radar = obs_array_radar
        self.observer_codes_radar = observer_codes_radar
        return None

    def flatten_and_clean(self, arr):
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

    def assemble_observation_arrays(self):
        """
        Assemble the optical and/or radar observation arrays into a
        single array for orbit fitting.

        Returns
        -------
        None : NoneType
            None
        """
        if self.fit_optical and self.fit_radar:
            self.merge_observation_arrays()
        elif self.fit_optical:
            self.obs_array, sort_idx = self.sort_array_by_another(self.obs_array_optical,
                                                                    self.obs_array_optical[:, 0])
            self.observer_codes = tuple(np.array(self.observer_codes_optical,
                                                    dtype=tuple)[sort_idx])
        elif self.fit_radar:
            self.obs_array, sort_idx = self.sort_array_by_another(self.obs_array_radar,
                                                                    self.obs_array_radar[:, 0])
            self.observer_codes = tuple(np.array(self.observer_codes_radar,
                                                    dtype=tuple)[sort_idx])
        # number of observations is the number of non-nan values
        # in the second and third columns of the observation array
        self.n_obs = np.count_nonzero(~np.isnan(self.obs_array[:, 1:3]))
        self.rejection_flag = [False]*len(self.obs_array)
        self.create_obs_weight()
        self.past_obs_idx = np.where(self.obs_array[:, 0] < self.t_sol)[0]
        self.past_obs_exist = len(self.past_obs_idx) > 0
        self.future_obs_idx = np.where(self.obs_array[:, 0] >= self.t_sol)[0]
        self.future_obs_exist = len(self.future_obs_idx) > 0
        return None

    def create_obs_weight(self):
        """
        Assembles the weight matrix for the orbit fit based
        on observation uncertainties and correlations.

        Returns
        -------
        None : NoneType
            None
        """
        self.sigmas = self.obs_array[:, 3:5]
        sigma_corr = self.obs_array[:, 5]
        obs_cov = np.zeros((self.sigmas.size, self.sigmas.size))
        idx_to_remove = []
        for i in range(self.sigmas.shape[0]):
            sig_1, sig_2 = self.sigmas[i]
            sig_1_nan = np.isnan(sig_1)
            sig_2_nan = np.isnan(sig_2)
            sig_corr = sigma_corr[i]
            sig_corr_nan = np.isnan(sig_corr)
            off_diag = 0.0 if (sig_1_nan or sig_2_nan or sig_corr_nan) else sig_corr*sig_1*sig_2
            sub_cov = np.array([[sig_1**2, off_diag],
                                    [off_diag, sig_2**2]])
            obs_cov[2*i:2*i+2, 2*i:2*i+2] = sub_cov
            if sig_1_nan:
                idx_to_remove.append(2*i)
            if sig_2_nan:
                idx_to_remove.append(2*i+1)
            if sig_1_nan and sig_2_nan:
                raise ValueError("Both sigmas cannot be NaN.")
        # remove all rows and columns with NaN values from the obs_cov
        obs_cov = np.delete(obs_cov, idx_to_remove, axis=0)
        obs_cov = np.delete(obs_cov, idx_to_remove, axis=1)
        self.obs_cov = obs_cov
        self.obs_weight = np.linalg.inv(self.obs_cov)
        return None

    def sort_array_by_another(self, array, sort_by):
        """
        Sort an array by another array.

        Parameters
        ----------
        array : array
            Array to sort.
        sort_by : array
            Array to sort by.

        Returns
        -------
        array : array
            Sorted array.
        sort_idx : vector
            Indices of the sorted array.
        """
        sort_idx = np.argsort(sort_by)
        return array[sort_idx], sort_idx

    def merge_observation_arrays(self):
        """
        Merge the optical and radar observation arrays into a single
        array for orbit fitting.

        Returns
        -------
        None : NoneType
            None
        """
        # merge optical and radar arrays
        obs_array = np.vstack((self.obs_array_optical, self.obs_array_radar))
        observer_codes = self.observer_codes_optical + self.observer_codes_radar
        # sort by time
        self.obs_array, sort_idx = self.sort_array_by_another(obs_array, obs_array[:, 0])
        self.observer_codes = tuple(np.array(observer_codes, dtype=tuple)[sort_idx])
        return None

    def get_prop_sim_past(self, name, t_eval_utc, eval_apparent_state,
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
        prop_sim_past : prop.propSimulation object
            propSim object for the past observations.
        """
        t_eval_past = self.obs_array[self.past_obs_idx, 0]
        tf_past = np.min(t_eval_past)
        prop_sim_past = prop.propSimulation(name, self.t_sol, self.de_kernel, self.de_kernel_path)
        prop_sim_past.tEvalMargin = 1.0
        prop_sim_past.set_integration_parameters(tf_past, t_eval_past, t_eval_utc,
                                                    eval_apparent_state, converged_light_time,
                                                    observer_info)
        return prop_sim_past

    def get_prop_sim_future(self, name, t_eval_utc, eval_apparent_state,
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
        prop_sim_future : prop.propSimulation object
            propSim object for the future observations.
        """
        t_eval_future = self.obs_array[self.future_obs_idx, 0]
        tf_future = np.max(t_eval_future)
        prop_sim_future = prop.propSimulation(name, self.t_sol, self.de_kernel, self.de_kernel_path)
        prop_sim_future.tEvalMargin = 1.0
        prop_sim_future.set_integration_parameters(tf_future, t_eval_future, t_eval_utc,
                                                    eval_apparent_state, converged_light_time,
                                                    observer_info)
        return prop_sim_future

    def get_prop_sims(self, name):
        """
        Get propSim objects for the past and future observations.

        Parameters
        ----------
        name : str
            Name of the propSim objects.

        Returns
        -------
        prop_sim_past : prop.propSimulation object
            propSim object for the past observations.
        prop_sim_future : prop.propSimulation object
            propSim object for the future observations.
        """
        t_eval_utc = True
        eval_apparent_state = True
        converged_light_time = True
        observer_info = np.array(get_observer_info(self.observer_codes), dtype=tuple)
        observer_info_past = tuple(observer_info[self.past_obs_idx])
        observer_info_future = tuple(observer_info[self.future_obs_idx])
        prop_sim_past = None
        prop_sim_future = None
        if self.past_obs_exist:
            prop_sim_past = self.get_prop_sim_past(f"{name}_past", t_eval_utc,
                                                    eval_apparent_state, converged_light_time,
                                                    observer_info_past)
        if self.future_obs_exist:
            prop_sim_future = self.get_prop_sim_future(f"{name}_future", t_eval_utc,
                                                        eval_apparent_state, converged_light_time,
                                                        observer_info_future)
        return prop_sim_past, prop_sim_future

    def x_dict_to_state(self, x_dict):
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
                                    omega*np.pi/180.0, arg_peri*np.pi/180.0, inc*np.pi/180.0]
            state = cometary_elements
        else:
            raise ValueError("fit_cartesian or fit_cometary must be True")
        return state

    def x_dict_to_nongrav_params(self, x_dict):
        """
        Convert a dictionary of nominal state to non-gravitational parameters.

        Parameters
        ----------
        x_dict : dict
            Dictionary of the nominal state.

        Returns
        -------
        nongrav_params : prop.NongravParamaters object
            Non-gravitational parameters for the fitted body.
        """
        nongrav_params = prop.NongravParamaters()
        a1_val = x_dict['a1'] if 'a1' in x_dict.keys() else self.fixed_propsim_params['a1']
        a2_val = x_dict['a2'] if 'a2' in x_dict.keys() else self.fixed_propsim_params['a2']
        a3_val = x_dict['a3'] if 'a3' in x_dict.keys() else self.fixed_propsim_params['a3']
        nongrav_params.a1 = a1_val
        nongrav_params.a2 = a2_val
        nongrav_params.a3 = a3_val
        nongrav_params.alpha = self.fixed_propsim_params['alpha']
        nongrav_params.k = self.fixed_propsim_params['k']
        nongrav_params.m = self.fixed_propsim_params['m']
        nongrav_params.n = self.fixed_propsim_params['n']
        nongrav_params.r0_au = self.fixed_propsim_params['r0_au']
        return nongrav_params

    def x_dict_to_events(self, x_dict):
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
            event = self.fixed_propsim_params['events'][i]
            event[1] = x_dict[f"dvx{i}"] if f"dvx{i}" in x_dict.keys() else event[1]
            event[2] = x_dict[f"dvy{i}"] if f"dvy{i}" in x_dict.keys() else event[2]
            event[3] = x_dict[f"dvz{i}"] if f"dvz{i}" in x_dict.keys() else event[3]
            event[4] = x_dict[f"mult{i}"] if f"mult{i}" in x_dict.keys() else event[4]
            events.append(tuple(event))
        return events

    def check_and_add_events(self, prop_sim_past, prop_sim_future, integ_body, events):
        """
        Check if events are in the past or future and add them to the appropriate prop_sim.

        Parameters
        ----------
        prop_sim_past : prop.propSimulation object
            propSimulation object for the past.
        prop_sim_future : prop.propSimulation object
            propSimulation object for the future.
        integ_body : prop.IntegBody object
            IntegBody object for the body being fitted.
        events : list
            List of events for the fitted body.

        Returns
        -------
        prop_sim_past : prop.propSimulation object
            propSimulation object for the past.
        prop_sim_future : prop.propSimulation object
            propSimulation object for the future.
        """
        for event in events:
            t_event = event[0]
            dvx = event[1]
            dvy = event[2]
            dvz = event[3]
            multiplier = event[4]
            if t_event < self.t_sol:
                prop_sim_past.add_event(integ_body, t_event, [dvx, dvy, dvz], multiplier)
            else:
                prop_sim_future.add_event(integ_body, t_event, [dvx, dvy, dvz], multiplier)
        return prop_sim_past, prop_sim_future

    def get_perturbed_state(self, key):
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
        ng_params_plus : prop.NongravParamaters object
            Perturbed non-gravitational parameters in the positive direction.
        events_plus : list
            Perturbed events in the positive direction.
        state_minus : list
            Perturbed state in the negative direction.
        ng_params_minus : prop.NongravParamaters object
            Perturbed non-gravitational parameters in the negative direction.
        events_minus : list
            Perturbed events in the negative direction.
        fd_delta : float
            Finite difference perturbation.
        """
        if self.fit_cartesian:
            if key in ['x', 'y', 'z']:
                fd_pert = 1e-8
            elif key in ['vx', 'vy', 'vz']:
                fd_pert = 1e-10
        elif self.fit_cometary:
            fd_pert = 1e-8
            # if key in ['peri_dist']:
            #     fd_pert = 5e-5
            # elif key in ['ecc']:
            #     fd_pert = 5e-4
            # elif key in ['time_peri']:
            #     fd_pert = 1e-5
            # elif key in ['omega']:
            #     fd_pert = 1e-6
            # elif key in ['w']:
            #     fd_pert = 1e-6
            # elif key in ['i']:
            #     fd_pert = 1e-6
        if key in ['a1', 'a2', 'a3']:
            fd_pert = 1e-3
        if key[:4] == 'mult':
            fd_pert = 1e0
        if key[:3] in ['dvx', 'dvy', 'dvz']:
            fd_pert = 1e-4
        x_plus = self.x_nom.copy()
        x_minus = self.x_nom.copy()
        # fd_pert = finite difference perturbation to nominal state for calculating derivatives
        fd_delta = self.x_nom[key]*fd_pert
        x_plus[key] = self.x_nom[key]+fd_delta
        state_plus = self.x_dict_to_state(x_plus)
        ng_params_plus = self.x_dict_to_nongrav_params(x_plus)
        events_plus = self.x_dict_to_events(x_plus)
        x_minus[key] = self.x_nom[key]-fd_delta
        state_minus = self.x_dict_to_state(x_minus)
        ng_params_minus = self.x_dict_to_nongrav_params(x_minus)
        events_minus = self.x_dict_to_events(x_minus)
        return (state_plus, ng_params_plus, events_plus, state_minus,
                    ng_params_minus, events_minus, fd_delta)

    def get_perturbation_info(self):
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
            pert_result = self.get_perturbed_state(key)
            perturbation_info.append(tuple(pert_result))
        return perturbation_info

    def assemble_and_propagate_bodies(self, perturbation_info):
        """
        Assemble and propagate the body for the orbit fit.

        Parameters
        ----------
        perturbation_info : list
            List of tuples containing the perturbation information
            for each nominal state parameter.

        Returns
        -------
        prop_sim_past : prop.propSimulation object
            propagated propSimulation object for the past.
        prop_sim_future : prop.propSimulation object
            propagated propSimulation object for the future.
        """
        # sourcery skip: low-code-quality
        # get propagated states
        prop_sim_past, prop_sim_future = self.get_prop_sims("orbit_fit_sim")
        # create nominal integ_body object
        consts = prop.Constants()
        state_nom = self.x_dict_to_state(self.x_nom)
        ng_params_nom = self.x_dict_to_nongrav_params(self.x_nom)
        events_nom = self.x_dict_to_events(self.x_nom)
        cov_nom = self.covariance
        if self.fit_cartesian:
            integ_body_nom = prop.IntegBody("integ_body_nom", self.t_sol,
                                            self.fixed_propsim_params['mass'],
                                            self.fixed_propsim_params['radius'],
                                            state_nom[:3], state_nom[3:6],
                                            cov_nom, ng_params_nom, consts)
        elif self.fit_cometary:
            integ_body_nom = prop.IntegBody(self.de_kernel_path, "integ_body_nom",
                                            self.t_sol, self.fixed_propsim_params['mass'],
                                            self.fixed_propsim_params['radius'],
                                            state_nom, cov_nom, ng_params_nom, consts)
        # add the nominal integ_body for the residuals
        if self.past_obs_exist:
            prop_sim_past.add_integ_body(integ_body_nom)
        if self.future_obs_exist:
            prop_sim_future.add_integ_body(integ_body_nom)
        prop_sim_past, prop_sim_future = self.check_and_add_events(prop_sim_past, prop_sim_future,
                                                                    integ_body_nom, events_nom)
        # add the perturbed IntegBodies for numerical derivatives
        if not self.analytic_partials and perturbation_info is not None:
            for i in range(self.n_fit):
                key = list(self.x_nom.keys())[i]
                (state_plus, ng_params_plus, events_plus, state_minus,
                        ng_params_minus, events_minus, _) = perturbation_info[i]
                if self.fit_cartesian:
                    integ_body_plus = prop.IntegBody(f"integBody_pert_{key}_plus", self.t_sol,
                                                        self.fixed_propsim_params['mass'],
                                                        self.fixed_propsim_params['radius'],
                                                        state_plus[:3], state_plus[3:6],
                                                        cov_nom, ng_params_plus, consts)
                    integ_body_minus = prop.IntegBody(f"integBody_pert_{key}_minus", self.t_sol,
                                                        self.fixed_propsim_params['mass'],
                                                        self.fixed_propsim_params['radius'],
                                                        state_minus[:3], state_minus[3:6],
                                                        cov_nom, ng_params_minus, consts)
                elif self.fit_cometary:
                    integ_body_plus = prop.IntegBody(self.de_kernel_path,
                                                        f"integBody_pert_{key}_plus", self.t_sol,
                                                        self.fixed_propsim_params['mass'],
                                                        self.fixed_propsim_params['radius'],
                                                        state_plus, cov_nom,
                                                        ng_params_plus, consts)
                    integ_body_minus = prop.IntegBody(self.de_kernel_path,
                                                        f"integBody_pert_{key}_minus", self.t_sol,
                                                        self.fixed_propsim_params['mass'],
                                                        self.fixed_propsim_params['radius'],
                                                        state_minus, cov_nom,
                                                        ng_params_minus, consts)
                if self.past_obs_exist:
                    prop_sim_past.add_integ_body(integ_body_plus)
                    prop_sim_past.add_integ_body(integ_body_minus)
                if self.future_obs_exist:
                    prop_sim_future.add_integ_body(integ_body_plus)
                    prop_sim_future.add_integ_body(integ_body_minus)
                prop_sim_past, prop_sim_future = self.check_and_add_events(prop_sim_past,
                                                                            prop_sim_future,
                                                                            integ_body_plus,
                                                                            events_plus)
                prop_sim_past, prop_sim_future = self.check_and_add_events(prop_sim_past,
                                                                            prop_sim_future,
                                                                            integ_body_minus,
                                                                            events_minus)
        if self.past_obs_exist:
            prop_sim_past.preprocess()
            prop_sim_past.integrate()
        if self.future_obs_exist:
            prop_sim_future.preprocess()
            prop_sim_future.integrate()
        return prop_sim_past, prop_sim_future

    def get_computed_obs(self, prop_sim_past, prop_sim_future, integ_body_idx):
        """
        Computes the optical and radar observations from the propagated states.

        Parameters
        ----------
        prop_sim_past : prop.propSimulation object
            The propagated propSimulation object for the past.
        prop_sim_future : prop.propSimulation object
            The propagated propSimulation object for the future.
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
            apparent_states_past = np.array(prop_sim_past.xIntegEval)
            apparent_states_future = np.array(prop_sim_future.xIntegEval)
            apparent_states = np.vstack((apparent_states_past, apparent_states_future))
            radar_observations_past = np.array(prop_sim_past.radarObsEval)
            radar_observations_future = np.array(prop_sim_future.radarObsEval)
            radar_observations = np.vstack((radar_observations_past, radar_observations_future))
        elif self.past_obs_exist:
            apparent_states_past = np.array(prop_sim_past.xIntegEval)
            apparent_states = apparent_states_past
            radar_observations_past = np.array(prop_sim_past.radarObsEval)
            radar_observations = radar_observations_past
        elif self.future_obs_exist:
            apparent_states_future = np.array(prop_sim_future.xIntegEval)
            apparent_states = apparent_states_future
            radar_observations_future = np.array(prop_sim_future.radarObsEval)
            radar_observations = radar_observations_future
        integ_body_start_col = 6*integ_body_idx
        integ_body_end_col = 6*integ_body_idx+6
        apparent_states = apparent_states[:,integ_body_start_col:integ_body_end_col]
        radar_observations = radar_observations[:,integ_body_idx]
        measured_obs = self.obs_array[:, 1:3]
        computed_obs = np.nan*np.ones_like(measured_obs)
        observer_info = get_observer_info(self.observer_codes)
        for i in range(len(self.obs_array)):
            obs_info_len = len(observer_info[i])
            if obs_info_len in {4, 7}:
                computed_obs[i, :] = get_radec(apparent_states[i])
            elif obs_info_len == 9: # delay measurement
                computed_obs[i, 0] = radar_observations[i]
            elif obs_info_len == 10: # dopper measurement
                computed_obs[i, 1] = radar_observations[i]
            else:
                raise ValueError("Observer info length not recognized.")
        return computed_obs

    def get_analytic_partials(self, prop_sim_past, prop_sim_future):
        """
        Computes the analytic partials of the observations with respect to the
        initial nominal state.

        Parameters
        ----------
        prop_sim_past : prop.propSimulation object
            The propagated propSimulation object for the past.
        prop_sim_future : prop.propSimulation object
            The propagated propSimulation object for the future.

        Raises
        ------
        NotImplementedError
            Because analytic partials are not yet implemented.
        """
        raise NotImplementedError("Analytic partials not yet implemented.",
                                        "Please use numeric partials.")

    def get_numeric_partials(self, prop_sim_past, prop_sim_future, perturbation_info):
        """
        Computes the numeric partials of the observations with respect to the
        initial nominal state.

        Parameters
        ----------
        prop_sim_past : prop.propSimulation object
            The propagated propSimulation object for the past.
        prop_sim_future : prop.propSimulation object
            The propagated propSimulation object for the future.
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
            computed_obs_plus = self.get_computed_obs(prop_sim_past, prop_sim_future,
                                                        integ_body_idx=2*i+1)
            computed_obs_minus = self.get_computed_obs(prop_sim_past, prop_sim_future,
                                                        integ_body_idx=2*i+2)
            computed_obs_plus = self.flatten_and_clean(computed_obs_plus)
            computed_obs_minus = self.flatten_and_clean(computed_obs_minus)
            # get partials
            partials[:, i] = (computed_obs_plus - computed_obs_minus)/(2*fd_delta)
        return partials

    def get_partials(self, prop_sim_past, prop_sim_future, perturbation_info):
        """
        Computes the partials of the observations with respect to the
        initial nominal state.

        Parameters
        ----------
        prop_sim_past : prop.propSimulation object
            The propagated propSimulation object for the past.
        prop_sim_future : prop.propSimulation object
            The propagated propSimulation object for the future.
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
            return self.get_analytic_partials(prop_sim_past, prop_sim_future)
        else:
            return self.get_numeric_partials(prop_sim_past, prop_sim_future, perturbation_info)

    def get_residuals_and_partials(self):
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
        perturbation_info = None if self.analytic_partials else self.get_perturbation_info()
        prop_sim_past, prop_sim_future = self.assemble_and_propagate_bodies(perturbation_info)
        # get residuals
        computed_obs = self.get_computed_obs(prop_sim_past, prop_sim_future, integ_body_idx=0)
        residuals = self.obs_array[:, 1:3] - computed_obs
        # get partials
        partials = self.get_partials(prop_sim_past, prop_sim_future, perturbation_info)
        return residuals, partials

    def reject_outliers(self, partials, weights, residuals):
        """
        Outlier rejection algorithm for the residuals.

        Parameters
        ----------
        partials : array
            The partials of the observations with respect to the
            initial nominal state.
        weights : array
            The weights of the observations.
        residuals : array
            The residuals of the observations

        Returns
        -------
        partials : array
            The partials of the observations with respect to the
            initial nominal state after outlier rejection.
        weights : array
            The weights of the observations after outlier rejection.
        residuals : array
            The residuals of the observations after outlier rejection.
        residual_chi_squared : array
            The chi squared values of the residuals after outlier rejection.

        Raises
        ------
        ValueError
            If the observer information is not well-defined.
        """
        chi_reject = 3.0
        chi_recover = 2.8
        full_cov = np.linalg.inv(partials.T @ self.obs_weight @ partials)
        observer_info = get_observer_info(self.observer_codes)
        j = 0
        residual_chi_squared = np.zeros(len(self.obs_array))
        rejected_indices = []
        resid_cov_full = self.obs_cov - partials @ full_cov @ partials.T
        for i in range(len(self.obs_array)):
            obs_info_len = len(observer_info[i])
            if obs_info_len in {4, 7}:
                size = 2
            elif obs_info_len in {9, 10}:
                size = 1
                j += size
                continue
            else:
                raise ValueError("Observer info length not recognized.")
            # calculate chi-squared for each residual
            resid = residuals[j:j+size].reshape((1, size))
            # resid_cov = covs - partials @ full_cov @ partials.T
            resid_cov = resid_cov_full[j:j+size, j:j+size]
            residual_chi_squared[i] = resid @ np.linalg.inv(resid_cov) @ resid.T
            # outlier rejection, only reject RA/Dec measurements
            if residual_chi_squared[i] > chi_reject**2 and size == 2:
                self.rejection_flag[i] = True
            if self.rejection_flag[i] and residual_chi_squared[i] < chi_recover**2:
                self.rejection_flag[i] = False
                print("Recovered observation at time", self.obs_array[i, 0])
            # # sigma clipping
            # if (abs(resid[0, 0])/self.obs_array[i,3] >= chi_reject
            #         or abs(resid[0, 1])/self.obs_array[i,4] >= chi_reject):
            #     self.rejection_flag[i] = True
            # if (self.rejection_flag[i] and (abs(resid[0, 0])/self.obs_array[i,3] < chi_recover
            #         and abs(resid[0, 1])/self.obs_array[i,4] < chi_recover)):
            #     self.rejection_flag[i] = False
            #     print("Recovered observation at time", self.obs_array[i, 0])
            if self.rejection_flag[i]:
                rejected_indices.extend((j, j+1))
            j += size
        # remove rejected residuals from partials, weights, and residuals
        partials = np.delete(partials, rejected_indices, axis=0)
        weights = np.delete(weights, rejected_indices, axis=0)
        weights = np.delete(weights, rejected_indices, axis=1)
        residuals = np.delete(residuals, rejected_indices, axis=0)
        # print(rejected_indices)
        num_rejected = len(rejected_indices)//2
        num_total = len(self.obs_array)
        if num_rejected/num_total > 0.25:
            print("WARNING: More than 25% of observations rejected. Consider changing chi_reject",
                    "and chi_recover values, or turning off outlier rejection altogether.")
        self.residual_chi_squared = residual_chi_squared
        return partials, weights, residuals, residual_chi_squared

    def add_iteration(self, iter_number, residuals, obs_weight_rej):
        """
        Adds an iteration to the list of iterations in the FitSimulation object.

        Parameters
        ----------
        iter_number : int
            The iteration number.
        residuals : array
            The residuals of the observations.
        obs_weight_rej : array
            Observation weight matrix after outlier rejection.

        Returns
        -------
        None : NoneType
            None
        """
        self.iters.append(IterationParams(iter_number, self.x_nom, self.covariance, residuals,
                                            self.obs_array, obs_weight_rej, self.observer_codes,
                                            self.rejection_flag))
        return None

    def check_convergence(self):
        """
        Checks if the orbit fit has converged.

        Returns
        -------
        None : NoneType
            None
        """
        if self.n_iter > 1:
            del_rms_convergence = 1e-3
            curr_rms = self.iters[-1].weighted_rms
            prev_rms = self.iters[-2].weighted_rms
            del_rms = abs(prev_rms - curr_rms)/prev_rms
            if del_rms < del_rms_convergence:
                self.converged = True
        return None

    def filter_lsq(self, verbose=True):
        """
        Performs a least-squares fit on the observations.

        Parameters
        ----------
        verbose : bool, optional
            FLag for printing the iteration information while fitting, by default True.

        Returns
        -------
        None : NoneType
            None
        """
        start_rejecting = False
        spice.furnsh(self.de_kernel_path)
        if verbose:
            print("Iteration\t\tUnweighted RMS\t\tWeighted RMS",
                    "\t\tChi-squared\t\tReduced Chi-squared")
        for i in range(self.n_iter_max):
            self.n_iter = i+1
            # get residuals and partials
            residuals, partials = self.get_residuals_and_partials()
            weights = self.obs_weight
            if i == 0:
                # add prefit iteration
                prefit_residuals = residuals.copy()
                self.add_iteration(0, prefit_residuals, weights)
            clean_residuals = self.flatten_and_clean(residuals)
            # reject outliers here
            if start_rejecting:
                a_rej, w_rej, b_rej, _ = self.reject_outliers(partials, weights, clean_residuals)
            else:
                a_rej, w_rej, b_rej = partials, weights, clean_residuals
            # get initial guess
            curr_state = np.array(list(self.x_nom.values()))
            # get state correction
            atwa = a_rej.T @ w_rej @ a_rej
            cov = np.linalg.inv(atwa)
            atwb = a_rej.T @ w_rej @ b_rej
            delta_x = cov @ atwb
            # delta_x = np.linalg.solve(atwa, atwb)
            # get new state
            next_state = curr_state + delta_x
            self.x_nom = dict(zip(self.x_nom.keys(), next_state))
            # get new covariance
            self.covariance = cov
            # add iteration
            self.add_iteration(i+1, residuals, w_rej)
            if verbose:
                print(f"{self.iters[-1].iter_number}\t\t\t",
                        f"{self.iters[-1].unweighted_rms:.3f}\t\t\t",
                        f"{self.iters[-1].weighted_rms:.3f}\t\t\t",
                        f"{self.iters[-1].chi_squared:.3f}\t\t\t",
                        f"{self.iters[-1].reduced_chi_squared:.3f}")
            self.check_convergence()
            if self.converged:
                if start_rejecting:
                    print("Converged after rejecting outliers.")
                    break
                else:
                    print("Converged without rejecting outliers. Starting outlier rejection now.")
                    start_rejecting = True
                    self.converged = False
        if self.n_iter == self.n_iter_max and not self.converged:
            print("WARNING: Maximum number of iterations reached without converging.")
        spice.unload(self.de_kernel_path)
        return None

    def print_summary(self, iter_idx=-1):
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
        data = self.iters[iter_idx]
        print("Summary of the orbit fit calculations at iteration",
                    f"{data.iter_number} (of {self.n_iter}):")
        print("=======================================================")
        print(f"RMS unweighted: {data.unweighted_rms}")
        print(f"RMS weighted: {data.weighted_rms}")
        print(f"chi-squared: {data.chi_squared}")
        print(f"reduced chi-squared: {data.reduced_chi_squared}")
        print(f"square root of reduced chi-squared: {np.sqrt(data.reduced_chi_squared)}")
        print("=======================================================")
        print(f"t: MJD {self.t_sol} TDB")
        print("Fitted Variable\t\tInitial Value\t\t\tUncertainty\t\t\tFitted Value",
                    "\t\t\tUncertainty\t\t\tChange\t\t\t\tChange (sigma)")
        init_variance = np.sqrt(np.diag(self.covariance_init))
        final_variance = np.sqrt(np.diag(self.covariance))
        init_sol = self.iters[0].x_nom
        final_sol = data.x_nom
        with np.errstate(divide='ignore'):
            for i, key in enumerate(init_sol.keys()):
                print(f"{key}\t\t\t{init_sol[key]:.11e}\t\t{init_variance[i]:.11e}",
                        f"\t\t{final_sol[key]:.11e}\t\t{final_variance[i]:.11e}",
                        f"\t\t{final_sol[key]-init_sol[key]:+.11e}"
                        f"\t\t{(final_sol[key]-init_sol[key])/init_variance[i]:+.3f}")
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
        ticks = np.arange(1, self.n_iter+1, 1)
        start_idx = 1
        iters_for_plot = self.iters[start_idx:]
        final_iter = iters_for_plot[-1]
        plt.figure(figsize=(21,10), dpi=150)
        plt.subplot(2, 2, 1)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot],
                        [iteration.unweighted_rms for iteration in iters_for_plot],
                        label=f"Final Unweighted RMS={final_iter.unweighted_rms:.3e}")
        plt.xticks(ticks)
        plt.xlabel("Iteration #")
        plt.ylabel("Unweighted RMS")
        plt.legend()
        plt.subplot(2, 2, 2)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot],
                        [iteration.weighted_rms for iteration in iters_for_plot],
                        label=f"Final Weighted RMS={final_iter.weighted_rms:.3e}")
        plt.xticks(ticks)
        plt.xlabel("Iteration #")
        plt.ylabel("Weighted RMS")
        plt.legend()
        plt.subplot(2, 2, 3)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot],
                        [iteration.chi_squared for iteration in iters_for_plot],
                        label=fr"Final $\chi^2$={final_iter.chi_squared:.3e}")
        plt.xticks(ticks)
        plt.xlabel("Iteration #")
        plt.ylabel(r"$\chi^2$")
        plt.legend()
        plt.subplot(2, 2, 4)
        plt.semilogy([iteration.iter_number for iteration in iters_for_plot],
                        [iteration.reduced_chi_squared for iteration in iters_for_plot],
                        label=fr"Final Reduced $\chi^2$={final_iter.reduced_chi_squared:.3e}")
        plt.xticks(ticks)
        plt.xlabel("Iteration #")
        plt.ylabel(r"Reduced $\chi^2$")
        plt.legend()
        block = not auto_close
        plt.show(block=block)
        if auto_close:
            plt.close()
        return None

def _generate_simulated_obs(ref_sol, ref_cov, ref_ng_info, events, modified_obs_arrays,
                            simulated_optical_obs_idx, optical_obs_types,
                            simulated_radar_obs_idx, radar_obs_types, de_kernel,
                            de_kernel_path, noise):
    """
    Generates simulated observations for a given reference solution.

    Parameters
    ----------
    ref_sol : dict
        Dictionary containing the reference solution.
    ref_cov : array
        Covariance matrix for the reference solution.
    ref_ng_info : dict
        Dictionary containing the reference non-gravitational acceleration information.
    events : list
        List of lists containing the event information.
    modified_obs_arrays : tuple
        Tuple containing the modified optical and radar observation arrays and observer codes.
    simulated_optical_obs_idx : vector
        Vector containing the indices of the simulated optical observations.
    optical_obs_types : list
        List of optical observation types.
    simulated_radar_obs_idx : vector
        Vector containing the indices of the simulated radar observations.
    radar_obs_types : list
        List of radar observation types.
    de_kernel : int
        SPICE kernel version.
    de_kernel_path : str
        Path to SPICE kernel.
    noise : bool
        Flag to add noise to the simulated observations.

    Returns
    -------
    obs_array_optical : array
        Simulated optical observation data for the given reference solution.
    observer_codes_optical : tuple
        Observer locations for each observation in sim_obs_array_optical
    obs_array_radar : array
        Simulated radar observation data for the given reference solution.
    observer_codes_radar : tuple
        Observer locations for each observation in sim_obs_array_radar

    Raises
    ------
    ValueError
        If the length of simulated optical and radar observation times are both zero.
    ValueError
        If the length of simulated optical observation times and
        observation types are not the same.
    ValueError
        If the length of simulated radar observation times and
        observation types are not the same.
    ValueError
        If the reference solution does not contain the time.
    ValueError
        If the reference solution does not contain a full cartesian
        or cometary state.
    """
    obs_sigma_dict = {  'astrometry': 1, # arcsec
                        'occultation': 2.2e-3, # arcsec
                        'delay': 2, # microseconds
                        # 15 meter (conservative since it is 1m random + 10m systematic)
                        'delay_hera': 15*2 /299792458*1e6, # microseconds
                        'doppler': 0.5, # Hz
                        # 'doppler_hera': 0.1 mm/s to Hz
    }
    # check length of times is not zero
    (obs_array_optical, observer_codes_optical,
                    obs_array_radar, observer_codes_radar) = modified_obs_arrays
    optical_times = tuple(obs_array_optical[:, 0])
    radar_times = tuple(obs_array_radar[:, 0])
    if len(simulated_optical_obs_idx) == 0 and len(simulated_radar_obs_idx) == 0:
        raise ValueError("Must provide at least one observation time.")
    # check length of times and obs_types are the same
    if len(optical_times) != len(optical_obs_types):
        raise ValueError("Must provide the same number of optical",
                            "observation times and observation types.")
    if len(radar_times) != len(radar_obs_types):
        raise ValueError("Must provide the same number of radar",
                            "observation times and observation types.")
    # check that the reference solution has the requisite information and create integration body
    if 't' not in ref_sol:
        raise ValueError("Must provide a time for the initial solution.")
    consts = prop.Constants()
    nongrav_params = prop.NongravParamaters()
    nongrav_params.a1 = ref_ng_info['a1']
    nongrav_params.a2 = ref_ng_info['a2']
    nongrav_params.a3 = ref_ng_info['a3']
    nongrav_params.alpha = ref_ng_info['alpha']
    nongrav_params.k = ref_ng_info['k']
    nongrav_params.m = ref_ng_info['m']
    nongrav_params.n = ref_ng_info['n']
    nongrav_params.r0_au = ref_ng_info['r0_au']
    if all(key in ref_sol for key in ("x", "y", "z", "vx", "vy", "vz")):
        pos = [ref_sol['x'], ref_sol['y'], ref_sol['z']]
        vel = [ref_sol['vx'], ref_sol['vy'], ref_sol['vz']]
        target_body = prop.IntegBody("body_simulated_obs", ref_sol['t'], ref_sol['mass'],
                                        ref_sol['radius'], pos, vel,
                                        ref_cov, nongrav_params, consts)
    elif all(key in ref_sol for key in ("e", "q", "tp", "om", "w", "i")):
        ecc = ref_sol['e']
        peri_dist = ref_sol['q']
        time_peri = ref_sol['tp']
        omega = ref_sol['om']
        arg_peri = ref_sol['w']
        inc = ref_sol['i']
        cometary_elements = [ecc, peri_dist, time_peri,
                                omega*np.pi/180.0, arg_peri*np.pi/180.0, inc*np.pi/180.0]
        target_body = prop.IntegBody(de_kernel_path, "body_simulated_obs", ref_sol['t'],
                                        ref_sol['mass'], ref_sol['radius'], cometary_elements,
                                        ref_cov, nongrav_params, consts)
    else:
        raise ValueError("Must provide either a full cartesian or cometary",
                                "state for the initial solution.")
    # initialize past and future times and observer codes
    obs_times = np.array(optical_times + radar_times)
    observer_codes = observer_codes_optical + observer_codes_radar
    obs_types = tuple(optical_obs_types + radar_obs_types)
    sort_idx = np.argsort(obs_times)
    obs_times = obs_times[sort_idx]
    observer_codes = tuple(np.array(observer_codes, dtype=tuple)[sort_idx])
    obs_types = tuple(np.array(obs_types, dtype=tuple)[sort_idx])
    past_obs_idx = np.where(obs_times < ref_sol['t'])[0]
    future_obs_idx = np.where(obs_times > ref_sol['t'])[0]
    past_obs_exist = len(past_obs_idx) > 0
    future_obs_exist = len(future_obs_idx) > 0
    observer_info = get_observer_info(observer_codes)
    if past_obs_exist:
        t_eval_past = obs_times[past_obs_idx]
        tf_past = np.min(t_eval_past)
        observer_info_past = tuple(np.array(observer_info, dtype=tuple)[past_obs_idx])
    if future_obs_exist:
        t_eval_future = obs_times[future_obs_idx]
        tf_future = np.max(t_eval_future)
        observer_info_future = tuple(np.array(observer_info, dtype=tuple)[future_obs_idx])
    # initialize the propagator
    t_eval_utc = True
    eval_apparent_state = True
    converged_light_time = True
    prop_sim_past = prop.propSimulation("simulated_obs_past", ref_sol['t'],
                                        de_kernel, de_kernel_path)
    prop_sim_past.tEvalMargin = 1.0
    prop_sim_future = prop.propSimulation("simulated_obs_future", ref_sol['t'],
                                        de_kernel, de_kernel_path)
    prop_sim_future.tEvalMargin = 1.0
    if past_obs_exist:
        prop_sim_past.set_integration_parameters(tf_past, t_eval_past, t_eval_utc,
                                                eval_apparent_state, converged_light_time,
                                                observer_info_past)
        prop_sim_past.add_integ_body(target_body)
    if future_obs_exist:
        prop_sim_future.set_integration_parameters(tf_future, t_eval_future, t_eval_utc,
                                                    eval_apparent_state, converged_light_time,
                                                    observer_info_future)
        prop_sim_future.add_integ_body(target_body)
    # add events
    if events is not None:
        for event in events:
            t_event = event[0]
            dvx = event[1]
            dvy = event[2]
            dvz = event[3]
            multiplier = event[4]
            if t_event < ref_sol['t']:
                prop_sim_past.add_event(target_body, t_event, [dvx, dvy, dvz], multiplier)
            else:
                prop_sim_future.add_event(target_body, t_event, [dvx, dvy, dvz], multiplier)
    # propagate
    if past_obs_exist:
        prop_sim_past.preprocess()
        prop_sim_past.integrate()
    if future_obs_exist:
        prop_sim_future.preprocess()
        prop_sim_future.integrate()
    # get the propagated solution
    if past_obs_exist and future_obs_exist:
        apparent_states_past = np.array(prop_sim_past.xIntegEval)
        apparent_states_future = np.array(prop_sim_future.xIntegEval)
        apparent_states = np.vstack((apparent_states_past, apparent_states_future))
        radar_observations_past = np.array(prop_sim_past.radarObsEval)
        radar_observations_future = np.array(prop_sim_future.radarObsEval)
        radar_observations = np.vstack((radar_observations_past, radar_observations_future))
    elif past_obs_exist:
        apparent_states_past = np.array(prop_sim_past.xIntegEval)
        apparent_states = apparent_states_past
        radar_observations_past = np.array(prop_sim_past.radarObsEval)
        radar_observations = radar_observations_past
    elif future_obs_exist:
        apparent_states_future = np.array(prop_sim_future.xIntegEval)
        apparent_states = apparent_states_future
        radar_observations_future = np.array(prop_sim_future.radarObsEval)
        radar_observations = radar_observations_future
    optical_idx = 0
    radar_idx = 0
    for idx in range(len(obs_times)):
        typ = obs_types[idx]
        info_len = len(observer_info[idx])
        if info_len in {4, 7}:
            if typ != 'actual_obs_optical':
                if typ == 'actual_obs_to_sim_optical':
                    r_asc, dec = get_radec(apparent_states[idx])
                    obs_array_optical[optical_idx, 1] = r_asc
                    obs_array_optical[optical_idx, 2] = dec
                else:
                    r_asc, dec = get_radec(apparent_states[idx])
                    obs_array_optical[optical_idx, 1] = r_asc
                    obs_array_optical[optical_idx, 2] = dec
                    obs_array_optical[optical_idx, 3] = obs_sigma_dict[typ]
                    obs_array_optical[optical_idx, 4] = obs_sigma_dict[typ]
                    obs_array_optical[optical_idx, 5] = 0.0
                if noise:
                    ra_noise = np.random.normal(0, obs_array_optical[optical_idx, 3])
                    obs_array_optical[optical_idx, 1] += ra_noise
                    dec_noise = np.random.normal(0, obs_array_optical[optical_idx, 4])
                    obs_array_optical[optical_idx, 2] += dec_noise
            optical_idx += 1
        elif info_len in {9, 10}:
            if typ != 'actual_obs_radar':
                if typ == 'actual_obs_to_sim_radar':
                    if info_len == 9:
                        obs_array_radar[radar_idx, 1] = radar_observations[idx]
                    elif info_len == 10:
                        obs_array_radar[radar_idx, 2] = radar_observations[idx]
                else:
                    if info_len == 9:
                        obs_array_radar[radar_idx, 1] = radar_observations[idx]
                        obs_array_radar[radar_idx, 3] = obs_sigma_dict[typ]
                        obs_array_radar[radar_idx, 5] = np.nan
                    elif info_len == 10:
                        obs_array_radar[radar_idx, 2] = radar_observations[idx]
                        obs_array_radar[radar_idx, 4] = obs_sigma_dict[typ]
                        obs_array_radar[radar_idx, 5] = np.nan
                if noise:
                    if info_len == 9:
                        delay_noise = np.random.normal(0, obs_array_radar[radar_idx, 3])
                        obs_array_radar[radar_idx, 1] += delay_noise
                    elif info_len == 10:
                        doppler_noise = np.random.normal(0, obs_array_radar[radar_idx, 4])
                        obs_array_radar[radar_idx, 2] += doppler_noise
            radar_idx += 1
    return obs_array_optical, observer_codes_optical, obs_array_radar, observer_codes_radar

def create_simulated_obs_arrays(simulated_traj_info, real_obs_arrays, simulated_obs_start_time,
                                add_extra_simulated_obs, extra_simulated_obs_info, noise):
    """
    Creates observation arrays with simulated observations added to the real observations.

    Parameters
    ----------
    simulated_traj_info : tuple
        Tuple containing the nominal trajectory, covariance matrix, events, target radius,
        non-gravitational acceleration information, DE kernel, and DE kernel path.
    real_obs_arrays : tuple
        Tuple containing the real optical and radar observation arrays and observer codes.
    simulated_obs_start_time : float
        Time at which to start simulating observations.
    add_extra_simulated_obs : bool
        Flag indicating whether to add extra simulated observations.
    extra_simulated_obs_info : tuple
        Tuple containing the extra simulated optical and radar observation times and types.
    noise : bool
        Flag to add noise to the simulated observations.

    Returns
    -------
    obs_array_optical : array
        Real and simulated optical observation data for the given body
    observer_codes_optical : tuple
        Real and simulated observer locations for each observation in obs_array_optical
    obs_array_radar : array
        Real and simulated radar observation data for the given body
    observer_codes_radar : tuple
        Real and simulated observer locations for each observation in obs_array_radar
    """
    (x_nom, covariance, events, target_radius, nongrav_info,
        de_kernel, de_kernel_path) = simulated_traj_info
    (obs_array_optical, observer_codes_optical,
        obs_array_radar, observer_codes_radar) = real_obs_arrays
    if add_extra_simulated_obs:
        (extra_simulated_optical_obs_times, extra_simulated_optical_obs_types,
            extra_simulated_radar_obs_times,
            extra_simulated_radar_obs_types) = extra_simulated_obs_info
    optical_obs_times = obs_array_optical[:,0]
    radar_obs_times = obs_array_radar[:,0]
    simulated_optical_obs_idx = np.where(optical_obs_times >= simulated_obs_start_time)[0]
    simulated_optical_obs_times = tuple(optical_obs_times[simulated_optical_obs_idx])
    optical_obs_types = ['actual_obs_optical']*(len(optical_obs_times)-
                                                len(simulated_optical_obs_times))
    optical_obs_types += ['actual_obs_to_sim_optical']*len(simulated_optical_obs_times)
    if (add_extra_simulated_obs and extra_simulated_optical_obs_times is not None
            and extra_simulated_optical_obs_types is not None
            and len(extra_simulated_optical_obs_times) == len(extra_simulated_optical_obs_types)):
        num_extra_optical_obs = len(extra_simulated_optical_obs_times)
        # add extra simulated optical obs times and types
        simulated_optical_obs_times = simulated_optical_obs_times+extra_simulated_optical_obs_times
        optical_obs_types += extra_simulated_optical_obs_types
        # add extra rows to obs_array_optical
        extra_simulated_optical_obs_array = np.nan*np.ones((num_extra_optical_obs, 6))
        extra_simulated_optical_obs_array[:,0] = extra_simulated_optical_obs_times
        obs_array_optical = np.vstack((obs_array_optical, extra_simulated_optical_obs_array))
        # add indices to simulated_optical_obs_idx
        simulated_optical_obs_idx = np.hstack((simulated_optical_obs_idx,
                                                np.arange(len(optical_obs_times),
                                                len(optical_obs_times)+
                                                num_extra_optical_obs)))
        # add extra rows to observer_codes_optical
        extra_simulated_optical_observer_codes = tuple(['500']*num_extra_optical_obs)
        observer_codes_optical = observer_codes_optical + extra_simulated_optical_observer_codes
    # sort optical observations by time
    sort_idx = np.argsort(obs_array_optical[:,0])
    obs_array_optical = obs_array_optical[sort_idx]
    observer_codes_optical = tuple(observer_codes_optical[i] for i in sort_idx)
    optical_obs_types = tuple(optical_obs_types[i] for i in sort_idx)
    simulated_optical_obs_idx = []
    for i, typ in enumerate(optical_obs_types):
        if typ != 'actual_obs_optical':
            simulated_optical_obs_idx.append(i)
    simulated_radar_obs_idx = np.where(radar_obs_times >= simulated_obs_start_time)[0]
    simulated_radar_obs_times = tuple(radar_obs_times[simulated_radar_obs_idx])
    radar_obs_types = ['actual_obs_radar']*(len(radar_obs_times)-len(simulated_radar_obs_times))
    radar_obs_types += ['actual_obs_to_sim_radar']*len(simulated_radar_obs_times)
    if (add_extra_simulated_obs and extra_simulated_radar_obs_times is not None
            and extra_simulated_radar_obs_types is not None
            and len(extra_simulated_radar_obs_times) == len(extra_simulated_radar_obs_types)):
        num_extra_radar_obs = len(extra_simulated_radar_obs_times)
        # add extra simulated radar obs times and types
        simulated_radar_obs_times = simulated_radar_obs_times + extra_simulated_radar_obs_times
        radar_obs_types += extra_simulated_radar_obs_types
        # add extra rows to obs_array_radar
        extra_simulated_radar_obs_array = np.nan*np.ones((num_extra_radar_obs, 6))
        extra_simulated_radar_obs_array[:,0] = extra_simulated_radar_obs_times
        obs_array_radar = np.vstack((obs_array_radar, extra_simulated_radar_obs_array))
        # add indices to simulated_radar_obs_idx
        simulated_radar_obs_idx = np.hstack((simulated_radar_obs_idx,
                                                np.arange(len(radar_obs_times),
                                                len(radar_obs_times)+
                                                num_extra_radar_obs)))
        # add extra rows to observer_codes_radar
        extra_simulated_radar_observer_codes = []
        for i, typ in enumerate(extra_simulated_radar_obs_types):
            if typ == 'doppler':
                extra_simulated_radar_observer_codes.append((('-14','-14'),0,8560e6))
            else:
                extra_simulated_radar_observer_codes.append((('-14','-14'),0))
        extra_simulated_radar_observer_codes = tuple(extra_simulated_radar_observer_codes)
        observer_codes_radar = observer_codes_radar + extra_simulated_radar_observer_codes
    # sort radar observations by time
    sort_idx = np.argsort(obs_array_radar[:,0])
    obs_array_radar = obs_array_radar[sort_idx]
    observer_codes_radar = tuple(observer_codes_radar[i] for i in sort_idx)
    radar_obs_types = tuple(radar_obs_types[i] for i in sort_idx)
    simulated_radar_obs_idx = []
    for i, typ in enumerate(radar_obs_types):
        if typ != 'actual_obs_radar':
            simulated_radar_obs_idx.append(i)
    simulated_obs_ref_sol = x_nom.copy()
    simulated_obs_ref_sol['mass'] = 0.0
    simulated_obs_ref_sol['radius'] = target_radius
    simulated_obs_ref_cov = covariance.copy()
    simulated_obs_event = events if events is not None else None
    modified_obs_arrays = (obs_array_optical, observer_codes_optical,
                            obs_array_radar, observer_codes_radar)
    simulated_obs_info = _generate_simulated_obs(simulated_obs_ref_sol, simulated_obs_ref_cov,
                                                    nongrav_info, simulated_obs_event,
                                                    modified_obs_arrays,
                                                    simulated_optical_obs_idx, optical_obs_types,
                                                    simulated_radar_obs_idx, radar_obs_types,
                                                    de_kernel, de_kernel_path, noise)
    return simulated_obs_info
