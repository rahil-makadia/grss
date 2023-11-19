"""Unscented transformations for the GRSS orbit propagation code"""
import numpy as np

class SigmaPoints:
    """
    Class for representing sigma points for the unscented transformation.
    """

    def __init__(self, x_dict, cov, sp_type, alpha=None, beta=None, kappa=None):
        """
        Initialize a sigma point object.

        Parameters
        ----------
        x_dict : dict
            solution dictionary containing cometary/cartesian elements
        P : numpy.ndarray
            covariance matrix of the solution
        type : str
            type of sigma point. Choose from merwe or julier
        alpha : float
            scaling parameter for sigma points
        beta : float
            scaling parameter for sigma points
        kappa : float
            scaling parameter for sigma points
        """
        # References
        # .. [1] R. Van der Merwe "Sigma-Point Kalman Filters for Probabilitic
        #       Inference in Dynamic State-Space Models" (Doctoral dissertation)
        # .. [2] Julier, Simon J.; Uhlmann, Jeffrey "A New Extension of the Kalman
        #       Filter to Nonlinear Systems". Proc. SPIE 3068, Signal Processing,
        #       Sensor Fusion, and Target Recognition VI, 182 (July 28, 1997)
        self.x_dict = x_dict
        self.x = np.array([x_dict[key] for key in x_dict if key != "t"])
        self.cov = cov
        self.type = sp_type
        self.alpha = alpha
        self.beta = beta
        self.kappa = kappa
        self.n = len(self.x)

        self.w_m, self.w_c = self._get_weights()
        self.sigma_points, self.diff = self._get_points()
        self.sigma_points_dict = self._get_points_dict()
        return None

    def _get_weights(self):
        if self.type == "merwe":
            w_m, w_c = self._merwe_weights()
        elif self.type == "julier":
            w_m, w_c = self._julier_weights()
        else:
            raise ValueError("Invalid sigma point type")
        return w_m, w_c

    def _merwe_weights(self):
        if self.alpha is None or self.beta is None or self.kappa is None:
            raise ValueError(
                "Must specify alpha, beta, and kappa for Merwe sigma points"
            )
        lambda_ = self.alpha**2 * (self.n + self.kappa) - self.n
        scale = self.n + lambda_
        w_m = np.zeros(2 * self.n + 1)
        w_m[0] = lambda_ / scale
        w_m[1:] = 0.5 / scale + np.zeros((2 * self.n))
        w_c = w_m.copy()
        w_c[0] = w_c[0] + (1 - self.alpha**2 + self.beta)
        return w_m, w_c

    def _julier_weights(self):
        if self.kappa is None:
            raise ValueError("Must specify kappa for Julier sigma points")
        scale = self.n + self.kappa
        w_m = np.ones(2 * self.n + 1) * 1 / (2 * scale)
        w_m[0] = self.kappa / scale
        return w_m, w_m

    def _get_points(self):
        if self.type == "merwe":
            lambda_ = self.alpha**2 * (self.n + self.kappa) - self.n
            fac = self.n + lambda_
        elif self.type == "julier":
            fac = self.n + self.kappa
        else:
            raise ValueError("Invalid sigma point type")
        sqrt_cov = np.linalg.cholesky(fac * self.cov)
        sigma_points = np.zeros((2 * self.n + 1, self.n))
        diff = np.zeros((2 * self.n + 1, self.n))
        sigma_points[0] = self.x
        for i in range(self.n):
            plus = self.x + sqrt_cov.T[i]
            minus = self.x - sqrt_cov.T[i]
            sigma_points[i + 1] = plus
            diff[i + 1] = sqrt_cov.T[i]
            sigma_points[i + 1 + self.n] = minus
            diff[i + 1 + self.n] = -sqrt_cov.T[i]
        return sigma_points, diff

    def _get_points_dict(self):
        dict_list = [None] * (2 * self.n + 1)
        for i in range(2 * self.n + 1):
            dict_i = {"t": self.x_dict["t"]}
            for j, key in enumerate(self.x_dict):
                if key != "t":
                    dict_i[key] = self.sigma_points[i][j - 1]
            dict_list[i] = dict_i
        return dict_list

    def reconstruct(self, tr_sigma_points):
        """
        Reconstruct the mean and covariance from transformed sigma points.

        Parameters
        ----------
        tr_sigma_points : np.ndarray
            transformed sigma points

        Returns
        -------
        new_x : np.ndarray
            reconstructed mean
        new_cov : np.ndarray
            reconstructed covariance
        """
        new_x = np.zeros(self.n)
        new_cov = np.zeros((self.n, self.n))
        for i in range(2 * self.n + 1):
            new_x += self.w_m[i] * tr_sigma_points[i]
        for i in range(2 * self.n + 1):
            diff = tr_sigma_points[i] - new_x
            new_cov += self.w_c[i] * np.outer(diff, diff)
        return new_x, new_cov
