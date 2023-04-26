from sherpa.models import model
import numpy as np
import astropy.units as u
from astropy.constants import c
from agnpy.spectra import BrokenPowerLaw
from agnpy.synchrotron import Synchrotron
from agnpy.compton import SynchrotronSelfCompton


class AgnpySSC(model.RegriddableModel1D):
    """Wrapper of agnpy's synchrotron and SSC classes. 
    A broken power law is assumed for the electron spectrum.
    To limit the span of the parameters space, we fit the log10 of the parameters 
    whose range is expected to cover several orders of magnitudes (normalisation, 
    gammas, size and magnetic field of the blob). 
    """

    def __init__(self, name="ssc"):

        # EED parameters
        self.log10_k_e = model.Parameter(name, "log10_k_e", -2.0, min=-20.0, max=10.0)
        self.p1 = model.Parameter(name, "p1", 2.1, min=-2.0, max=5.0)
        self.p2 = model.Parameter(name, "p2", 3.1, min=-2.0, max=5.0)
        self.log10_gamma_b = model.Parameter(name, "log10_gamma_b", 3, min=1, max=10)
        self.log10_gamma_min = model.Parameter(name, "log10_gamma_min", 1, min=0, max=10)
        self.log10_gamma_max = model.Parameter(name, "log10_gamma_max", 5, min=0, max=10)
        # source general parameters
        self.z = model.Parameter(name, "z", 0.1, min=0.01, max=5)
        self.d_L = model.Parameter(name, "d_L", 1e27, min=1e25, max=1e33, units="cm")
        # emission region parameters
        self.delta_D = model.Parameter(name, "delta_D", 10, min=0, max=100)
        self.log10_B = model.Parameter(name, "log10_B", -2, min=-4, max=2)
        self.t_var = model.Parameter(name, "t_var", 600, min=10, max=np.pi * 1e7, units="s")

        model.RegriddableModel1D.__init__(
            self,
            name,
            (   self.log10_k_e,
                self.p1,
                self.p2,
                self.log10_gamma_b,
                self.log10_gamma_min,
                self.log10_gamma_max,
                self.z,
                self.d_L,
                self.delta_D,
                self.log10_B,
                self.t_var,
            ),
        )

    def calc(self, pars, x):
        """evaluate the model calling the agnpy functions"""
        (   log10_k_e,
            p1,
            p2,
            log10_gamma_b,
            log10_gamma_min,  
            log10_gamma_max,
            z,
            d_L,
            delta_D,
            log10_B,
            t_var,
        ) = pars
        # add units, scale quantities
        x *= u.Hz
        k_e = 10 ** log10_k_e * u.Unit("cm-3")
        gamma_b = 10 ** log10_gamma_b
        gamma_min = 10 ** log10_gamma_min
        gamma_max = 10 ** log10_gamma_max
        B = 10 ** log10_B * u.G
        d_L *= u.cm
        R_b = c.to_value("cm s-1") * t_var * delta_D / (1 + z) * u.cm
        sed_synch = Synchrotron.evaluate_sed_flux(
            x,
            z,
            d_L,
            delta_D,
            B,
            R_b,
            BrokenPowerLaw,
            k_e,
            p1,
            p2,
            gamma_b,
            gamma_min,
            gamma_max,
        )
        sed_ssc = SynchrotronSelfCompton.evaluate_sed_flux(
            x,
            z,
            d_L,
            delta_D,
            B,
            R_b,
            BrokenPowerLaw,
            k_e,
            p1,
            p2,
            gamma_b,
            gamma_min,
            gamma_max,
        )
        return sed_synch + sed_ssc  #  return sed_synch + sed_ssc