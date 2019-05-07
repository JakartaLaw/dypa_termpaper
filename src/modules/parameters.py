# Singleton Pattern
import numpy as np
from modules.utils import Struct, hermgauss_lognorm, create_mortality

from numba import jitclass
from numba import int32, float32, double  # import the types

spec = [
    ("beta", float32,),
    ("delta", float32,),
    ("start_age", int32,),
    ("retire_age", int32,),
    ("max_age", int32,),
    ("mortality", float32[:],),
    ("rho_u", float32,),
    ("m_min", float32,),
    ("m_max", float32,),
    ("Nm", int32,),
    ("m_tuning", float32,),
    ("p_min", float32,),
    ("p_max", float32,),
    ("Np", int32,),
    ("f_min", float32,),
    ("f_max", float32,),
    ("Nf", int32,),
    ("r_max", float32,),
    ("r_min", float32,),
    ("R_bar", float32,),
    ("G", float32,),
    ("tolerance", float32,),
    ("psi", float32[:, :],),
    ("xi", float32[:, :],),
    ("eps", float32[:, :],),
    ("n", float32[:]),
    ("age_poly", float32[:]),
]


@jitclass(spec)
class Parameters(object):
    def __init__(
        self,
        beta,
        delta,
        start_age,
        retire_age,
        max_age,
        mortality,
        rho_u,
        m_min,
        m_max,
        Nm,
        m_tuning,
        p_min,
        p_max,
        Np,
        f_min,
        f_max,
        Nf,
        r_max,
        r_min,
        R_bar,
        G,
        tolerance,
        psi,
        xi,
        eps,
        n,
        age_poly,
    ):
        self.beta = beta
        self.delta = delta
        self.start_age = start_age
        self.retire_age = retire_age
        self.max_age = max_age
        self.mortality = mortality
        self.rho_u = rho_u
        self.m_min = m_min
        self.m_max = m_max
        self.Nm = Nm
        self.m_tuning = m_tuning
        self.p_min = p_min
        self.p_max = p_max
        self.Np = Np
        self.f_min = f_min
        self.f_max = f_max
        self.Nf = Nf
        self.r_max = r_max
        self.r_min = r_min
        self.R_bar = R_bar
        self.G = G
        self.tolerance = tolerance
        self.psi = psi
        self.xi = xi
        self.eps = eps
        self.n = n
        self.age_poly = age_poly
