# Singleton Pattern
import numpy as np
from modules.utils import Struct, hermgauss_lognorm, create_mortality

from numba import jitclass
from numba import int64, double, float32 # import the types

spec = [
    ("beta", double,),
    ("delta", double,),
    ("start_age", int64,),
    ("retire_age", int64,),
    ("max_age", int64,),
    ("mortality", double[:],),
    ("rho_u", double,),
    ("m_min", double,),
    ("m_max", double,),
    ("Nm", int64,),
    ("m_tuning", double,),
    ("p_min", double,),
    ("p_max", double,),
    ("Np", int64,),
    ("f_min", double,),
    ("f_max", double,),
    ("Nf", int64,),
    ("r_max", double,),
    ("r_min", double,),
    ("R_bar", double,),
    ("G", double,),
    ("tolerance", double,),
    ("psi", double[:, :],),
    ("xi", double[:, :],),
    ("eps", double[:, :],),
    ("n", double[:]),
    ("age_poly", double[:]),
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
