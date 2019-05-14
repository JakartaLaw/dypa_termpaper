# Singleton Pattern
import numpy as np
from modules.parameters import Parameters
from modules.utils import Struct, hermgauss_lognorm, create_mortality
from modules.consumptionpreference import create_consumption_preference_array
from modules.agepolynomial import create_age_poly_array

from numba import int32, float32    # import the types

N_GH = 4
mortality = create_mortality()
psi = hermgauss_lognorm(n=N_GH, sigma=0.1)
xi = hermgauss_lognorm(n=N_GH, sigma=0.1)
eps = hermgauss_lognorm(n=N_GH, sigma=0.1)
max_age, start_age = 90, 87

age_poly_hs = create_age_poly_array('HS')
n_hs = create_consumption_preference_array('HS')

parameters_initial = {
    "beta" : 0.98,
    "delta" : 0.04,
    "start_age" : start_age,
    "retire_age" : 65,
    "max_age" : max_age,
    "mortality" : mortality,
    "rho_u" : 0.96,
    "m_min" : 2.0,
    "m_max" : 500000.0, # Arbitrary set by Ditlev
    "Nm" : 4, # From paper (40)
    "m_tuning" : 0.3, # From paper
    "p_min" : -10.0,
    "p_max" : 10.0,
    "Np" : 4, # From paper 25
    "f_min" : 0.0,
    "f_max" : max_age - start_age,
    "Nf" : 4, # From paper 25
    "r_max": 0.1,
    "r_min": 0,
    "R_bar": 1.02,
    "G" : 1.03,
    "tolerance" : 10**(-6),
    "psi" : psi,
    "xi" : xi,
    "eps" : eps
}

hs_params = {key : value for key, value in parameters_initial.items()}
hs_params["n"] = n_hs
hs_params["age_poly"] = age_poly_hs

parameters = Parameters(**hs_params)
