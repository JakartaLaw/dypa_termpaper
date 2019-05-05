# Singleton Pattern
import numpy as np
from modules.parameters import Parameters
from modules.utils import Struct, hermgauss_lognorm, create_mortality

from numba import int32, float32    # import the types

mortality = create_mortality
psi = hermgauss_lognorm(n=8, sigma=0.1)
xi = hermgauss_lognorm(n=8, sigma=0.1)
eps = hermgauss_lognorm(n=8, sigma=0.1)
max_age, start_age = 90, 25

parameters_initial = {
    "beta" : 0.98,
    "delta" : 0.04,
    "start_age" : 87,
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
    "eps" : eps,
    "xi" : xi,
    "psi" : psi,
    "r_max" : 0.08,
    "r_min" : 0.00,
    "R_bar" : 1.02,
    "G" : 1.03,
    "tolerance" : 10**(-6),
}

parameters = Parameters(**parameters_initial)










#
# parameters.beta = 0.98
# parameters.delta = 0.04 # value not from paper (depreciation on financial knowledge)
#
# #parameters.start_age = 25
# parameters.start_age = 87 #ensuring model works
# parameters.retire_age = 65
# parameters.max_age = 90
# parameters.mortality = create_mortality()
# parameters.rho_u = 0.96
#
# # Asset grid
# parameters.m_min = 2
# parameters.m_max = 500000 # Arbitrary set by Ditlev
# parameters.Nm = 4 # From paper (40)
# parameters.m_tuning = 0.3 # From paper
#
# # Schock grid (permanent income p)
# parameters.p_min = -10
# parameters.p_max = 10
# parameters.Np = 4 # From paper 25
#
# # Financial grid
# ## We let i be binary, so one can only accumulate max_age-start_age divided by the corresponding delta
# parameters.f_min = 0
# parameters.f_max = parameters.max_age - parameters.start_age
# parameters.Nf = 4 # From paper 25
#
# # Gauss Hermite
#
# # Number of nodes
# parameters.Npsi = 8
# parameters.Nxi = 8
# parameters.Neps = 8
#
# # Weights and nodes
# parameters.psi_sigma = 0.1
# parameters.xi_sigma = 0.1
# parameters.eps_sigma = 0.1
#
# parameters.eps, parameters.eps_w = hermgauss_lognorm(parameters.Neps, parameters.eps_sigma)
# parameters.xi, parameters.xi_w = hermgauss_lognorm(parameters.Nxi, parameters.xi_sigma)
# parameters.psi, parameters.psi_w = hermgauss_lognorm(parameters.Npsi, parameters.psi_sigma)
#
# parameters.r_max = 0.08
# parameters.r_min = 0.00
# parameters.R_bar = 1.02
# parameters.G = 1.03
#
# parameters.tolerance = 10**(-6)
