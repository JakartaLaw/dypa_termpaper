# Singleton Pattern
import numpy as np
from modules.utils import Struct, hermgauss_lognorm, create_mortality

parameters = Struct()

parameters.beta = 0.98
parameters.delta = 0.04 # value not from paper (depreciation on financial knowledge)
parameters.c_min = 0 #value not from paper (max cash_on_hand for government transfer)
parameters.sigma_y_eps = 1 #value not from paper
parameters.sigma_y_v = 1 #value not from paper
parameters.sigma_o_eps = 1 #value not from paper
parameters.sigma_o_v = 1 #value not from paper
parameters.rho_y_e = 0.5 #value not from paper
parameters.rho_o_e = 0.5 ##value not from paper
parameters.c_d = 4 #value not from paper, Fixed cost of investing in financial knowledge
parameters.start_age = 25
parameters.retire_age = 65
parameters.max_age = 90
#parameters.mortality = np.array([0.03 for i in range(parameters.max_age + 1)])
parameters.mortality = create_mortality()


#utility function parameters
# is initialized within agent

# Asset grid
parameters.m_min = 2
parameters.m_max = 500000 # Arbitrary set by Ditlev
parameters.Nm = 40 # From paper
parameters.m_tuning = 0.3 # From paper

# Consumption grid
parameters.Ca = 15

# Schock grid (permanent income p)
parameters.p_min = -10
parameters.p_max = 10
parameters.Np = 25

parameters.rho_u = 0.96
# Financial grid
## We let i be binary, so one can only accumulate max_age-start_age divided by the corresponding delta
parameters.f_min = 0
parameters.f_max = parameters.max_age - parameters.start_age
parameters.Nf = 25 # From paper

# Gauss Hermite

# Number of nodes
parameters.Npsi = 8
parameters.Nxi = 8
parameters.Neps = 8

# Weights and nodes
parameters.psi_sigma = 0.1
parameters.xi_sigma = 0.1
parameters.eps_sigma = 0.1

parameters.eps, parameters.eps_w = hermgauss_lognorm(parameters.Neps, parameters.eps_sigma)
parameters.xi, parameters.xi_w = hermgauss_lognorm(parameters.Nxi, parameters.xi_sigma)
parameters.psi, parameters.psi_w = hermgauss_lognorm(parameters.Npsi, parameters.psi_sigma)

parameters.Nnu_y = 8
parameters.sigma_nu_y = 1 # Set arbitrary for now
parameters.r_max = 0.08
parameters.r_min = 0.00
parameters.R_bar = 1.02
parameters.G = 1.03

parameters.M_max = 100000
parameters.NM = 100
parameters.grid_M = np.linspace(0, parameters.M_max, parameters.NM)

parameters.tolerance = 10**(-6)
