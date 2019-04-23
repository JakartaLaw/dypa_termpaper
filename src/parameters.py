# Singleton Pattern
import numpy as np
from modules.utils import Struct

parameters = Struct()

parameters.delta = 0.04 # value not from paper (depreciation on financial knowledge)
parameters.mortatility = 0.03 # value not from paper ('p' in paper)
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

#utility function parameters
# is initialized within agent

# Asset grid
parameters.a_min = 0
parameters.a_max = 500000 # Arbitrary set by Ditlev
parameters.Na = 40 # From paper
parameters.a_tuning = 0.3 # From paper

# Consumption grid
parameters.Ca = 15

# Schock grid (mu)
parameters.mu_min = -10
parameters.mu_max = 10
parameters.Nmu = 25

parameters.rho_u = 0.96
# Financial grid
## We let i be binary, so one can only accumulate max_age-start_age divided by the corresponding delta
parameters.f_min = 0
parameters.f_max = parameters.max_age - parameters.start_age
parameters.Nf = 25 # From paper

# Gauss Hermite
parameters.Nnu_y = 8
parameters.sigma_nu_y = 1 # Set arbitrary for now
