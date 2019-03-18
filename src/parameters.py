# Singleton Pattern

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
