
##############################################################################################################
# Import libraries
import agent
from model import Model
from parameters import parameters as par
from modules.utils import Struct # GH functions
from collections import namedtuple
from modules.stategrid import create_statespace
import pandas as pd
from modules.consumptionpreference import create_consumption_preference_array
from modules.agepolynomial import create_age_poly_array
import numpy as np
from modules import stategrid

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")

# age_poly_hs = create_age_poly_array('HS')
# n_hs = create_consumption_preference_array('HS')
#
# len(age_poly_hs)
# len(n_hs)
# len(np.ones(90))
#
# n_hs
# age_poly_hs[89]
# len(age_poly_hs)

# Investigate state-space
statespace = create_statespace(par)
statespace

# # Run Model
m = Model()
V_sol, C_sol = m.solve(par)
# c, i, k
C_sol[80]


# Test R return
#
agent.R_riskful(20, par, 1)


agent.R_riskfree(par)
test_sp = create_statespace(par)
test_sp[30:32]

test_sp[30:32]
