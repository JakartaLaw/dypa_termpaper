
##############################################################################################################
# Import libraries
from model import Model
from parameters import parameters as par
from modules.utils import Struct # GH functions
from collections import namedtuple
from modules.stategrid import create_statespace
import pandas as pd


# Run Model
m = Model()
V_sol, C_sol = m.solve(par)


V_sol[88]
C_sol[87]
