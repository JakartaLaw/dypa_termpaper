
##############################################################################################################
# Import libraries
from model import Model
from parameters import parameters as par
from modules.utils import Struct # GH functions
from collections import namedtuple
from modules.stategrid import create_statespace
import pandas as pd
import time

# Run Model

start = time.time()

m = Model()
V_sol, C_sol = m.solve(par)

end = time.time()
print(end - start)
