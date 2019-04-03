# Libraries
import numpy as np

# Own modules
from agent import Agent
from parameters import parameters as par
from modules.utils import hermgauss_lognorm



class Model(Agent):

    def __init__(self, par, state):
        super().__init__(par=par, state=state)


    @staticmethod
    def create_a_grid():
        # Grid assets
        grid_a_temp = np.linspace(par.a_min, par.a_max**par.a_tuning, par.Na)
        grid_a = grid_a_temp ** (1/par.a_tuning)

    @staticmethod
    def create_f_grid():
        # Grid financial knowledge
        grid_f = np.linspace(par.f_min, par.f_max, par.Nf)

    @staticmethod
    def create_gauss_hermite():
        










# %% For Jeppe

@classmethod
def solve(cls, par):
    # 1. allocation solution struct and cells
    sol = Struct()
    sol.m = dict()
    sol.c = dict()

    # 2. last period (= consume all)

    sol.m[par.T] = np.linspace(0, par.a_max, par.Na)
    sol.c[par.T] = np.linspace(0, par.a_max, par.Na)

    # 2 Before last Period
    for t in reversed(range(1,par.T)): # Start in period T-1

        #a) Interpolant
        par.c_plus_interp = interpolate.interp1d(sol.m[t+1], sol.c[t+1], kind='linear', fill_value = "extrapolate")

        #b) EGM
        sol_c, sol_m = cls.EGM(sol, t, par.c_plus_interp, par)

        #c) Add zero Consumption
        sol.m[t] = np.append(par.a_min[t], sol_m)
        sol.c[t] = np.append(0, sol_c)

    return(sol)

@staticmethod
def simulate(par,sol):
    sim = Struct()

    #1) Allocate
    sim.m = np.empty((par.simN,par.simT)) * np.nan
    sim.c = np.empty((par.simN,par.simT)) * np.nan
    sim.a = np.empty((par.simN,par.simT)) * np.nan
    sim.p = np.empty((par.simN,par.simT)) * np.nan
    sim.y = np.empty((par.simN,par.simT)) * np.nan

    #2) Shocks
    shocki = np.random.choice(par.Nshocks, par.simN * par.simT, True, par.w)
    shocki = np.reshape(shocki, (par.simN, par.simT)) # Reshape into matrix format
    sim.psi = par.psi_vec[shocki] # Allocate the shocks to the different indices in the shocki matrix
    sim.xi = par.xi_vec[shocki]

    #3) Initialize values
    sim.m[:,0] = par.sim_mini
    sim.p[:,0] = 0

    #4) Simulation
    for t in range(1,par.T):
        if t < par.T - 1:
            print(t, end = "\r")
        else:
            print ('done') # Done

        if par.simlifecycle == 0: # Infinite horizon
            c_interp = interpolate.interp1d(sol.m[1], sol.c[1], kind='linear', fill_value = "extrapolate")
        else:
            c_interp = interpolate.interp1d(sol.m[t], sol.c[t], kind='linear', fill_value = "extrapolate")

        sim.c[:,t-1] = c_interp(sim.m[:,t-1])
        sim.a[:,t-1] = sim.m[:,t-1] - sim.c[:,t-1]
        # print(sim.c[:,t-1])

        if t < par.simT:
            if t > par.TR: # Retired
                sim.m[:,t] = par.R * sim.a[:,t-1] / (par.G * par.L[t-1]) + 1
                sim.p[:,t] = np.log(par.G) + np.log(par.L[t-1]) + sim.p[:,t-1]
                sim.y[:,t] = sim.p[:,t]
            else:
                sim.m[:,t] = par.R * sim.a[:,t-1] / (par.G * par.L[t-1] * sim.psi[:,t]) + sim.xi[:,t]
                sim.p[:,t] = np.log(par.G) + np.log(par.L[t-1]) + sim.p[:,t-1] + np.log(sim.psi[:,t])
                sim.y[:,t] = sim.p[:,t] + np.log(sim.xi[:,t])

    # 5. Re-normalize
    sim.P = np.exp(sim.p)
    sim.Y = np.exp(sim.y)
    sim.M = sim.m * sim.P
    sim.C = sim.c * sim.P
    sim.A = sim.a * sim.P

    return sim
