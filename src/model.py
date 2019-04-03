# Libraries
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

# Own modules
from agent import Agent
from parameters import parameters as par
from modules.utils import hermgauss_lognorm

class Model(Agent):

    def __init__(self, par, state):
        super().__init__(par=par, state=state)


    def create_a_grid(self):
        # Grid assets
        grid_a_temp = np.linspace(self.par.a_min, self.par.a_max**self.par.a_tuning, self.par.Na)
        grid_a = grid_a_temp ** (1/self.par.a_tuning)

    def create_f_grid(self):
        # Grid financial knowledge
        grid_f = np.linspace(self.par.f_min, self.par.f_max, self.par.Nf)

    def create_c_grid(self):
        pass

    @staticmethod
    def create_gauss_hermite():
        par.nu_y, par.nu_y_w = hermgauss_lognorm(par.Nnu_y, par.sigma_nu_y)

    @staticmethod
    def create_age_poly_dict():
        '''Returns a dictionary for the 3 education groups with income profiles (in $1,000) for ages between
        25 and 65, where the first element of the list corresponds to the income at age 25. This data is
        based on the graph on p. 438'''

        # Load data
        age = np.array([25,30,35,40,45,50,55,60,65]).reshape(-1,1)
        col = np.array([40, 50, 57, 62, 66, 67, 65, 61, 54])
        hs = np.array([32, 37, 42, 46, 49, 50, 48, 44, 37])
        lhs = np.array([26, 31, 34, 38, 39, 39, 37, 34, 31])
        educ = [lhs, hs, col]

        # Setup degree 3 regression and output dict format
        lm = LinearRegression()
        poly = PolynomialFeatures(3)
        age_poly =  poly.fit_transform(age)

        names = ['<HS', 'HS', 'College']
        pred_income = dict()

        # Regression and prediction
        for i in range(len(educ)):
            # Fit
            fit = lm.fit(age_poly, educ[i])
            pred_income[names[i]] = []
            for t in range(25, 65):
                age_t = np.array(t).reshape(-1,1)
                age_t_poly = poly.fit_transform(age_t)
                pred_t = float(lm.predict(age_t_poly))
                pred_income[names[i]].append(pred_t)
        return(pred_income)


    # For Jeppe

    def solve(self, par):
        # note: possibly use numpy 2d array for solving this with numba integration
        # note: we use j, as index variable in loops

        # 1. allocation solution struct and cells
        sol = Struct()
        sol.m = dict()
        sol.c = dict()

            # State Space
            for a_j in a_grid:
                for f_j in f_grid:

                    # control variable space
                    c_grid, j_gti = create_c_grid(a_i)
                    j_grid = create_j_grid()
                    for c_i in c_grid:

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
