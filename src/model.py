from agent import Agent

class Model(Agent):

    def __init__(self, par, state):
        super().__init__(par=par, state=state)


    @staticmethod
    def create_grids():

        # Grid 1: Assets


    # For Jeppe
    def create_c_grid(self, a_i):
        """
        Parameters
        ==========
        a_i : float
            the maximum of the consumption grid, for given amount of assets
        """
        return np.linspace(0, a_i, self.par.Nc)

    @staticmethod
    def consumption_interp(sol, t):
        return interpolate.interp1d(sol.a[t+1], sol.c[t+1], kind='linear', fill_value = "extrapolate")

    @staticmethod
    def financial_interp(sol, t):
        return interpolate.interp1d(sol.f[t+1], sol.c[t+1], kind='linear', fill_value = "extrapolate")

    def solve(self, cls, par):
        # note: possibly use numpy 2d array for solving this with numba integration

        # 1. allocation solution struct and cells
        sol = Struct()
        sol.a = dict()
        sol.c = dict()
        sol.f = dict()

        # 2. last period (= consume all)

        sol.a[self.par.T] = np.linspace(0, self.par.a_max, self.par.Na)
        sol.c[self.par.T] = np.linspace(0, self.par.a_max, self.par.Na)

        # 2 Before last Period
        a_grid, f_grid = create_a_grid(), create_f_grid()
        for t in reversed(range(1,par.T)): # Start in period T-1

            # Interpolants
            c_plus_interp = self.consumption_interp(sol, t)
            f_plus_interp = self.financial_interp(sol, t)

            # State Space
            for a_i in a_grid:
                for f_i in f_grid()

                # control variable space
                c_grid, j_gti = create_c_grid(a_i)
                j_grid = create_j_grid()
                for c_i in c_grid:



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
