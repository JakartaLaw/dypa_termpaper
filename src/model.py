# Libraries
import numpy as np
from scipy import interpolate
from scipy.optimize import minimize_scalar
from collections import namedtuple

# Own modules
from agent import Agent
from parameters import parameters as par
from modules.utils import hermgauss_lognorm

StateTuple = namedtuple('statetuple', ['m', 'f', 'p'])
ChoiceTuple = namedtuple('choicetuple', ['kappa', 'i'])

StatePolicyTuple = namedtuple('statepolicytuple', ['m', 'f', 'p', 't'])
ChoicePolicyTuple = namedtuple('choicepolicytuple', ['kappa', 'i', 'c'])

class Model(Agent):

    def __init__(self, par, state, education_lvl):
        super().__init__(par=par, state=state, education_lvl=education_lvl)
        self.statespace = list()
        self.choicespace = list()
        self.policy = dict()

    def create_m_grid(self):
        # Grid assets
        grid_m_temp = np.linspace(self.par.m_min, self.par.m_max**self.par.m_tuning, self.par.Nm)
        grid_m = grid_m_temp ** (1/self.par.m_tuning)
        return grid_m

    def create_f_grid(self):
        # Grid financial knowledge
        grid_f = np.linspace(self.par.f_min, self.par.f_max, self.par.Nf)
        return grid_f

    def create_p_grid(self):
        grid_p = np.linspace(self.par.p_min, self.par.p_max, self.par.Np)
        return grid_p

    def create_statespace(self):
        '''grid over m, f, p'''
        m_grid = self.create_m_grid()
        f_grid = self.create_f_grid()
        p_grid = self.create_p_grid()

        for m in m_grid:
            for f in f_grid:
                for p in p_grid:
                    self.statespace.append(StateTuple(m, f, p))

    def create_choicespace(self):
        '''grid over i, kappa'''
        for i in [0, 1]:
            for kappa in [0, 55]:
                self.choicespace.append(ChoiceTuple(kappa = kappa, i = i))

    @staticmethod
    def create_gauss_hermite():
        par.nu_y, par.nu_y_w = hermgauss_lognorm(par.Nnu_y, par.sigma_nu_y)

    # For Jeppe
    @staticmethod
    def create_interp(x_vals, y_vals):
        """Returns function which interpolates"""
        return interpolate.interp1d(x_vals, y_vals, kind='linear', fill_value = "extrapolate")

    def create_grids():
        pass

    def V_integrate(self, c, choice, t):
        '''Calculates E_t(V_t+1) via brute force looping'''
        V_fut = 0
        for j_psi in range(1, self.par.Npsi):
            for i_xi in range(1, self.par.Nxi):
                for k_eps in range(1, self.par.Neps):
                    self.update_f(choice.i) # updateting to f_t+1
                    r = self.r() # Calculate return today

                    #interest factor
                    interest_factor = self.R_tilde(choice.kappa, shock=self.par.eps[k_eps])
                    assets = self.a(c, choice.i, choice.kappa, t)

                    income = par.xi[i_xi] * (par.G * self.state.p *  par.psi[j_psi] + self.par.age_poly[t+1] +  self.par.age_poly[t])

                    integrand = interest_factor * assets + income

                    V = self.par.psi_w[j_psi] * self.par.xi_w[i_xi] * self.par.eps_w[k_eps] * self.V_plus_interp(integrand) # GH weighting

                    V_fut += V

        return(V_fut)

    def find_V(self, choice, t):
        '''Find optimal c for all states for given choices (i,kappa) in period t'''

        if t == self.par.max_age:
            Vfunc = lambda c: self.utility(c, t)
        else:
            Vfunc = lambda c: self.utility(c, t) + self.par.beta * self.par.mortality[t] * self.V_integrate(c, choice, t)

        # Optimizer
        # Convert function to negative for minimization
        Vfunc_neg = lambda x: -Vfunc(x)

        # c) Find optimum
        res = minimize_scalar(Vfunc_neg, tol = self.par.tolerance
                              , bounds = [1,self.state.m], method = "bounded")

        Vstar, Cstar = float(-res.fun), float(res.x)

        return Vstar, Cstar

    def find_V_for_choices(t):
        '''Find optimal V for all i, k in period t'''

        #initialize V and C
        Vstar, Cstar = - np.inf, None

        for choice in self.choicespace:
            V, C = self.find_V(choice, t)
            if V > Vstar:
                _c = choice
                Cstar = Cstar

        spt = StatePolicyTuple(m = self.state.m, f = self.state.f, p = self.state.p, t=t)
        cpt = ChoicePolicyTuple(kappa = _c.kappa, i = _c.i, c = Cstar)
        self.policy[spt] = cpt


    def create_V_interp(self, Vstar, t):
        self.V_plus_interp = self.create_interp(self.par.grid_M, Vstar[t+1])

    def initialize_Vstar(self):
        Vstar = dict()
        Vstar[self.par.max_age] = np.array([self.utility(m, self.par.max_age) for m in self.par.grid_M])
        return Vstar

    def solve(self):
        # create state_space grid

        # create choice_space grid (over i og kappa)

        # V in period T
        Vstar = self.initialize_Vstar()
        # backwards loop:

        # 1) (V_star_interpolant) interpolant over næste periode mellem m_grid og v_star_t+1
        for t in reversed(range(par.start_age, par.max_age)):
            self.create_V_interp(Vstar)
            for s in statespace:
                self.state = s
                self.find_V_for_choices(t)

        # 2) optimér mht c
