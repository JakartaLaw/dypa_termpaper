# Libraries
import numpy as np
from scipy import interpolate

# Own modules
from agent import Agent
from parameters import parameters as par
from modules.utils import hermgauss_lognorm

class Model(Agent):

    def __init__(self, par, state, education_lvl):
        super().__init__(par=par, state=state, education_lvl=education_lvl)

    def create_a_grid(self):
        # Grid assets
        grid_a_temp = np.linspace(self.par.a_min, self.par.a_max**self.par.a_tuning, self.par.Na)
        grid_a = grid_a_temp ** (1/self.par.a_tuning)

    def create_f_grid(self):
        # Grid financial knowledge
        grid_f = np.linspace(self.par.f_min, self.par.f_max, self.par.Nf)

    def create_c_grid(self):
        pass

    def create_mu_grid(self):
        grid_mu = np.linspace(self.par.mu_min, self.par.mu_max, self.par.Nmu)

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

    def find_V(choice, t):

        for state in self.statespace:

            self.state = state

            if self.state.t == self.par.max_age:
                Vfunc = self.utility
            else:
                Vfunc = lambda c: choice.self.utility(c, t) + self.par.beta * self.par.mortality[t] * V_integrate(c, choice, t)

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

            find_V()

        # 2) optimér mht c
