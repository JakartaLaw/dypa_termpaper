# Libraries
import numpy as np
from scipy import interpolate
from scipy.optimize import minimize_scalar
from scipy.interpolate import LinearNDInterpolator
from collections import namedtuple
import datetime
# Own modules
from agent import utility, update_f, calc_a, R_tilde
from parameters import parameters as par
from modules.utils import hermgauss_lognorm, Struct
from modules.stategrid import create_statespace
from modules.namedtuples import StateTuple, ChoiceTuple
from modules.consumptionpreference import create_consumption_preference_array
from modules.agepolynomial import create_age_poly_array
# order: ['m', 'f', 'p', 't']
# order: ['c', 'kappa', 'i']

class Model():

    def __init__(self, par, education_lvl):

        assert education_lvl in ['<HS', 'HS', 'College']
        # create consumption preference
        par.n = create_consumption_preference_array(education_lvl)
        par.age_poly = create_age_poly_array(education_lvl)
        #instantiation
        self.par = par

    # For Jeppe
    @staticmethod
    def create_interp(statespace, Vstar):
        """Returns function which interpolates"""
        return LinearNDInterpolator(statespace, Vstar)


    def V_integrate(self, choice, state, par, interpolant):
        '''Calculates E_t(V_t+1) via brute force looping'''
        V_fut = 0.0
        for j_psi in range(1, par.Npsi):
            for i_xi in range(1, par.Nxi):
                for k_eps in range(1, par.Neps):

                    state = update_f(choice, state, par) # updateting to f_t+1

                    #interest factor
                    interest_factor = R_tilde(choice, shock=par.eps[k_eps])
                    assets = calc_a(choice, state, par)

                    income = par.xi[i_xi] * (par.G * state.p * par.psi[j_psi] + par.age_poly[t+1] +  par.age_poly[t])

                    integrand = interest_factor * assets + income

                    V = par.psi_w[j_psi] * par.xi_w[i_xi] * par.eps_w[k_eps] * self.V_plus_interp(integrand) # GH weighting

                    V_fut += V

        return V_fut

    def find_V(self, i, kappa, state, interpolant):
        '''Find optimal c for all states for given choices (i,kappa) in period t'''

        if state.t == self.par.max_age - 1:
            Vfunc = lambda c: utility(ChoiceTuple(c, kappa, i), state, self.par)
        else:
            Vfunc = lambda c: utility(ChoiceTuple(c, kappa, i), state, self.par) + \
            self.par.beta * self.par.mortality[state.t] * \
            self.V_integrate(ChoiceTuple(c, kappa, i), state, self.par, interpolant)

        # Optimizer
        # Convert function to negative for minimization
        Vfunc_neg = lambda x: -Vfunc(x)

        # c) Find optimum
        res = minimize_scalar(Vfunc_neg, tol = self.par.tolerance
                              , bounds = [1, state.m], method = "bounded")

        Vstar, Cstar = float(-res.fun), float(res.x)

        return Vstar, Cstar

    def find_V_for_choices(self, state, interpolant):
        '''Find optimal V for all i, k in period t'''

        #initialize V and C
        V, C = - np.inf, None
        choices = np.array(((0.0, 0.0), (1., 0.), (0., 0.55), (1., 0.55)))

        for i, kappa in choices:

            #using notation _V, _C for best V, C for given kappa, i
            _V, _C = self.find_V(i, kappa, state, interpolant)
            if _V > V:
                C = _C

        return V, C


    @staticmethod
    def create_Vstar(statespace):
        """Creates data container for Vstar"""
        return np.empty(statespace.shape[0])

    @staticmethod
    def create_Cstar(statespace):
        """Creates data container for Vstar"""
        return np.empty((statespace.shape[0], 3))

    def initialize_Vstar(self, statespace):
        Vstar = self.create_Vstar(statespace)
        for state_index, s in enumerate(statespace):
            m, f, p, t = s[0], s[1], s[2], self.par.max_age
            state, choice = StateTuple(m, f, p, t), ChoiceTuple(m, 0.0, 0.0) #consuming all
            Vstar[state_index] = utility(choice, state, self.par)
        return Vstar

    def initialize_Cstar(self, statespace):
        Cstar = self.create_Cstar(statespace)
        for state_index, s in enumerate(statespace):
            m, f, p, t = s[0], s[1], s[2], self.par.max_age
            choice = ChoiceTuple(m, 0.0, 0.0) #consuming all
            Cstar[state_index] = choice
        return Cstar

    def solve(self):
        # create state_space grid values. order (m, f, p)
        statespace = create_statespace(self.par)

        V_solution, C_solution = dict(), dict()

        # V, C in period T
        Vstar = self.initialize_Vstar(statespace)
        Cstar = self.initialize_Cstar(statespace)
        # backwards loop:

        # 1) (V_star_interpolant) interpolant over næste periode mellem m_grid og v_star_t+1
        for t in reversed(range(par.start_age, par.max_age)):
            print('Solution at time step t: ', t, ', time is: ', datetime.datetime.utcnow())

            V_solution[t], C_solution = Vstar, Cstar

            Vstar_plus, Cstar_plus = self.create_Vstar(statespace), self.create_Cstar(statespace)
            interpolant = self.create_interp(statespace, Vstar)

            for s_ix, s in enumerate(statespace):
                m, f, p = s[0], s[1], s[2]
                state = StateTuple(m, f, p, t)
                V, C = self.find_V_for_choices(state, interpolant)
                Vstar_plus[s_ix], Cstar_plus[s_ix] = V, C

            Vstar, Cstar = Vstar_plus, Cstar_plus

        return V_solution, C_solution
#%%
