# Libraries
import numpy as np
from scipy import interpolate
from scipy.optimize import minimize_scalar, minimize
from scipy.interpolate import LinearNDInterpolator
from collections import namedtuple
import datetime
import time
# Own modules
from agent import utility, update_f, calc_a, R_tilde
from parameters import parameters as par
from modules.utils import hermgauss_lognorm, Struct
from modules.stategrid import create_grids, create_statespace
from modules.namedtuples import StateTuple, ChoiceTuple
from modules.consumptionpreference import create_consumption_preference_array
from modules.agepolynomial import create_age_poly_array
from modules.meshmapping import create_lookup_dict, create_mesh #used for interpolater
from modules.interp import Interpolate3D

# order: ['m', 'f', 'p', 't']
# order: ['c', 'kappa', 'i']

from numba import jit, njit

def V_integrate(c, kappa, i, m, f, p, t, par, interpolant):
    '''Calculates E_t(V_t+1) via brute force looping'''

    # global NUMBER_OF_ITERATIONS
    # NUMBER_OF_ITERATIONS = NUMBER_OF_ITERATIONS + 1

    V_fut = 0.0

    # Calculations that can be moved outside the loop
    state = update_f(i, f, par) # updateting to f_t+1
    assets = calc_a(c, i, kappa, m, par)

    # refactored loop to return a value and weight from gaus_hermite
    for psi, psi_w in par.psi:
        for xi, xi_w in par.xi:
            for eps, eps_w in par.eps:

                interest_factor = R_tilde(kappa, f, par, shock=eps)
                income = xi * (par.G * p * psi) + par.age_poly[t+1]

                # Future state values
                m_fut = interest_factor * assets + income
                p_fut = par.G * p * psi
                f_fut = np.float64(f)

                interp_ = interpolant.interpolate(m_fut, f_fut, p_fut)
                V = psi_w * xi_w * eps_w * interp_ # GH weighting
                if np.isnan(V) == True:
                    import pdb; pdb.set_trace()
                    print('interp', interp_, type(interp_))
                    print('types', type(m_fut), type(f_fut), type(p_fut))
                    print('m, f, p',m_fut, f_fut, p_fut)
                    print('types', type(psi), type(xi_w), type(eps_w))
                    print('psi, xi_w, eps_w', psi_w * xi_w * eps_w)
                    raise ValueError('errors')
                V_fut += V

    return V_fut


class Model():

    def find_V(self, i, kappa, m, f, p, t, par, interpolant):
        '''Find optimal c for all states for given choices (i,kappa) in period t'''

        #OLD CODE:
        if t == par.max_age - 1:
            n, rho_u = par.n[t], par.rho_u
            Vfunc = lambda c: utility(c, n, rho_u)
        else:
            n, rho_u = par.n[t], par.rho_u
            Vfunc = lambda c: utility(c, n, rho_u) + \
            par.beta * par.mortality[t] * \
            V_integrate(c, kappa, i, m, f, p, t, par, interpolant)

        #NEW CODE:
        # n, rho_u = par.n[t], par.rho_u
        # Vfunc = lambda c: utility(c, n, rho_u) + \
        # par.beta * par.mortality[t] * \
        # V_integrate(c, kappa, i, m, f, p, t, par, interpolant)

        # Optimizer
        # Convert function to negative for minimization
        Vfunc_neg = lambda x: -Vfunc(x)

        # c) Find optimum
        res = minimize_scalar(Vfunc_neg, tol = par.tolerance, bounds = [1, m], method = "bounded")
        if t == 89:
            print(res.fun, res.x)

        Vstar, Cstar = float(-res.fun), float(res.x)

        return Vstar, Cstar

    #@jit
    def find_V_for_choices(self, state, par, interpolant):
        '''Find optimal V for all i, k in period t'''

        #initialize V and C
        V, C = - np.inf, None
        choices = np.array(((0.0, 0.0), (1., 0.), (0., 0.55), (1., 0.55)))
        m, f, p, t = state.m, state.f, state.p, state.t
        for i, kappa in choices:

            #using notation _V, _C for best V, C for given kappa, i
            _V, _C = self.find_V(i, kappa, m, f, p, t, par, interpolant)
            if np.isnan(_V) == True:
                print(ChoiceTuple(_C, kappa, i), t)

            if _V > V:
                C, V = ChoiceTuple(_C, kappa, i), _V

        return V, C


    @staticmethod
    def create_Vstar(statespace):
        """Creates data container for Vstar"""
        return np.empty(statespace.shape[0])

    @staticmethod
    def create_Cstar(statespace):
        """Creates data container for Vstar"""
        return np.empty((statespace.shape[0], 3))

    def initialize_Vstar(self, statespace, par):
        '''Only applicable for final period T'''
        Vstar = self.create_Vstar(statespace)
        for state_index, s in enumerate(statespace):
            m, n, rho_u = s[0], par.n[89], par.rho_u
            try:
                _u = utility(m, n, rho_u)
                Vstar[state_index] = _u
                if np.isnan(_u) is True:
                    print("initialize v star", m, n, rho_u)
            except Exception as e:
                print(m, t, par.rho_u)
                raise Exception('the values was:', m, t, par.rho_u, e)
        return Vstar

    def initialize_Cstar(self, statespace, par):
        '''Only applicable for final period T'''
        Cstar = self.create_Cstar(statespace)
        for state_index, s in enumerate(statespace):
            m = s[0]
            choice = ChoiceTuple(m, 0.0, 0.0) #consuming all
            Cstar[state_index] = choice
        return Cstar

    def solve(self, par):
        # create state_space grid values. order (m, f, p)
        statespace = create_statespace(par)
        m_grid, f_grid, p_grid = create_grids(par)

        V_solution, C_solution = dict(), dict()

        # V, C in period T
        Vstar = self.initialize_Vstar(statespace, par)
        Cstar = self.initialize_Cstar(statespace, par)
        # backwards loop:
        V_solution[par.max_age], C_solution[par.max_age] = Vstar, Cstar

        # 1) (V_star_interpolant) interpolant over næste periode mellem m_grid og v_star_t+1
        for t in reversed(range(par.start_age, par.max_age)):

            print("======== BEGINNING ============")
            print('Solution at time step t: ', t, ', time is: ', datetime.datetime.utcnow())
            tic = time.time()

            V_solution[t], C_solution[t] = Vstar, Cstar

            Vstar_plus, Cstar_plus = self.create_Vstar(statespace), self.create_Cstar(statespace)
            ld = create_lookup_dict(statespace, Vstar)
            #import pdb; pdb.set_trace()
            Vstar_mesh = create_mesh(m_grid, f_grid, p_grid, ld, par)

            interpolant = Interpolate3D(m_grid, f_grid, p_grid, Vstar_mesh)

            for s_ix, s in enumerate(statespace):
                #print(NUMBER_OF_ITERATIONS)
                if s_ix % 100 == 0:
                    print(s_ix,  'time is: ', datetime.datetime.utcnow())
                m, f, p = s[0], s[1], s[2]
                state = StateTuple(m, f, p, t)
                V, C = self.find_V_for_choices(state, par, interpolant)
                Vstar_plus[s_ix], Cstar_plus[s_ix] = V, C

            Vstar, Cstar = Vstar_plus, Cstar_plus
            V_solution[t], C_solution[t] = Vstar, Cstar
            print("======== ENDING ===========")
            toc = time.time()
            print(f' t = {t} solved in {toc-tic:.1f} secs')

        return V_solution, C_solution
