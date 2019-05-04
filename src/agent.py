import numpy as np
from copy import copy

from modules.namedtuples import StateTuple

# Utility function, CRRA with equivalence scale adjuting
def utility(choice, state, par):
    # not conditioning on education level differences yet
    u = par.n[state.t] * ( (choice.c/par.n[state.t])**(1-par.rho_u) ) /(1-par.rho_u)
    return u

# Updating state variables
def update_f(choice, state, par):
    """for given choice returns new state (with updated f)"""
    f = (1-par.delta) * state.f + choice.i
    return StateTuple(f=f, p=p, m=m, t=t)


def calc_a(choice, state, par):
    return state.m - choice.c - pi(choice.i) - kappa_cost(choice.kappa) + par.tr[state.t]

# Evolution of asset
## Remember to add stochastic element to the return factor evolution, p. 446
def R_tilde(self, choice, state, shock):
    return (1-choice.kappa)*R_riskfree(par) + choice.kappa*R_riskful(state, par, shock)

def R_riskfree(par):
    # Might be changed don't know if r_min is correct
    return par.r_min

def R_riskful(self, state, par, shock):
    return par.r_min + r(state, par) + shock

# Return function of financial literacy assumed linear (p. 450)
def r(state, par):
    # Calculate slope
    slope = (par.r_max - par.r_min)/(par.f_max - par.f_min)
    return par.r_min + slope * state.f

def update_mu(self):
    pass
    #self.state.mu = self.par.rho * self.state.mu + self.par.sigma_psi * np.random.normal()

def reset_f(self):
    pass
    #self.state.f = copy(self._f_old)

def Y_next(self, t):

    if t<65:
        return state.mu + par.g[t + 1] + par.sigma_xi * np.random.normal()
    else:
        return 0




# Government transfer
def tr(self, x):
    return max(par.cmin - x, 0)



# Cost function for financial knowledge (p. 451, bottom). Combining both fixed and variable
@staticmethod
def pi(i):
    # constants are derived from the paper
    return 50*(i**1.75)

@staticmethod
def kappa_cost(kappa):
    # constants are derived from the paper
    return 750 * (kappa > 0)
