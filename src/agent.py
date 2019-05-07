import numpy as np
from copy import copy
from numba import jit

from modules.namedtuples import StateTuple

# Utility function, CRRA with equivalence scale adjuting

@jit
def utility(choice, state, par):
    # not conditioning on education level differences yet
    u = par.n[state.t] * ( (choice.c/par.n[state.t])**(1-par.rho_u) ) /(1-par.rho_u)
    return u

# Updating state variables
@jit
def update_f(choice, state, par):
    """for given choice returns new state (with updated f)"""
    f = (1-par.delta) * state.f + choice.i
    return StateTuple(f=f, p=state.p, m=state.m, t=state.t)

@jit
def calc_a(choice, state, par):
    return state.m - choice.c - pi(choice) - kappa_cost(choice) + tr(state)

# Evolution of asset
## Remember to add stochastic element to the return factor evolution, p. 446
@jit
def R_tilde(choice, state, par, shock):
    return (1-choice.kappa)*R_riskfree(par) + choice.kappa*R_riskful(state, par, shock)

@jit
def R_riskfree(par):
    # Might be changed don't know if r_min is correct
    return 1 + par.r_min

@jit
def R_riskful(state, par, shock):
    return np.exp(par.r_min + r(state, par)) * shock

# Return function of financial literacy assumed linear (p. 450)

@jit
def r(state, par):
    # Calculate slope
    slope = (par.r_max - par.r_min)/(par.f_max - par.f_min)
    return par.r_min + slope * state.f

# Government transfer

@jit
def tr(state):
    return 0
    #return max(par.cmin - x, 0)

# Cost function for financial knowledge (p. 451, bottom). Combining both fixed and variable

@jit
def pi(choice):
    # constants are derived from the paper
    return 50*(choice.i**1.75)


@jit
def kappa_cost(choice):
    # constants are derived from the paper
    return 750 * (choice.kappa > 0)
