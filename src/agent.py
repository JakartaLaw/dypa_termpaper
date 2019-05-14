import numpy as np
from copy import copy
from numba import jit, njit
from numba.types import float64, int64, double

# Utility function, CRRA with equivalence scale adjuting

#@njit
#@njit((double, double, double))
def utility(c, n, rho_u):
    # not conditioning on education level differences yet
    return n * ( (c/n)**(1-rho_u) ) /  (1-rho_u)
#
# @njit
# def T_max_utility(m, f, p, t, par):
#     """specifically made for when instantiating V"""
#     return utility(m, t, par)
#

# Updating state variables
#@njit
def update_f(i, f, par):
    """for given choice returns new state (with updated f)"""
    f = (1-par.delta) * f + i
    return f

#@njit
def calc_a(c, i, kappa, m, par):
    return m - c - pi(i) - kappa_cost(kappa) + tr()

# Evolution of asset
## Remember to add stochastic element to the return factor evolution, p. 446
#@njit
def R_tilde(kappa, f, par, shock):
    return (1-kappa)*R_riskfree(par) + kappa*R_riskful(f, par, shock)

#@njit
def R_riskfree(par):
    # Might be changed don't know if r_min is correct
    return par.R_bar

#@njit
def R_riskful(f, par, shock):
    return np.exp(par.r_min + r(f, par)) * shock

# Return function of financial literacy assumed linear (p. 450)
#@njit
def r(f, par):
    # Calculate slope
    slope = (par.r_max - par.r_min)/(par.f_max - par.f_min)
    return par.r_min + slope * f

# Government transfer
#@njit
def tr():
    return 0
    #return max(par.cmin - x, 0)

# Cost function for financial knowledge (p. 451, bottom). Combining both fixed and variable
#@njit
def pi(i):
    # constants are derived from the paper
    return 50*(i**1.75)

#@njit
def kappa_cost(kappa):
    # constants are derived from the paper
    return 750 * (kappa > 0)
