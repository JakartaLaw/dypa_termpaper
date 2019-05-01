import numpy as np

from modules.consumptionpreference import create_consumption_preference_dict
from modules.agepolynomial import create_age_poly_dict



class Agent():

    def __init__(self, par, state, education_lvl):


        assert education_lvl in ['<HS', 'HS', 'College']
        # create consumption preference
        par.n = create_consumption_preference_dict(education_lvl, par.start_age, par.max_age)
        par.age_poly = create_age_poly_dict(education_lvl)

        #instantiation
        self.par = par
        self.state = state
        self.e = education_lvl


    # Utility function, CRRA with equivalence scale adjuting
    def utility(self, c, t):
        # not conditioning on education level differences yet
        u = self.par.n[t] * ( (c/self.par.n[t])**(1-self.par.rho_u) ) /(1-self.par.rho_u)
        return u

    # Updating state variables
    def update_f(self, i):
        self._f_old = copy(self.state.f)
        self.state.f = (1-self.par.delta) * self.state.f + i

    def update_m(self, kappa):
        self.state.m = R_tilde(kappa) * self.a() + self.Y_next() + self.tr()

    def update_mu(self):
        self.state.mu = self.par.rho * self.state.mu + self.par.sigma_psi * np.random.normal()

    def reset_f(self):
        self.state.f = copy(self_f_old)

    def a(self, c, i, kappa, t):
        return self.state.m - c - self.pi(i) - self.kappa_cost(kappa) + self.par.tr[t]

    def Y_next(self, t):

        if t<65:
            return self.state.mu + self.par.g[t + 1] + self.par.sigma_xi * np.random.normal()
        else:
            return 0

    # Evolution of asset
    ## Remember to add stochastic element to the return factor evolution, p. 446
    def R_tilde(self, kappa, shock):
        return (1-kappa)*self.R_riskfree() + kappa*self.R_riskful(shock)

    def R_riskful(self, shock):
        return self.par.r_min + self.r() + shock

    def R_riskfree(self):
        # Might be changed don't know if r_min is correct
        return self.par.r_min

    # Government transfer
    def tr(self, x):
        return max(self.par.cmin - x, 0)

    # Return function of financial literacy assumed linear (p. 450)
    def r(self):
        # Calculate slope
        slope = (self.par.r_max - self.par.r_min)/(self.par.f_max - self.par.f_min)
        return self.par.r_min + slope * self.state.f

    # Cost function for financial knowledge (p. 451, bottom). Combining both fixed and variable
    @staticmethod
    def pi(i):
        # constants are derived from the paper
        return 50*(i**1.75)

    @staticmethod
    def kappa_cost(kappa):
        # constants are derived from the paper
        return 750 * (kappa > 0)
