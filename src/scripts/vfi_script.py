
# Brute force method of solving (inspired by EGM in Exercise 2, model.py)
def find_V_bruteforce(f_t, m_t, p_t, tr_t, c_t, kappa_t, invest_t, t, par):
    '''Calculates E_t(V_t+1) via brute force looping'''
    V_fut = 0
    for j_psi in range(1, par.Npsi):
        for i_xi in range(1, par.Nxi):
            for k_eps in range(1, par.Neps):
                r = cls.r((1-par.delta) * f_t + invest_t) # Calculate return today
                interest_factor = (1-kappa_t) * par.R_bar + kappa_t * (par.r_bar + r + par.eps[k_eps])
                assets = m_t + tr_t - c_t - cls.pi_cost(invest_t) - cls.kappa_cost(kappa_t)
                income = par.xi[i_xi] * (par.G * p_t *  par.psi[j_psi] + cls.agepoly[t+1])

                integrand = interest_factor * assets + income

                V = par.psi_w * par.xi_w * par.eps_w * par.Vplus_interp(integrand) # GH weighting
                V_fut += V

    return(V_fut)


# Pre-compute integral for each state (to avoid looping 8^3 times for each grid in state space and choice)

# f_t in [0,100]
# m_t in [0, 250,000]
# p_t in [0, X]

# import numpy as np
# from scipy.interpolate import RegularGridInterpolator
# def f(x,y,z):
#     return 2 * x**3 + 3 * y**2 - z
# x = np.linspace(1, 4, 11)
# y = np.linspace(4, 7, 22)
# z = np.linspace(7, 9, 33)
# data = f(*np.meshgrid(x, y, z, indexing='ij', sparse=True))
#
# data
