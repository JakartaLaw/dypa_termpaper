

# Brute force method of solving (inspired by EGM in Exercise 2, model.py)
def find_V_bruteforce(f_t, m_t, tr_t, c_t, kappa_t, invest_t, t, par):
    '''Calculates E_t(V_t+1) via brute force looping'''
    integrand = 0
    for j_psi in range(1, par.Npsi):
        for i_xi in range(1, par.Nxi):
            for k_eps in range(1, par.Neps):
                r = cls.r((1-par.delta) * f_t + invest_t) # Calculate return today
                interest_factor = (1-kappa) * par.R_bar + kappa * (par.r_bar + r + par.eps[k_eps])
                assets = m_t + tr_t - c_t - cls.pi_cost(invest_t) - cls.kappa_cost(kappa_t)
                income = par.xi[i_xi] * par.G * par.psi[j_psi] + cls.agepoly[t+1]

                integrand += interest_factor * assets + income
    return(integrand)
