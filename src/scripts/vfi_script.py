
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


def V_integrate(self, choice, state, par, interpolant):
    '''Calculates E_t(V_t+1) via brute force looping'''

    # global NUMBER_OF_ITERATIONS
    # NUMBER_OF_ITERATIONS = NUMBER_OF_ITERATIONS + 1

    V_fut = 0.0

    # Calculations that can be moved outside the loop
    state = update_f(choice, state, par) # updateting to f_t+1
    assets = calc_a(choice, state, par)

    # refactored loop to return a value and weight from gaus_hermite
    for psi, psi_w in par.psi:
        for xi, xi_w in par.xi:
            for eps, eps_w in par.eps:

                #state = update_f(choice, state, par) # updateting to f_t+1
                interest_factor = R_tilde(choice, state, par, shock=eps)
                #assets = calc_a(choice, state, par)
                income = xi * (par.G * state.p * psi) + par.age_poly[state.t+1]

                # Future state values
                m_fut = interest_factor * assets + income
                p_fut = par.G * state.p * psi
                f_fut = state.f
                s_fut = (m_fut, f_fut, p_fut)

                V = psi_w * xi_w * eps_w * interpolant(s_fut) # GH weighting

                V_fut += V

    return V_fut


def pre_compute_Vfut(t,):
    for s_ix, s in enumerate(statespace): # Run thru post-decision grid (same as standard grid)
        assets, f_fut, p = s[0], s[1], s[2]

        # Calculate integral
        for psi, psi_w in par.psi:
            for xi, xi_w in par.xi:
                for eps, eps_w in par.eps:

                    income = xi * (par.G * p * psi) + par.age_poly[t+1]
                    interest_factor = R_tilde(choice, state, par, shock=eps)





        # Return V(m,f,p) for each grid point

        # Return interpolant
