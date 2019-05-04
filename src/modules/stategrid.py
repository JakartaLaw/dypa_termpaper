import numpy as np

def create_m_grid(par):
    # Grid assets
    grid_m_temp = np.linspace(par.m_min, par.m_max**par.m_tuning, par.Nm)
    grid_m = grid_m_temp ** (1/par.m_tuning)
    return grid_m

def create_f_grid(par):
    # Grid financial knowledge
    grid_f = np.linspace(par.f_min, par.f_max, par.Nf)
    return grid_f

def create_p_grid(par):
    grid_p = np.linspace(par.p_min, par.p_max, par.Np)
    return grid_p

def create_statespace(par):
    '''grid over m, f, p'''

    statespace = list()
    m_grid = create_m_grid(par)
    f_grid = create_f_grid(par)
    p_grid = create_p_grid(par)

    for m in m_grid:
        for f in f_grid:
            for p in p_grid:
                state = (m, f, p)
                statespace.append(state)

    return np.array(statespace)
