import numpy as np

"""
Module made for the sole purpose of mapping the statespace into proper format for the interpolater
"""


def create_lookup_dict(statespace, Vstar):
    d = dict()
    for i in range(statespace.shape[0]):
        key, value = statespace[i], Vstar[i]
        d[tuple(key)] = value
    return d

def create_mesh(m_grid, f_grid, p_grid, lookup_dict, par):
    res = np.empty(shape=(par.Nm, par.Nf, par.Np))

    for m_ix in range(par.Nm):
        for f_ix in range(par.Nf):
            for p_ix in range(par.Np):
                tup = (m_grid[m_ix], f_grid[f_ix], p_grid[p_ix])
                res[m_ix, f_ix, p_ix] = lookup_dict[tup]

    return res
