import numpy as np
from scipy import interpolate

def z(j, k):
    """Calculates equivalence scale. z(2,1) is the equivalence scale for a household with two adults and one child

    Parameters
    ----------
    j : float
        the number of adults in the household.
    k : float
        the number of children (under 18 years old).

    Returns
    -------
    type
        float

    """
    return (j + 0.7* k)**0.75

def calc_consumption_preference(j_et, k_et):
    """Calculates n_et for given age.

    Parameters
    ----------
    j_et : float
        the average number of adults in the household by age and education group.
    k_et : float
        the average number of children in the household by age and education group.

    Returns
    -------
    type
        float

    """

    return z(j_et,k_et) / z(2,1)

def create_consumption_preference_dict(e_lvl, start_age, end_age):



    age = [25,30,35,40,45,50,55,60,65, 90]
    col = [0.3, 1.56, 1.65, 1.68, 1.68, 1.60, 1.25, 0.5, 0.0, 0.0]
    hs = [val * 1.07 for val in col]
    lhs = [val * 1.15 for val in col]

    names = {
        '<HS': lhs,
        'HS': hs,
        'College': col
    }

    assert e_lvl in names.keys()
    f = interpolate.interp1d(age, names[e_lvl], kind='linear', fill_value = "extrapolate")

    return {_age: float(f(_age)) for _age in range(start_age, end_age + 1)}
