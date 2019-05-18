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

def create_consumption_preference_array(e_lvl):
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

    k_et = create_avg_children(e_lvl)
    j_et = create_avg_adults(e_lvl)

    return z(j_et,k_et) / z(2,1)


def create_avg_children(e_lvl):
    '''Data from https://www.census.gov/data/tables/2016/demo/families/cps-2016.html'''

    age = [25,30,35,40,45,50,55,60,65, 90]
    col = [0.70, 1.20, 1.45, 1.29, 0.90, 0.47, 0.23, 0.13, 0.09, 0.05]
    hs = [val * 1.00 for val in col]
    lhs = [val * 1.40 for val in col] # Notice the multiplying factor

    names = {
        '<HS': lhs,
        'HS': hs,
        'College': col
    }

    assert e_lvl in names.keys()
    f = interpolate.interp1d(age, names[e_lvl], kind='linear', fill_value = "extrapolate")

    output = np.array([float(f(_age)) for _age in range(0, 90 + 1)])
    output[0:25] = 0 # For t<25 set size to 0, else there will be problems in the z function
    return output

def create_avg_adults(e_lvl):
    '''Data from https://www.census.gov/data/tables/2016/demo/families/cps-2016.html'''

    age = [25,30,35,40,45,50,55,60,65, 90]
    col = [1.84, 1.82, 1.92, 2.04, 2.21, 2.19, 2.10, 1.94, 1.83, 1.60]
    hs = [val * 1.00 for val in col]
    lhs = [val * 1.00 for val in col]

    names = {
        '<HS': lhs,
        'HS': hs,
        'College': col
    }

    assert e_lvl in names.keys()
    f = interpolate.interp1d(age, names[e_lvl], kind='linear', fill_value = "extrapolate")

    output = np.array([float(f(_age)) for _age in range(0, 90 + 1)])
    output[0:25] = 0 # For t<25 set size to 0, else there will be problems in the z function
    return output


# Scripts from previous version

# def calc_consumption_preference(j_et, k_et):
#     """Calculates n_et for given age.
#
#     Parameters
#     ----------
#     j_et : float
#         the average number of adults in the household by age and education group.
#     k_et : float
#         the average number of children in the household by age and education group.
#
#     Returns
#     -------
#     type
#         float
#
#     """
#
#
#
#     return z(j_et,k_et) / z(2,1)

# def create_consumption_preference_array(e_lvl):
#
#     age = [25,30,35,40,45,50,55,60,65, 90]
#     col = [0.3, 1.56, 1.65, 1.68, 1.68, 1.60, 1.25, 0.8, 0.5, 0.5]
#     hs = [val * 1.00 for val in col]
#     lhs = [val * 1.40 for val in col]
#
#     names = {
#         '<HS': lhs,
#         'HS': hs,
#         'College': col
#     }
#
#     assert e_lvl in names.keys()
#     f = interpolate.interp1d(age, names[e_lvl], kind='linear', fill_value = "extrapolate")
#
#     return np.array([float(f(_age)) for _age in range(0, 90 + 1)])
