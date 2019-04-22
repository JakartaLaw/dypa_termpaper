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
