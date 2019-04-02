import numpy as np
from numpy.polynomial.hermite import hermgauss

class Struct():

    def __repr__(self):
        return str(vars(self))

    def __init__(self):
        pass


def hermgauss_lognorm(n, sigma):

    x, w = hermgauss(n)
    x = np.exp(x*np.sqrt(2)*sigma-0.5*sigma**2);
    w = w/np.sqrt(np.pi);

    return x, w
