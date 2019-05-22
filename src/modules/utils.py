import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy import interpolate


class Struct():

    def __repr__(self):
        return str(vars(self))

    def __init__(self):
        pass


def hermgauss_lognorm(n, sigma):
    """returns np.array of shape (n, 2)

    x in first column, w in second column.

    Easy to iterate over!

    for psi, psi_w in psi:
        <do loop>
    """

    x, w = hermgauss(n)
    x = np.exp(x*np.sqrt(2)*sigma-0.5*sigma**2);
    w = w/np.sqrt(np.pi);

    return np.array((x, w)).T

def create_mortality():
    '''Mortality from http://www.bandolier.org.uk/booth/Risk/dyingage.html'''
    age = [15,25,35,45,55,65,75]
    mortality_rates = [1/((1908+4143)/2),1/((1215+2488)/2), 1/((663+1106)/2), 1/((279+421)/2), 1/((112+178)/2),
           1/((42+65)/2), 1/((15+21+6+7)/4)]

    f = interpolate.interp1d(age, mortality_rates, kind='linear', fill_value = "extrapolate")

    output = np.array([float(f(_age)) for _age in range(0, 90 + 1)])
    output[0:25] = 0 # For t<25 set size to 0
    return output
