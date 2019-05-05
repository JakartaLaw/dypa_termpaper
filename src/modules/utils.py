import numpy as np
from numpy.polynomial.hermite import hermgauss

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

    return np.array((x, w), dtype=np.float32).T

def create_mortality():
    # Mortality from http://www.bandolier.org.uk/booth/Risk/dyingage.html
    mort_0_14 = 0 * np.ones(15)
    mort_15_24 = 1/((1908+4143)/2) * np.ones(10)
    mort_25_34 = 1/((1215+2488)/2) * np.ones(10)
    mort_35_44 = 1/((663+1106)/2) * np.ones(10)
    mort_45_54 = 1/((279+421)/2) * np.ones(10)
    mort_55_64 = 1/((112+178)/2) * np.ones(10)
    mort_65_74 = 1/((42+65)/2) * np.ones(10)
    mort_75 = 1/((15+21+6+7)/4) * np.ones(15)

    mort = np.append(mort_0_14, np.append(mort_15_24, np.append(mort_25_34, np.append(mort_35_44, np.append(mort_45_54, np.append(mort_55_64, np.append(mort_65_74, mort_75)))))))
    return np.array(mort, dtype=np.float32)
