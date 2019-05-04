import pytest
from parameters import parameters as par

from modules.stategrid import create_statespace

def test_create_statespace1():
    statespace = create_statespace(par)
    # its time 3 because it's a vector of three entries
    assert statespace.size == par.Nm * par.Nf * par.Np * 3
    # being sure how many entries we get we want to create the paring for the interpolation
    assert statespace.shape[0] == par.Nm * par.Nf * par.Np
