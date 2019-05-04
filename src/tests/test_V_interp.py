import pytest
from model import Model
from parameters import parameters as par
from modules.stategrid import create_statespace

def test_V_interp():

    s = create_statespace(par)
    m = Model(par, 'HS')
    init_V = m.initialize_Vstar(s)
    interp = m.create_interp(s, init_V)
    assert interp(s[100]) ==  init_V[100]
