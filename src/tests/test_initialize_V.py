import pytest
from model import Model
from parameters import parameters as par
from modules.stategrid import create_statespace

def test_initialize_V():

    s = create_statespace(par)
    m = Model(par, 'HS')
    init_V = m.initialize_Vstar(s)
    print(init_V)
    print(init_V.shape)
