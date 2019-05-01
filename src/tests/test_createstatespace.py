import pytest

from model import Model
from parameters import parameters as par
from modules.utils import Struct

def test_initalize_V_star():
    state = Struct()
    model = Model(par=par, state=state, education_lvl='HS')
    model.create_statespace()

    ss = model.statespace
    state = ss[0]
    print(state.f, state.p, state.m)
