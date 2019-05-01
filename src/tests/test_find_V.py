import pytest

from model import Model
from parameters import parameters as par
from modules.utils import Struct
from collections import namedtuple

choice_struct = namedtuple('choice', ['kappa', 'i'])

def test_can_call_find_V():
    state = Struct()
    m = Model(par=par, state=state, education_lvl='HS')

    choice = choice_struct(0.55, 1)
    t = 35

    state.f = 5.5
    state.m = 100000
    state.p = 500
    par.tr = [0 for _ in range(91)]

    Vstar = m.initialize_Vstar()
    Vstar[t+1] = Vstar[90]
    m.create_V_interp(Vstar, t=t)

    res = m.find_V(choice,t)
    print(res)

def test_can_call_find_V_90():
    state = Struct()
    m = Model(par=par, state=state, education_lvl='HS')

    choice = choice_struct(0.55, 1)
    t = 90

    state.f = 5.5
    state.m = 100000
    state.p = 500
    par.tr = [0 for _ in range(91)]

    res = m.find_V(choice,t)
    print(res)
