import pytest
from modules.utils import Struct
from model import Model

from collections import namedtuple

choice_struct = namedtuple('choice', ['kappa', 'i'])


from parameters import parameters as par

def test_can_call_v_integrate():
    state = Struct()
    state.f = 5.5
    state.m = 1000
    state.p = 500
    par.tr = [0 for _ in range(91)]

    choice = choice_struct(0.55, 1)

    m = Model(par=par, state=state, education_lvl='HS')

    t = 35

    Vstar = m.initialize_Vstar()
    Vstar[t+1] = Vstar[90]
    m.create_V_interp(Vstar, t=t)
    res = m.V_integrate(c=100, choice=choice, t=t)
    assert res == pytest.approx(11.7, abs=0.1)
