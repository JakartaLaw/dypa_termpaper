import pytest

from model import Model
from parameters import parameters as par
from modules.utils import Struct

@pytest.fixture
def dummy_model():
    state = Struct()
    return Model(par=par, state=state, education_lvl='HS')

def test_can_call_model():
    state = Struct()
    model = Model(par=par, state=state, education_lvl='HS')

def test_can_interpolate(dummy_model):
    f = dummy_model.create_interp([1, 2, 3], [5,6,7])
    assert f(2.5) == pytest.approx(6.5)

def test_initalize_V_star():
    state = Struct()
    model = Model(par=par, state=state, education_lvl='HS')
    Vstar = model.initialize_Vstar()

    for i in Vstar[90]:
        pass
