import pytest

from model import Model
from parameters import parameters as par
from modules.utils import Struct

@pytest.fixture
def mf():
    #model for test suite
    state = Struct()
    return Model(par=par, state=state, education_lvl='HS')
