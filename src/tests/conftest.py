import pytest

from agent import Agent
from model import Model
from parameters import parameters as par
from modules.utils import Struct

@pytest.fixture
def ma():
    #agent for test suite
    state = Struct()
    return Agent(par=par, state=state, education_lvl='HS')

@pytest.fixture
def mf():
    #model for test suite
    state = Struct()
    return Model(par=par, state=state, education_lvl='HS')
