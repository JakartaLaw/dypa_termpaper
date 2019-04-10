import pytest

from agent import Agent
from parameters import parameters as par
from modules.utils import Struct

state = Struct()

def test_can_call_agent():
    Agent(par=par, state=state, education_lvl='HS')

# ma an instantiated agent
def test_utility_func1(ma):
    ma.utility(c=10000, t=50)
