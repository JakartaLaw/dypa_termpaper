import pytest

from agent import Agent
from parameters import parameters as par
from modules.utils import Struct


state = Struct()

def test_can_call_agent():
    Agent(par=par, state=state)
