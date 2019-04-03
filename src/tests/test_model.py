import pytest

from model import Model
from parameters import parameters as par
from modules.utils import Struct

def test_can_call_model():
    state = Struct()
    model = Model(par=par, state=state)
