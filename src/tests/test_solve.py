import pytest

from model import Model
from parameters import parameters as par
from modules.utils import Struct
from collections import namedtuple

def test_can_call_find_V():
    m = Model()
    m.solve(par)
