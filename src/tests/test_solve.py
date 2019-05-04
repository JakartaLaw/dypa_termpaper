import pytest

from model import Model
from parameters import parameters as par
from modules.utils import Struct
from collections import namedtuple

choice_struct = namedtuple('choice', ['kappa', 'i'])

def test_can_call_find_V():
    state = Struct()
    m = Model(par=par, education_lvl='HS')

    m.solve()
