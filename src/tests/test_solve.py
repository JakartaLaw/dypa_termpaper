import pytest

from model import Model
from parameters import parameters as par
from modules.utils import Struct
from collections import namedtuple

#@pytest.mark.skip
def test_can_call_find_V():
    m = Model()
    res = m.solve(par)
    import pdb; pdb.set_trace()
    res
