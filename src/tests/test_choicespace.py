import pytest

from model import Model
from parameters import parameters as par
from modules.utils import Struct

def test_initalize_V_star():
    state = Struct()
    model = Model(par=par, state=state, education_lvl='HS')
    model.create_choicespace()

    cs = model.choicespace
    choice = cs[3]
    print(choice.kappa, choice.i)
