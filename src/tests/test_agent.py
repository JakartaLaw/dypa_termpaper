import pytest
from agent import utility, update_f, calc_a, R_tilde, T_max_utility
from parameters import parameters

@pytest.fixture
def par():
    return parameters

def test_utility(par):
    u = utility(4, 4, par)
    print(u)

def test_update_f(par):
    f = update_f(1, 14, par)
    print(f)

def test_calc_a(par):
    a = calc_a(1, 2, 3, 4, par)
    print(a)

def test_calc_R_tilde(par):
    r = R_tilde(0.55, 27, par, 2)
    print(r)

def test_t_max_utility(par):
    res1 = utility(4, 4, par)
    res2 = T_max_utility(4, 4, 4, 4, par)
    assert res1==res2
