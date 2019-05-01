import pytest

from modules.consumptionpreference import z, calc_consumption_preference, create_consumption_preference_array

def test_z():
    assert z(1, 1) == pytest.approx(1.4888, abs=1.0e-3)

def test_calc_consumption_preference():
    assert calc_consumption_preference(2,1) == pytest.approx(1)

def test_create_consumption_preference():
    d = create_consumption_preference_array('<HS')
    #from pprint import pprint; pprint(d)
    assert len(d) == 91
