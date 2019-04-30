import pytest
import numpy as np

from modules.multikeylookup import MultiKeyLookup

def test_can_init():
    mkl = MultiKeyLookup()

def test_can_save_and_retrieve_integers():
    choice = ['a', 'b', 'c']

    mkl = MultiKeyLookup()
    mkl.save(1, 2, 3, 4, choice)
    assert mkl.lookup(1, 2, 3, 4) == choice

def test_can_save_multiple_and_retrieve():

    statespace = list()
    g = np.linspace(0,1,10)
    for i in g:
        for j in g:
            statespace.append([i, j, 1, 1])

    mkl = MultiKeyLookup()
    for state in statespace:
        val = ('{state[0]}, {state[1]}')
        mkl.save(state[0], state[1], state[2], state[3], val)

    #checking random state (19)
    state = statespace[19]
    val = mkl.lookup(state[0], state[1], state[2], state[3])
    assert val == ('{state[0]}, {state[1]}')
