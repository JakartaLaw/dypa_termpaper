"""
Class and functions for multi-linear interpolation
"""

import numpy as np
from numba import njit, jitclass, prange, boolean, int64, double, void, float64

@njit(int64(int64,int64,double[:],double))
def binary_search(imin,Nx,x,xi):

    # a. checks
    if xi <= x[0]:
        return 0
    elif xi >= x[Nx-2]:
        return Nx-2

    # b. binary search
    half = Nx//2
    while half:
        imid = imin + half
        if x[imid] <= xi:
            imin = imid
        Nx -= half
        half = Nx//2

    return imin

@njit(double(double[:],double[:],double[:],double[:,:,:],double,double,double))
def interp_3d(grid1,grid2,grid3,value,xi1,xi2,xi3):
    """ raw 3D interpolation """

    # a. search in each dimension
    j1 = binary_search(0,grid1.size,grid1,xi1)
    j2 = binary_search(0,grid2.size,grid2,xi2)
    j3 = binary_search(0,grid3.size,grid3,xi3)

    # b. interpolation
    denom = (grid1[j1+1]-grid1[j1])*(grid2[j2+1]-grid2[j2])*(grid3[j3+1]-grid3[j3])
    nom = 0
    for k1 in range(2):
        nom_1 = (grid1[j1+1]-xi1) if k1 == 0 else (xi1-grid1[j1])
        for k2 in range(2):
            nom_2 = (grid2[j2+1]-xi2) if k2 == 0 else (xi2-grid2[j2])
            for k3 in range(2):
                nom_3 = (grid3[j3+1]-xi3) if k3 == 0 else (xi3-grid3[j3])
                nom += nom_1*nom_2*nom_3*value[j1+k1,j2+k2,j3+k3]

    return nom/denom


spec = [
    ('grid1', double[:]),
    ('grid2', double[:]),
    ('grid3', double[:]),
    ('values', double[:, :, :])
]

@jitclass(spec)
class Interpolate3D(object):
    def __init__(self, grid1, grid2, grid3, values):
        self.grid1 = grid1
        self.grid2 = grid2
        self.grid3 = grid3
        self.values = values

    def interpolate(self, x_i, y_i, z_i):
        return interp_3d(self.grid1, self.grid2, self.grid3, self.values, x_i, y_i, z_i)


@njit(double(double[:],double[:],double[:],double[:,:,:],double,double,double,double,double))
def majority_3d(grid1,grid2,grid3,value, low_val, high_val, xi1,xi2,xi3,):
    """Majority when only to classes

    cutoff is the parameter, that separates the two classes
    """

    cutoff = (low_val + high_val)/2
    res = low_val

    # a. search in each dimension
    j1 = binary_search(0,grid1.size,grid1,xi1)
    j2 = binary_search(0,grid2.size,grid2,xi2)
    j3 = binary_search(0,grid3.size,grid3,xi3)

    class_high_count = 0

    for g1 in range(2):
        for g2 in range(2):
            for g3 in range(2):

                v = value[j1+g1,j2+g2,j3+g3]
                if v > cutoff:
                    class_high_count += 1

    if class_high_count >= 4:

        res = high_val

    return res


spec_vote = [
    ('grid1', double[:]),
    ('grid2', double[:]),
    ('grid3', double[:]),
    ('values', double[:, :, :]),
    ('low_val', double),
    ('high_val', double)
]

@jitclass(spec_vote)
class Majority3D(object):
    def __init__(self, grid1, grid2, grid3, values, low_val, high_val):
        self.grid1 = grid1
        self.grid2 = grid2
        self.grid3 = grid3
        self.values = values
        self.low_val = low_val
        self.high_val = high_val

    def interpolate(self, x_i, y_i, z_i):
        return majority_3d(self.grid1, self.grid2, self.grid3, self.values, self.low_val, self.high_val, x_i, y_i, z_i)
