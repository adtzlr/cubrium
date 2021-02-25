# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 08:29:52 2021

@author: adutz
"""

import numpy as np


def gridvecns(F):
    "Normal and shear components of gridvectors (edges of cube)."

    Fn = np.zeros(3)
    Fs = np.zeros((3, 3))
    for i in range(3):
        Fn[i] = np.linalg.norm(F[:, i])

    Fs = F.copy()
    # np.fill_diagonal(Fs,np.nan)

    return Fn, Fs


def defgrd(H):
    "Calculates deformation gradient from displacement gradient tensor."
    return np.eye(3) + H.reshape(3, 3)
