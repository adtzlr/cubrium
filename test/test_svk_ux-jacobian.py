# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 20:51:05 2021

@author: adutz
"""

import matplotlib.pyplot as plt
import numpy as np

import pytest

import cubrium


def test_svk_ux_k2():
    ux = cubrium.loadcase.uniaxial
    base_svk(2, 1, 5000, 20, ux, 0)


def base_svk(k, mu, bulk, steps, lcase, ix):
    MDL = cubrium.init()

    MDL.GLO.constitution.matid = 1
    MDL.GLO.constitution.parameters = [mu, bulk, k]
    
    MDL.GLO.constitution.matid = 3
    MDL.GLO.constitution.parameters = [0.5,0.0,0.0,-0.05,0.005, bulk]

    # import loadcase settings
    MDL = lcase(MDL)

    # Update Problem Definition with loadcase settings
    MDL = cubrium.update(MDL)

    # initial solution
    x0 = np.zeros(9)
    lpf0 = 0.0

    # numeric continuation
    Res = cubrium.solve(MDL)(
        x0=x0,
        lpf0=lpf0,
        dxmax=0.05,
        dlpfmax=0.05,
        control0=10,
        jacmode=3,
        jaceps=1e-4,
        maxsteps=steps,
        maxcycles=2,
        maxiter=20,
        tol=1e-10,
    )

    # compile results
    Y = np.array([res.x for res in Res])
    # F = np.array([res.fun for res in Res])
    # J = np.array([np.linalg.det(res.x[:-1].reshape(3,3)+np.eye(3)) for res in Res])

    # recover internal quantities
    history = cubrium.recover(Y, MDL)

    # plot results
    #plt.figure()
    plt.plot(Y[:, ix], Y[:, -1], "-")
    plt.plot(Y[0, ix], Y[0, -1], "ko")

    comp = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]

    plt.xlabel("comp. %d%d of displacement gradient" % comp[ix])
    plt.ylabel("load-proportionality-factor LPF")

    return Res, MDL


if __name__ == "__main__":
    ux = cubrium.loadcase.uniaxial
    bx = cubrium.loadcase.biaxial

    mu = 1
    bulk = 5000
    steps = 40
    
    Res, MDL = base_svk(2, mu, bulk, steps, bx, 0)
    #Res, MDL = base_svk(2, mu, bulk, steps, ux, 0)
    Y = np.array([res.x for res in Res])
    history = cubrium.recover(Y, MDL)
    
    stiffness_shear = [res.jac[3,3] for res in Res[1:]]