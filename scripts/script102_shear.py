# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:58:52 2021

@author: adutz
"""

import cubrium

MDL = cubrium.init()

import numpy as np


def umat_svk(F, parameters):
    """(U)ser (MAT)erial Function.
    Returns First Piola-Kirchhoff stress tensor for a given
    deformation gradient tensor with a list of material parameters."""

    # expand list of material parameters
    mu, K = parameters[:2]
    gamma = K - 2 / 3 * mu

    C = F.T @ F
    E = 1 / 2 * (C - np.eye(3))

    S = 2 * mu * E + gamma * np.trace(E) * np.eye(3)

    return F @ S


#MDL.GLO.constitution.umat = umat_svk
MDL.GLO.constitution.matid = 1
MDL.GLO.constitution.parameters = [1.0, 5000.0, 2]


def shear(MDL):

    MDL.EXT.force.normal[1] = 0.0
    MDL.EXT.force.normal[2] = 0.0
    MDL.EXT.force.shear[0, 1] = 1

    MDL.EXT.gridvec.components[:, 0] = np.array([1, 0, 0])
    MDL.EXT.gridvec.components[:2, 2] = np.array([0, 0])

    MDL.EXT.gridvec.components[2, 1] = 0

    MDL.GLO.lpftype = 1
    MDL.GLO.title = "Shear (no normal force on surfaces 2 and 3)"
    return MDL


MDL = shear(MDL)
MDL = cubrium.update(MDL)

x0 = np.zeros(9)
lpf0 = 0.0

Res = cubrium.solve(MDL)(
    x0=x0,
    lpf0=lpf0,
    dlpfmax=0.02,
    dxmax=0.02,
    rebalance=True,
    maxsteps=50,
    maxiter=8,
    increase=1,
    high=20,
    
)

Y = np.array([res.x for res in Res])
history = cubrium.recover(Y, MDL)

import matplotlib.pyplot as plt

plt.plot(Y[:, 1], Y[:, -1], ".-")
plt.xlabel("amount of shear $F_{12}$")
plt.ylabel("load-proportionality-factor LPF")
plt.savefig(MDL.GLO.title + "_shear-lpf.svg")

cubrium.writer.xdmf(
    history,
    filename=MDL.GLO.title,
)
