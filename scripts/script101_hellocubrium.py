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


MDL.GLO.constitution.umat = umat_svk
MDL.GLO.constitution.parameters = [1.0, 5000.0]


def uniaxial(MDL):
    MDL.EXT.force.normal[0] = 1
    MDL.EXT.force.normal[1] = 0
    MDL.EXT.force.normal[2] = 0

    MDL.EXT.force.shear[0, 1] = 0
    MDL.EXT.force.shear[1, 2] = 0
    MDL.EXT.force.shear[0, 2] = 0

    MDL.EXT.gridvec.symmetry = [1, 1, 1]

    MDL.GLO.lpftype = 0
    MDL.GLO.title = "Uniaxial"
    return MDL


MDL = uniaxial(MDL)
MDL = cubrium.update(MDL)

x0 = np.zeros(9)
lpf0 = 0.0

Res = cubrium.solve(MDL)(
    x0=x0,
    lpf0=lpf0,
)

Y = np.array([res.x for res in Res])
history = cubrium.recover(Y, MDL)

import matplotlib.pyplot as plt

plt.plot(1 + Y[:, 0], Y[:, -1], ".-")
plt.xlabel("stretch $\lambda_1$")
plt.ylabel("load-proportionality-factor LPF")
plt.savefig(MDL.GLO.title + "_stretch-lpf.svg")

cubrium.writer.xdmf(
    history,
    filename=MDL.GLO.title,
)
