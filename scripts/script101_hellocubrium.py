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
    rebalance=True,
    maxsteps=200,
    high=10,
    
)

Y = np.array([res.x for res in Res])
history = cubrium.recover(Y, MDL)
J = np.array([h.INT.gridvec.volumeratio for h in history])

a, b = np.where(J <= 0)[0][[0,-1]]

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.figure()
plt.plot([0],[0],'C0o')
plt.plot(Y[:a, 0], Y[:a, -1], "C0-")
plt.plot(Y[a:b, 0], Y[a:b, -1], "C0--")
plt.plot(Y[b:, 0], Y[b:, -1], "C0-")
plt.xlabel("$\lambda_1 - 1$")
plt.ylabel("load-proportionality-factor LPF")
plt.savefig(MDL.GLO.title + "_stretch-normal-lpf.svg")

plt.figure()
plt.plot([0],[0],'C0o')
plt.plot(Y[:a, 4], Y[:a, -1], "C0-")
plt.plot(Y[a:b, 4], Y[a:b, -1], "C0--")
plt.plot(Y[b:, 4], Y[b:, -1], "C0-")
plt.xlabel("$\lambda_2 - 1$")
plt.ylabel("load-proportionality-factor LPF")
plt.savefig(MDL.GLO.title + "_stretch-transv-lpf.svg")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.plot([0],[0],[0],'C0o')
plt.plot(Y[:a, 0], Y[:a, 4], Y[:a, -1], "C0-")
plt.plot(Y[a:b, 0], Y[a:b, 4], Y[a:b, -1], "C0--")
plt.plot(Y[b:, 0], Y[b:, 4], Y[b:, -1], "C0-")
ax.set_xlabel("$\lambda_1 - 1$")
ax.set_ylabel("$\lambda_2 - 1$")
ax.set_zlabel("LPF")
plt.savefig(MDL.GLO.title + "_stretch1-stretch2-lpf.svg")

cubrium.writer.xdmf(
    history,
    filename=MDL.GLO.title,
)
