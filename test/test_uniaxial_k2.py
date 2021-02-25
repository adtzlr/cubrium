# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 16:50:18 2021

@author: Andreas
"""

import matplotlib.pyplot as plt
import numpy as np

import pytest

import cubrium
import contique


def test_uniaxial_2():
    test_uniaxial(2, 1, 5000, 250)
    return


def test_uniaxial_1():
    test_uniaxial(1, 1, 5000, 250)
    return


def test_uniaxial_1m():
    test_uniaxial(-1, 1, 5000, 250)
    return


def test_uniaxial_0():
    test_uniaxial(0, 1, 5000, 250)
    return


def test_uniaxial(k, mu, bulk, steps):
    MDL = cubrium.system.init()

    MDL.GLO.constitution.matid = 1
    MDL.GLO.constitution.parameters = [mu, bulk, k]

    # import loadcase settings
    MDL = cubrium.loadcase.uniaxial(MDL)

    # Update Problem Definition with loadcase settings
    MDL = cubrium.system.update(MDL)

    # initial solution
    x0 = np.zeros(9)
    lpf0 = 0.0

    # numeric continuation
    Res = contique.solve(
        fun=cubrium.assembly.equilibrium,
        x0=x0,
        lpf0=lpf0,
        args=(MDL,),
        dxmax=0.05,
        dlpfmax=0.05,
        control0=10,
        jacmode=3,
        jaceps=1e-4,
        maxsteps=steps,
        maxcycles=4,
        maxiter=20,
        tol=1e-10,
    )

    # compile results
    Y = np.array([res.x for res in Res])
    # F = np.array([res.fun for res in Res])
    # J = np.array([np.linalg.det(res.x[:-1].reshape(3,3)+np.eye(3)) for res in Res])

    # recover internal quantities
    history = cubrium.assembly.recover(Y, MDL)

    # xdmf time-series writer
    fname = MDL.GLO.title + "_mu={0:g}_K={1:g}_k={2:g}".format(
        *MDL.GLO.constitution.parameters
    )
    cubrium.writer.xdmf(
        history,
        filename=fname,
    )

    # plot results
    plt.figure()
    plt.plot(Y[:, 0], Y[:, -1], "-")
    plt.plot(Y[0, 1], Y[0, -1], "ko")

    plt.xlabel(r"stretch $\lambda_1$")
    plt.ylabel(r"nondimensional normal force $N_1\ /\ (\mu A)$")

    plt.savefig(fname + ".svg")

    return


if __name__ == "__main__":
    test_uniaxial_2()
    test_uniaxial_1()
    test_uniaxial_0()
    test_uniaxial_1m()
