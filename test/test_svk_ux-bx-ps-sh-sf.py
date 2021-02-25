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

def test_base(k, mu, bulk, steps, lcase, ix):
    MDL = cubrium.system.init()

    MDL.GLO.constitution.matid = 1
    MDL.GLO.constitution.parameters = [mu, bulk, k]

    # import loadcase settings
    MDL = lcase(MDL)

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
    plt.plot(Y[:, ix], Y[:, -1], "-")
    plt.plot(Y[0, ix], Y[0, -1], "ko")
    
    comp = [(1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3)]

    plt.xlabel("comp. %d%d of displacement gradient" % comp[ix])
    plt.ylabel("load-proportionality-factor LPF")

    plt.savefig(fname + ".svg")

    return


if __name__ == "__main__":
    ux = cubrium.loadcase.uniaxial
    bx = cubrium.loadcase.biaxial
    ps = cubrium.loadcase.planarshear
    
    sh = cubrium.loadcase.simpleshear
    sf = cubrium.loadcase.simpleshearfree2free3
    
    mu = 1
    bulk = 5000
    steps = 250
    
    kk = [2,1,0,-1]
    
    for k in kk:
        test_base(k,mu,bulk,steps,sf,1)
        test_base(k,mu,bulk,steps,sh,1)
        test_base(k,mu,bulk,steps,ux,0)
        test_base(k,mu,bulk,steps,bx,0)
        test_base(k,mu,bulk,steps,ps,0)