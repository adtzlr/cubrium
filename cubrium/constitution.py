# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 08:33:10 2021

@author: adutz
"""

import numpy as np
import numpy.linalg as la

from .helpers import ddot, dya, cdya


def umat_nh_compr(F, parameters):
    """(U)ser (MAT)erial Function.
    Returns First Piola-Kirchhoff stress tensor for a given
    deformation gradient tensor with a list of material parameters."""

    # expand list of material parameters
    mu, K = parameters[:2]

    J = la.det(F)
    return mu * (F - J * la.inv(F).T) + K * np.log(J) * la.inv(F).T


def umat_svk(F, parameters):
    """(U)ser (MAT)erial Function.
    Returns First Piola-Kirchhoff stress tensor for a given
    deformation gradient tensor with a list of material parameters."""

    # expand list of material parameters
    mu, K, k = parameters[:3]
    gamma = K - 2 / 3 * mu

    C = F.T @ F
    wC, vC = la.eigh(C)
    w = np.sqrt(wC)
    if k == 0:
        Ekp = np.log(wC) / 2
    else:
        Ekp = 1 / k * (w ** k - 1)
    Ck = np.sum(np.array([wa ** k * dya(Na, Na) for wa, Na in zip(w, vC.T)]), 0)
    Ek = np.sum(np.array([Eka * dya(Na, Na) for Eka, Na in zip(Ekp, vC.T)]), 0)

    Sk = 2 * mu * Ek + gamma * np.trace(Ek) * np.eye(3)

    return F @ Sk @ Ck @ la.inv(C)


def umat_svk_multi(F, parameters):
    return np.sum(np.array([umat_svk(F, p) for p in parameters]), 0)


def umatdb(matid):
    "Internal umat switcher based on material id."
    if matid == 0:
        return umat_nh_compr
    elif matid == 1:
        return umat_svk
    elif matid == 2:
        return umat_svk_multi
    else:
        return
