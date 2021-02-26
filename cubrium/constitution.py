# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 08:33:10 2021

@author: adutz
"""

import numpy as np
import numpy.linalg as la

from .helpers import ddot, dya, cdya, dev


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
    mu, K = parameters[:2]
    gamma = K - 2 / 3 * mu

    C = F.T @ F
    E = 1 / 2 * (C - np.eye(3))
    S = 2 * mu * E + gamma * np.trace(E) * np.eye(3)

    return F @ S


def umat_ksvk(F, parameters):
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


def umat_ksvk_multi(F, parameters):
    return np.sum(np.array([umat_ksvk(F, p) for p in parameters]), 0)


def umat_tod(F, parameters):
    """(U)ser (MAT)erial Function.
    Returns First Piola-Kirchhoff stress tensor for a given
    deformation gradient tensor with a list of material parameters."""

    # expand list of material parameters
    C10, C01, C11, C20, C30, K = parameters[:6]

    J = la.det(F)
    I = np.eye(3)
    C = F.T @ F
    Cu = J ** (-2 / 3) * C

    I1u = np.trace(Cu)
    I2u = (I1u ** 2 - np.trace(Cu @ Cu)) / 2

    W1u = C10 + 2 * C20 * (I1u - 3) + 3 * C30 * (I1u - 3) ** 2 + C11 * (I2u - 3)
    W2u = C01 + C11 * (I1u - 3)

    Su = 2 * W1u * I + 2 * W2u * (I1u * I - Cu)

    p = K * (J - 1)

    S = dev(Su @ Cu) @ la.inv(C) + p * J * la.inv(C)
    return F @ S


def umatdb(matid):
    "Internal umat switcher based on material id."
    if matid == 0:
        return umat_svk
    elif matid == 1:
        return umat_ksvk
    elif matid == 2:
        return umat_ksvk_multi
    elif matid == 3:
        return umat_tod
    elif matid == 4:
        return umat_nh_compr
    else:
        return
