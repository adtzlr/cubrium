# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 08:32:06 2021

@author: adutz
"""

import numpy as np


def force(F, P, dX, dA):
    "Differential force on undeformed differential area elements."
    df = P(F) @ dA

    dfn = np.zeros(3)
    dfs = np.zeros((3, 3))
    for i in range(3):
        dfn[i] = np.dot(dA[:, i], df[:, i]) / np.linalg.norm(dA[:, i])
        for j in range(3):
            if j != i:
                dfs[i, j] = np.dot(dX[:, j], df[:, i]) / np.linalg.norm(dX[:, j])
    return df, dfn, dfs


def cauchy(F, P):
    return 1 / np.linalg.det(F) * P(F) @ F.T


def traction(F, P, dX, dA):
    "Traction vectors on undeformed differential are elements."
    df = P(F) @ dA

    tn = np.zeros(3)
    ts = np.zeros((3, 3))
    for i in range(3):
        tn[i] = np.dot(dA[:, i], df[:, i]) / np.linalg.norm(dA[:, i]) ** 2
        for j in range(3):
            if j != i:
                ts[i, j] = (
                    np.dot(dX[:, j], df[:, i])
                    / np.linalg.norm(dA[:, i])
                    / np.linalg.norm(dX[:, j])
                )
    return df / np.linalg.norm(dA, axis=0), tn, ts
