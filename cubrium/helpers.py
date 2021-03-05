# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 16:45:59 2021

@author: Andreas
"""

import numpy as np


def ddot(A, B):
    "Returns the double-dot product of two tensors."
    return np.tensordot(A, B, 2)


def dya(A, B):
    "Returns the dyadic product of two tensors."
    return np.tensordot(A, B, 0)


def cdya(A, B):
    "Returns the cross-dyadic product of two tensors."
    return (np.einsum("ij,kl->ikjl", A, B) + np.einsum("ij,kl->ilkj", A, B)) / 2


def dev(A):
    "Returns the deviator of a tensor."
    return A - np.trace(A) / 3 * np.eye(3)
