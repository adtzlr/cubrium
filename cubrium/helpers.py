# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 16:45:59 2021

@author: Andreas
"""

import numpy as np


def ddot(A, B):
    return np.tensordot(A, B, 2)


def dya(A, B):
    return np.tensordot(A, B, 0)


def cdya(A, B):
    return (np.einsum("ij,kl->ikjl", A, B) + np.einsum("ij,kl->ilkj", A, B)) / 2
