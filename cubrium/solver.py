# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:50:00 2021

@author: adutz
"""

from functools import partial
import contique

from .assembly import equilibrium


def solve(MDL):
    return partial(contique.solve, fun=equilibrium, args=(MDL,))
