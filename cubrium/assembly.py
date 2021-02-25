# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 08:39:06 2021

@author: adutz
"""

import numpy as np
from functools import partial

from copy import deepcopy as copy

from . import kinematics
from . import kinetics


def recover(Y, MDL):
    "Recover equilibrium of multiple solutions."
    return [copy(system(y[:-1], y[-1], MDL)[1]) for y in Y]


def equilibrium(H, lpf, MDL):
    "System equilibrium function returning only residuals."
    return system(H, lpf, MDL)[0]


def system(H, lpf, MDL):
    "Assemble system equilibrium equations."

    # H, lpf = y[:-1], y[-1]
    F = kinematics.defgrd(H)

    Lpf = np.ones(6)
    Lpf[MDL.GLO.lpftype] = lpf
    MDL.EXT.lpf = lpf

    P = partial(MDL.GLO.constitution.umat, parameters=MDL.GLO.constitution.parameters)
    (
        MDL.INT.force.components,
        MDL.INT.force.normal,
        MDL.INT.force.shear,
    ) = kinetics.force(F, P, MDL.GLO.cube.edges, MDL.GLO.cube.areas)

    (
        MDL.INT.traction.components,
        MDL.INT.traction.normal,
        MDL.INT.traction.shear,
    ) = kinetics.traction(F, P, MDL.GLO.cube.edges, MDL.GLO.cube.areas)

    MDL.INT.gridvec.length, MDL.INT.gridvec.components = kinematics.gridvecns(F)
    MDL.INT.gridvec.volumeratio = np.linalg.det(MDL.INT.gridvec.components)

    MDL.INT.cauchy = kinetics.cauchy(F, P)

    res_fn = (
        -MDL.INT.force.normal[MDL.GLO.dof.force.normal]
        + MDL.EXT.force.normal[MDL.GLO.dof.force.normal] * Lpf[0]
    )

    res_fs = (
        -MDL.INT.force.shear[MDL.GLO.dof.force.shear]
        + MDL.EXT.force.shear[MDL.GLO.dof.force.shear] * Lpf[1]
    )

    res_tn = (
        -MDL.INT.traction.normal[MDL.GLO.dof.traction.normal]
        + MDL.EXT.traction.normal[MDL.GLO.dof.traction.normal] * Lpf[2]
    )

    res_ts = (
        -MDL.INT.traction.shear[MDL.GLO.dof.traction.shear]
        + MDL.EXT.traction.shear[MDL.GLO.dof.traction.shear] * Lpf[3]
    )

    res_Fn = (
        -MDL.INT.gridvec.length[MDL.GLO.dof.gridvec.length]
        + MDL.EXT.gridvec.length[MDL.GLO.dof.gridvec.length] * Lpf[4]
    )

    res_Fc = (
        -MDL.INT.gridvec.components[MDL.GLO.dof.gridvec.components]
        + MDL.EXT.gridvec.components[MDL.GLO.dof.gridvec.components] * Lpf[5]
    )

    res_Sy = (
        -MDL.INT.gridvec.components[[0, 1, 2], [1, 2, 0]][MDL.GLO.dof.gridvec.symmetry]
        + MDL.INT.gridvec.components[[1, 2, 0], [0, 1, 2]][MDL.GLO.dof.gridvec.symmetry]
    )

    residuals = np.hstack((res_fn, res_fs, res_tn, res_ts, res_Fn, res_Fc, res_Sy))

    return residuals, MDL
