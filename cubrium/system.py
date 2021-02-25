# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 08:37:32 2021

@author: adutz
"""

from types import SimpleNamespace
import numpy as np

from copy import deepcopy as copy

from . import constitution


def init(dlpf=0.05, du=0.05):
    """Init problem namespaces: GLObal, INTernal and EXTernal quantities."""

    MDL = SimpleNamespace()

    GLO = SimpleNamespace()
    INT = SimpleNamespace()
    EXT = SimpleNamespace()

    GLO.dof = SimpleNamespace()
    GLO.dof.force = SimpleNamespace()
    GLO.dof.traction = SimpleNamespace()
    GLO.dof.gridvec = SimpleNamespace()

    GLO.dlpf = dlpf
    GLO.du = du
    GLO.lpftype = np.nan

    GLO.cube = SimpleNamespace()
    GLO.cube.edges = np.eye(3)
    GLO.cube.areas = np.eye(3)

    GLO.constitution = SimpleNamespace()
    GLO.constitution.umat = None

    force = SimpleNamespace()
    force.normal = np.nan * np.ones(3)
    force.shear = np.nan * np.ones((3, 3))
    force.components = np.nan * np.ones((3, 3))

    gridvec = SimpleNamespace()
    gridvec.length = np.nan * np.ones(3)
    gridvec.components = np.nan * np.ones((3, 3))
    gridvec.symmetry = np.nan * np.ones(3)

    INT.force = copy(force)
    EXT.force = copy(force)

    INT.traction = copy(force)
    EXT.traction = copy(force)

    INT.gridvec = copy(gridvec)
    EXT.gridvec = copy(gridvec)

    MDL.INT = INT
    MDL.EXT = EXT
    MDL.GLO = GLO
    return MDL


def update(MDL):
    "Update GLObal namespace with corresponding DOFs from external quantities."
    MDL.GLO.dof.force.normal = np.where(~np.isnan(MDL.EXT.force.normal))
    MDL.GLO.dof.force.shear = np.where(~np.isnan(MDL.EXT.force.shear))
    MDL.GLO.dof.force.components = np.where(~np.isnan(MDL.EXT.force.components))
    MDL.GLO.dof.traction.normal = np.where(~np.isnan(MDL.EXT.traction.normal))
    MDL.GLO.dof.traction.shear = np.where(~np.isnan(MDL.EXT.traction.shear))
    MDL.GLO.dof.gridvec.length = np.where(~np.isnan(MDL.EXT.gridvec.length))
    MDL.GLO.dof.gridvec.components = np.where(~np.isnan(MDL.EXT.gridvec.components))
    MDL.GLO.dof.gridvec.symmetry = np.where(~np.isnan(MDL.EXT.gridvec.symmetry))
    MDL.GLO.ndof = (
        len(MDL.GLO.dof.force.normal)
        + len(MDL.GLO.dof.force.shear[0])
        + len(MDL.GLO.dof.force.components[0])
        + len(MDL.GLO.dof.traction.normal)
        + len(MDL.GLO.dof.traction.shear[0])
        + len(MDL.GLO.dof.gridvec.length)
        + len(MDL.GLO.dof.gridvec.components[0])
        + len(MDL.GLO.dof.gridvec.symmetry)
    )

    MDL.GLO.areas = (
        np.linalg.det(MDL.GLO.cube.edges) * np.linalg.inv(MDL.GLO.cube.edges).T
    )

    if MDL.GLO.constitution.umat is None:
        MDL.GLO.constitution.umat = constitution.umatdb(MDL.GLO.constitution.matid)

    return MDL
