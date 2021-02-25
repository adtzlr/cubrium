# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 14:28:44 2021

@author: adutz
"""

import numpy as np


def simpleshear(MDL):

    # external shear load
    MDL.EXT.force.shear[0, 1] = 1

    # fix gridvector "1" and "3" to keep initial directions
    MDL.EXT.gridvec.components[:, 0] = np.array([1, 0, 0])
    MDL.EXT.gridvec.components[:, 2] = np.array([0, 0, 1])

    # keep height in drection "2"
    MDL.EXT.gridvec.components[1, 1] = 1

    # fix inactive shear component of gridvector "2"
    MDL.EXT.gridvec.components[2, 1] = 0

    ## lpf acting on external shear forces
    MDL.GLO.lpftype = 1
    MDL.GLO.title = "Simple Shear"
    return MDL


def simpleshearfree3(MDL):

    # no normal force on surface "3"
    MDL.EXT.force.normal[2] = 0

    # external shear load along direction "0" on surface "1"
    MDL.EXT.force.shear[0, 1] = 1

    # fix gridvector "1"
    # fix gridvector "3" to keep initial direction
    MDL.EXT.gridvec.components[:, 0] = np.array([1, 0, 0])
    MDL.EXT.gridvec.components[:2, 2] = np.array([0, 0])

    # keep height in drection "2"
    MDL.EXT.gridvec.components[1, 1] = 1

    # fix inactive shear component of gridvector "2"
    MDL.EXT.gridvec.components[2, 1] = 0

    ## lpf acting on external shear forces
    MDL.GLO.lpftype = 1
    MDL.GLO.title = "Simple Shear (no normal force on surface 3)"
    return MDL


def simpleshearfree2free3(MDL):

    # no normal forces on surfaces "2" and "3"
    MDL.EXT.force.normal[1] = 0
    MDL.EXT.force.normal[2] = 0

    # external shear load along direction "0" on surface "1"
    MDL.EXT.force.shear[0, 1] = 1

    # fix gridvector "1"
    # fix gridvector "3" to keep initial direction
    MDL.EXT.gridvec.components[:, 0] = np.array([1, 0, 0])
    MDL.EXT.gridvec.components[:2, 2] = np.array([0, 0])

    # fix inactive shear component of gridvector "2"
    MDL.EXT.gridvec.components[2, 1] = 0

    ## lpf acting on external shear forces
    MDL.GLO.lpftype = 1
    MDL.GLO.title = "Simple Shear (no normal force on surfaces 2 and 3)"
    return MDL


def uniaxial(MDL):
    MDL.EXT.force.normal[0] = 1
    MDL.EXT.force.normal[1] = 0
    MDL.EXT.force.normal[2] = 0

    MDL.EXT.force.shear[0, 1] = 0
    MDL.EXT.force.shear[1, 2] = 0
    MDL.EXT.force.shear[0, 2] = 0

    MDL.EXT.gridvec.symmetry = [1, 1, 1]

    MDL.GLO.lpftype = 0
    MDL.GLO.title = "Uniaxial"
    return MDL


def planarshear(MDL):
    MDL.EXT.force.normal[0] = 1
    MDL.EXT.force.normal[2] = 0

    MDL.EXT.force.shear[0, 1] = 0
    MDL.EXT.force.shear[1, 2] = 0
    MDL.EXT.force.shear[0, 2] = 0

    MDL.EXT.gridvec.length[1] = 1

    MDL.EXT.gridvec.symmetry = [1, 1, 1]

    MDL.GLO.lpftype = 0
    MDL.GLO.title = "Planar Shear"
    return MDL


def biaxial(MDL):
    MDL.EXT.force.normal[0] = 1
    MDL.EXT.force.normal[1] = 1
    MDL.EXT.force.normal[2] = 0

    MDL.EXT.gridvec.symmetry = [1, 1, 1]

    MDL.EXT.gridvec.components[1, 0] = 0
    MDL.EXT.gridvec.components[2, 0] = 0
    MDL.EXT.gridvec.components[2, 1] = 0

    MDL.GLO.lpftype = 0
    MDL.GLO.title = "Biaxial"
    return MDL
