# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 08:37:32 2021

@author: Andreas
"""

import numpy as np
import meshio


def xdmf(history, filename="timeseries"):

    p = np.zeros(3)
    dX = history[0].GLO.cube.edges

    X = np.array(
        [
            p,
            p + dX[:, 0],
            p + dX[:, 0] + dX[:, 1],
            p + dX[:, 1],
            p + dX[:, 2],
            p + dX[:, 0] + dX[:, 2],
            p + dX[:, 0] + dX[:, 1] + dX[:, 2],
            p + dX[:, 1] + dX[:, 2],
        ]
    )

    pts = np.array(
        [
            p + dX[:, 1] / 2 + dX[:, 2] / 2,
            p + dX[:, 2] / 2 + dX[:, 0] / 2,
            p + dX[:, 0] / 2 + dX[:, 1] / 2,
            #
            p + dX[:, 1] / 2 + dX[:, 2] / 2 + dX[:, 0],
            p + dX[:, 2] / 2 + dX[:, 0] / 2 + dX[:, 1],
            p + dX[:, 0] / 2 + dX[:, 1] / 2 + dX[:, 2],
        ]
    )

    cells = [("hexahedron", np.arange(0, 8).reshape(1, 8))]
    verts = [("vertex", np.arange(0, 6).reshape(6, 1))]

    filename1 = filename + "_cube.xdmf"
    with meshio.xdmf.TimeSeriesWriter(filename1) as writer:
        writer.write_points_cells(X, cells)
        for i, h in enumerate(history):
            F = h.INT.gridvec.components
            u = (F @ X.T).T - X
            s = h.INT.cauchy
            lpf = h.EXT.lpf
            s6 = np.array([s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[1, 2], s[2, 0]])
            writer.write_data(
                i,
                point_data={"Displacement": u},
                cell_data={
                    "Cauchy Stress": s6.reshape(1, 6),
                    "Load-Proportionality-Factor (LPF)": np.array([[lpf]]),
                },
            )

    filename2 = filename + "_points.xdmf"
    with meshio.xdmf.TimeSeriesWriter(filename2) as writer:
        writer.write_points_cells(pts, verts)
        for i, h in enumerate(history):
            F = h.INT.gridvec.components
            u = (F @ pts.T).T - pts

            r = h.INT.force.components
            t = h.INT.traction.components

            rr = np.vstack((-r.T, r.T))
            tt = np.vstack((-t.T, t.T))
            writer.write_data(
                i, point_data={"Displacement": u, "Reaction Force": rr, "Traction": tt}
            )

    return
