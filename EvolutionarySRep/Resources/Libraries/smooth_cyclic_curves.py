"""
A few tools for dealing with cyclic curves as points and PolyDatas
"""

import numpy as np
import pyvista as pv
import vtk
from scipy import interpolate

def smooth_curve(points: np.ndarray)-> pv.PolyData:
    points = np.unique(points, axis=0)

    tck, u = interpolate.splprep(points.T, s=0, per=1)
    xi, yi, zi = interpolate.splev(np.linspace(0,1,1000), tck)
    points = np.array(list(zip(xi,yi,zi)))
    curve = curve_pd_from_points(points)
    return curve

def curve_pd_from_points(points: np.ndarray) -> pv.PolyData:

    # make a polydata out of the max_curv_pts
    pd = vtk.vtkPolyData()
    pts = vtk.vtkPoints()
    lines = vtk.vtkCellArray()

    first_id = None
    id_b = None
    for i in range(points.shape[0] - 1):
        id_s = pts.InsertNextPoint(points[i])
        if first_id is None:
            first_id = id_s
        id_b = pts.InsertNextPoint(points[i+1])

        cell = vtk.vtkLine()
        cell.GetPointIds().SetId(0, id_s)
        cell.GetPointIds().SetId(1, id_b)
        lines.InsertNextCell(cell)
    
    # remaining link to close the curve
    cell = vtk.vtkLine()
    cell.GetPointIds().SetId(0, id_b)
    cell.GetPointIds().SetId(1, first_id)
    lines.InsertNextCell(cell)

    pd.SetPoints(pts)
    pd.SetLines(lines)
    
    return pd

def get_points_from_cyclic_pd(pd: pv.PolyData)-> np.ndarray:
    """returns the points array in order (around the cycle)"""

    lines = pd.lines
    # print(lines)

    print(pd.points.shape)
    print(np.unique(pd.points, axis=0).shape)
    quit()
    return