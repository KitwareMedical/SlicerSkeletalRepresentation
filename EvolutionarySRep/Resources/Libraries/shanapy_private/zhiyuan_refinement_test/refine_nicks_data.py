"""
This file is testing the refiner implemented in shanapy/models/sreps
Run this script in the root directory, e.g., ~/shanapy/
"""
from shanapy.models.sreps.refiner import Refiner
import vtk
from shanapy.models.sreps import Initializer, refiner
from shanapy.visualization import SrepViewer
import numpy as np
import pyvista as pv
def rotate_vector(vector, angle, axis=np.array([1, 1, 1])):
    """
    Rotate a 3D vector by an angle around an axis using the rotation matrix.

    Args:
        vector (np.ndarray): 3D vector to rotate.
        angle (float): Angle to rotate the vector by in radians.
        axis (np.ndarray): 3D axis to rotate the vector around.

    Returns:
        np.ndarray: Rotated vector.
    """

    # Normalize the rotation axis
    axis = axis / np.linalg.norm(axis)

    # Calculate the rotation matrix
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    rot_matrix = np.array([
        [cos_theta + axis[0]**2 * (1 - cos_theta),
         axis[0] * axis[1] * (1 - cos_theta) - axis[2] * sin_theta,
         axis[0] * axis[2] * (1 - cos_theta) + axis[1] * sin_theta],
        [axis[1] * axis[0] * (1 - cos_theta) + axis[2] * sin_theta,
         cos_theta + axis[1]**2 * (1 - cos_theta),
         axis[1] * axis[2] * (1 - cos_theta) - axis[0] * sin_theta],
        [axis[2] * axis[0] * (1 - cos_theta) - axis[1] * sin_theta,
         axis[2] * axis[1] * (1 - cos_theta) + axis[0] * sin_theta,
         cos_theta + axis[2]**2 * (1 - cos_theta)]
    ])

    # Apply the rotation matrix to the vector
    return np.dot(rot_matrix, vector)


## Read the input surface mesh (produced by SPHARM-PDM)
## Test data is located in ~/shanapy/data/
reader = vtk.vtkPolyDataReader()
reader.SetFileName('data/pe_warped.vtk')
reader.Update()
input_mesh = reader.GetOutput()

## Initialize an s-rep for the input mesh
srep_reader = vtk.vtkPolyDataReader()
srep_reader.SetFileName('data/reordered_srep.vtk')
srep_reader.Update()
srep_poly = srep_reader.GetOutput()
# num_crest_points=24
# initializer = Initializer(num_crest_points)
# srep_poly = initializer.fit(input_mesh)

# for i in range(srep_poly.GetNumberOfPoints() // 2):
#     base_pt = np.array(srep_poly.GetPoint(i * 2))
#     bdry_pt = np.array(srep_poly.GetPoint(i * 2 + 1))

#     plt = pv.Plotter()
#     plt.add_mesh(input_mesh, color='white', opacity=0.2)
#     plt.add_mesh(srep_poly, color='blue', opacity=0.8)
#     plt.add_mesh(pv.Line(base_pt, bdry_pt), color='orange', line_width=4)
#     plt.show()

## Refine the above s-rep
refiner = Refiner(input_mesh)
refined_srep = refiner.refine(srep_poly)

## Collect fold points from refined srep to form a loop (i.e., fold curve)
fold_pts = []
num_pts = refined_srep.GetNumberOfPoints()
for fold_pt_id in range(num_pts - 1*2, num_pts - 25 * 2, -2):
    fold_pts.append(refined_srep.GetPoint(fold_pt_id))
fold_pts.append(fold_pts[0])

## Visualize the s-rep and the input mesh
viewer = SrepViewer()
viewer.srep_with_fold_in_surface(fold_pts, refined_srep, input_mesh, srep_poly)
