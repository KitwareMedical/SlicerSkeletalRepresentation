import vtk
import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib
from vtk.util import numpy_support
from shanapy.utils import viz
from shanapy.models import srep_fitter as sf, mean_curvature_flow as mcf
import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
from glob import glob
import pyvista as pv
import vtk
import scipy
import os
import re
from shanapy.utils import obtain_class_labels as ocl
import pyacvd
from glob import glob

# refined_srep = pv.read("/home/nick/dev/libigl/data/mandible/refined_ellipsoid_srep_65_pv.vtk")
# refined_srep = pv.read(f"./data/obj_srep_mandible_{25}.vtk")
refined_srep = pv.read(f"./data/obj_srep_mandible_{60}_spine_terminals_1_using_srep_endpoints.vtk")

num_points = refined_srep.GetNumberOfPoints()

# print(num_points)

# i = 4
i = 217

point = refined_srep.GetPoint(i)
point4 = refined_srep.GetPoint(4)
point289 = refined_srep.GetPoint(289)
point313 = refined_srep.GetPoint(313)
point76 = refined_srep.GetPoint(76)

test_points = [refined_srep.GetPoint(i) for i in range(313, 314, 2)]
test_points.append(point76)
test_points = np.array(test_points)

# point 313 is the vertex corresponding
# to point76
points = np.array([point76, point313])

# point 289 is the vertex (spoke endpoint) 
# corresponding to point4
# points = np.array([point4, point289])

# points = np.array([point4, point76])

p = pv.Plotter()
p.add_mesh(refined_srep, color='red')
p.add_points(points, render_points_as_spheres=True, point_size=30.)
p.show()
