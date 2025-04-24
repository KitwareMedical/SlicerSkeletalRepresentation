"""
Test the functionalities of the Srep class by visualizing.
"""
import pyvista as pv
from nick_srep import Srep
import os
import numpy as np

# radii for ellipsoid:
radii = (5,3.1,3)
# radii = (0.5317556115855036, 0.18166444118218897, 0.1793808089017737)

my_srep = Srep()
pe = pv.ParametricEllipsoid(*radii)

# my_srep.fit_ellipsoid_radii(radii, num_spine_points=15, spine_distribution=(lambda x: abs(x)), num_steps=4)
my_srep.fit_ellipsoid_radii(radii, num_spine_points=10, spine_distribution=None, num_steps=4)
# print(my_srep.polydata.GetNumberOfPoints())

# my_srep.discretize_spokes(8)

# my_srep.load_polydata_file("/home/nick/dev/deformetrica_hotdog/shootings-with-regularization/outputPE_60/Shooting__GeodesicFlow__srep__tp_10__age_1.00.vtk")
# my_srep.load_polydata_file("/home/nick/dev/deformetrica_hotdog/data/sreps/PE_srep_regular.vtk")

"""
Looking at skeleton points here.
Trying to get the "onion skins" inside the skeleton.
first 10 are the spine, next 20 are the first onion skin, next 20 are the second onion skin, etc...
That continues until you get to the onion skin which is not yet at the fold curve.
Then the next 10 are the spine again, the next 20 are the onion skins again, etc..
Again that continues but you don't get all the way to the fold curve.
Finally, the last 2*(num spine points) - 2 points are the fold curve.

So in total there will be 2*(num_spine_points) + 2*(num_steps - 2)*(2*(num_spine_points) - 2) + 2*(num_spine_points)-2
Those three terms accounting for:
    both sets of spine points
    both sets of all intermediate "onion skins" on the skeleton
    the single set of fold curve points
"""
s = len(my_srep.skeleton_point_indices)
d = my_srep.num_spine_points

# spine_pts = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[:d])
# tau1_pts = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[d:3*d-2])
# tau2_pts = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[3*d-2:5*d-4])

# z = 5*d-4
# spine_pts_2 = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[z:z+d])
# tau1_pts_2 = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[z+d:z+3*d-2])
# tau2_pts_2 = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[z+3*d-2:z+5*d-4])
# spine_pts = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[z+3*d-2:z+5*d-4])

# z = z + 5*d-4
# test = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[z:z + 2*d-2])

# plt.add_points(test, render_points_as_spheres=True, point_size=10, color="white")
# plt.add_points(spine_pts_2, render_points_as_spheres=True, point_size=10, color="white")
# plt.add_points(tau1_pts_2, render_points_as_spheres=True, point_size=10, color="grey")
# plt.add_points(tau2_pts_2, render_points_as_spheres=True, point_size=10, color="green")

# # demonstrates how to change a point
# change_ind = my_srep.skeleton_point_indices[5]
# old_point = my_srep.get_points_from_index_list([change_ind])
# print(old_point)
# newpt = old_point + np.array([[0,0,0.5]])
# my_srep.polydata.GetPoints().SetPoint(change_ind, newpt[0])
# my_srep.polydata.Modified()

# plt.add_points(spine_pts, render_points_as_spheres=True, point_size=10, color="white")
# plt.add_points(spine_pts[:2], render_points_as_spheres=True, point_size=10, color="green")
# plt.add_points(spine_pts[-2:], render_points_as_spheres=True, point_size=10, color="red")

# ##########################################3


# # warp skelly
# for i in range(my_srep.polydata.GetNumberOfPoints()):
#     pt = my_srep.polydata.GetPoint(i)
#     x, y, z = pt
#     newpt = np.array([x, y, 0.1*x**2 + 0.01*y**2 + z])
#     my_srep.polydata.GetPoints().SetPoint(i, newpt)

# warp skelly
for i in range(my_srep.polydata.GetNumberOfPoints()):
    pt = my_srep.polydata.GetPoint(i)
    x, y, z = pt
    if z > 0.1:
        newpt = np.array([x, y, 2 * z])
    else:
        newpt = np.array([x, y, z])
    my_srep.polydata.GetPoints().SetPoint(i, newpt)

for i in range(pe.GetNumberOfPoints()):
    pt = pe.GetPoint(i)
    x, y, z = pt
    if z > 0.1:
        newpt = np.array([x, y, 2 * z])
    else:
        newpt = np.array([x, y, z])
    pe.GetPoints().SetPoint(i, newpt)

# warp skelly
# for i in range(my_srep.polydata.GetNumberOfPoints()):
#     pt = my_srep.polydata.GetPoint(i)
#     x, y, z = pt
#     yp = y*np.cos(x/3) - z*np.sin(x/3)
#     zp = y*np.sin(x/3) + z*np.cos(x/3)
#     newpt = np.array([(1 + np.exp(-abs(x)))*x, y, np.sin(x)/5 + 0.01*y**2 + z])
#     newpt = np.array([(1 + np.exp(-abs(x)))*x, y + np.cos(y)/5, z + np.sin(x)/5])
#     newpt = np.array([1.5*x, yp, zp])
#     my_srep.polydata.GetPoints().SetPoint(i, newpt)

# my_srep.polydata.Modified()

spine_pts = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[:d])
interior_skel_inds = []
n = my_srep.num_spine_points
m = my_srep.num_skeletal_onionskins
for i in range(m-2):
    interior_skel_inds += my_srep.skeleton_point_indices[
        (1+2*(i))*n-(2*i):(1+2*(i+1))*n-(2*(i+1))
    ]
interior_skel_pts = my_srep.get_points_from_index_list(interior_skel_inds)
start = (1+2*(m-2))*n - 2*(m-2)
fold_inds = my_srep.skeleton_point_indices[-(2*n - 2):]
fold_points = my_srep.get_points_from_index_list(fold_inds)

plt = pv.Plotter()
plt.add_axes()
plt.add_mesh(my_srep.polydata, color='orange',line_width=2, opacity=0.9)
plt.add_mesh(pe, opacity=0.2)
plt.add_points(spine_pts, render_points_as_spheres=True, point_size=5, color="blue")
plt.add_points(interior_skel_pts, render_points_as_spheres=True, point_size=5, color="red")
plt.add_points(fold_points, render_points_as_spheres=True, point_size=5, color="cornflowerblue")
plt.show()

my_srep.regularize_skeleton()
for _ in range(5):
    my_srep.orthogonalize_spokes_to_boundary(pe)
my_srep.refine_spoke_lengths(input_mesh=pe)
my_srep.medialize_skeleton()

spine_pts = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[:d])
interior_skel_inds = []
n = my_srep.num_spine_points
m = my_srep.num_skeletal_onionskins
for i in range(m-2):
    interior_skel_inds += my_srep.skeleton_point_indices[
        (1+2*(i))*n-(2*i):(1+2*(i+1))*n-(2*(i+1))
    ]
interior_skel_pts = my_srep.get_points_from_index_list(interior_skel_inds)
start = (1+2*(m-2))*n - 2*(m-2)
fold_inds = my_srep.skeleton_point_indices[-(2*n - 2):]
fold_points = my_srep.get_points_from_index_list(fold_inds)

plt = pv.Plotter()
plt.add_axes()
plt.add_mesh(my_srep.polydata, color='orange',line_width=2, opacity=0.9)
plt.add_mesh(pe, opacity=0.2,)
plt.add_points(spine_pts, render_points_as_spheres=True, point_size=5, color="blue")
plt.add_points(interior_skel_pts, render_points_as_spheres=True, point_size=5, color="red")
plt.add_points(fold_points, render_points_as_spheres=True, point_size=5, color="cornflowerblue")
plt.show()