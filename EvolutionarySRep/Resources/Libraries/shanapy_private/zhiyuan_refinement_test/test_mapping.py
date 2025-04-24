import numpy as np
import pyvista as pv

osrp = pv.PolyData("srep_before_refinement_after_regularization.vtk")
rsrp = pv.PolyData("reordered_srep.vtk")
nsrp = pv.PolyData("srep_ordered_back.vtk")

# skel_pts = np.array(rsrp.points)[::2]
# redundant_r_pts = []
# for i, pt1 in enumerate(skel_pts):
#     for j, pt2 in enumerate(skel_pts):
#         for k, pt3 in enumerate(skel_pts):
#             if (pt1 == pt2).all() and (pt2 == pt3).all() and i != j and i != k and k != j:
#                 redundant_r_pts.append(pt1)

# redundant_r_pts = np.array(redundant_r_pts)

# for i in range(5):
#     plt = pv.Plotter()
#     plt.add_mesh(nsrp, color="blue")
#     plt.add_points(np.array([nsrp.points[::2][i]]), render_points_as_spheres=True, color="red")
#     plt.add_axes()
#     plt.show()


# plt = pv.Plotter()
# plt.add_mesh(osrp, color="blue")
# plt.add_mesh(rsrp, color='green')
# plt.add_mesh(nsrp, color="red")
# plt.add_points(redundant_r_pts, render_points_as_spheres=True, color="red")
# plt.add_axes()
# plt.show()

try:
    assert osrp.GetNumberOfPoints() == nsrp.GetNumberOfPoints()
except:
    print(osrp.GetNumberOfPoints())
    print(rsrp.GetNumberOfPoints())
    print(nsrp.GetNumberOfPoints())
    quit()

tol = 1e-5
for i in range(osrp.GetNumberOfPoints()):
    opt = np.array(osrp.GetPoint(i))
    npt = np.array(nsrp.GetPoint(i))
    try:
        d = np.linalg.norm(opt - npt)
        assert d < tol
    except:
        print(i)
        print(d)
        print(osrp.GetPoint(i))
        print(nsrp.GetPoint(i))
        quit()

