import pyvista as pv
import numpy as np

z_srep = pv.PolyData("reordered_srep.vtk")
n_srep = pv.PolyData("srep_before_refinement_after_regularization.vtk")

z_pts = np.array(z_srep.points[::2][:2])
n_pts = np.array(n_srep.points[::2][:2])

z_pts += np.array([[0.05,0,0]])

plt = pv.Plotter()
plt.add_axes()
plt.add_mesh(z_srep, color='orange',line_width=1, opacity=0.9)
plt.add_points(z_pts, render_points_as_spheres=True, point_size=5, color="orange")
plt.add_mesh(n_srep, color='blue',line_width=1, opacity=0.9)
plt.add_points(n_pts, render_points_as_spheres=True, point_size=5, color="blue")
plt.show()
