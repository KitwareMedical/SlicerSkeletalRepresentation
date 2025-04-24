"""
Test fitting an srep to an ellipsoid.
"""
from re import L
import pyvista as pv
import numpy as np
from srep_fitter import fit_ellipsoid_nick
from scipy.stats import gaussian_kde
from matplotlib import pyplot as plt

# points = np.random.randn(100000)
# points = np.concatenate([np.cos(np.random.randn(1000)), -np.cos(np.random.randn(1000))])
# dist = gaussian_kde(points[np.abs(points) <= 1])
# dist2 = lambda x: (dist(x) + dist(2 - x) + dist(-2 - x)) * (x <= 1) * (x >= -1)

dist = lambda x: x**2 * 1.5

pos = np.linspace(-1,1,500)
plt.plot(pos, dist(pos))
plt.title("Distribution of points on spine")
plt.show()

radii = (3,2,1)

ellipsoid = pv.ParametricEllipsoid(*radii)
srep, coc_indices, vertex_indices = fit_ellipsoid_nick(radii, num_spine_points=100, spine_distribution=dist, spine_distribution_is_from_kde=False)

coc_points = [srep.GetPoint(coc_indices[0]), srep.GetPoint(coc_indices[1])]
coc_points = np.array(coc_points)
vertex_points = [srep.GetPoint(vertex_indices[0]), srep.GetPoint(vertex_indices[1])]
vertex_points = np.array(vertex_points)

# srep = fit_ellipsoid_nick(radii, num_spine_points=10)

p=pv.Plotter()
p.add_mesh(srep, color='blue')
p.add_mesh(ellipsoid, opacity=0.2, color="white")
p.add_points(coc_points, render_points_as_spheres=True, point_size=10.0, color='white')
p.add_points(vertex_points, render_points_as_spheres=True, point_size=10.0, color='red')
p.show()
