from refiner_nick import Refiner
from nick_srep import Srep
import pyvista as pv
import numpy as np
from fitting_srep_to_quasi import get_perfect_ellipsoid_and_srep_from_quasi

radii = (5,4,3)
srep = Srep()
srep.fit_ellipsoid_radii(radii, num_spine_points=13, num_steps=4, spine_distribution=lambda x: np.abs(x))
pe = pv.ParametricEllipsoid(*radii)

plt = pv.Plotter()
plt.add_mesh(srep.polydata, color="blue")
plt.add_mesh(pe, opacity=0.5)
plt.add_points(srep.get_points_from_index_list(srep.skeleton_point_indices), render_points_as_spheres=True, color="red")
plt.add_axes()
plt.show()


# get number of crest points

# num = len(srep.crest_point_indices)
# print(num) # 24

refiner = Refiner(pe)
refiner.refine(srep)


plt = pv.Plotter()
plt.add_mesh(srep.polydata, color="blue")
plt.add_mesh(pe, opacity=0.5)
plt.add_points(srep.get_points_from_index_list(srep.skeleton_point_indices), render_points_as_spheres=True, color="red")
plt.add_axes()
plt.show()
