import pyvista as pv
import numpy as np
from fitting_srep_to_quasi import get_perfect_ellipsoid_and_srep_from_quasi
from nick_srep import Srep
from refiner_wrapper_nick import RefinerWrapper

fname = "/home/nick/dev/deformetrica_hotdog/data/hotdogs/deformed_hotdog60_remesh.vtk"

radii = (5,4,3)
srep = Srep()
srep.fit_ellipsoid_radii(radii, num_spine_points=13, num_steps=4, spine_distribution=lambda x: np.abs(x))
pe = pv.ParametricEllipsoid(*radii)


# warp 
for i in range(srep.polydata.GetNumberOfPoints()):
    pt = srep.polydata.GetPoint(i)
    x, y, z = pt
    newpt = np.array([x, y, 0.1*x**2 + 0.01*y**2 + z])
    srep.polydata.GetPoints().SetPoint(i, newpt)
srep.polydata.Modified()

for i in range(pe.GetNumberOfPoints()):
    pt = pe.GetPoint(i)
    x, y, z = pt
    newpt = np.array([x, y, 0.1*x**2 + 0.01*y**2 + z])
    pe.GetPoints().SetPoint(i, newpt)
pe.Modified()

pv.wrap(srep.polydata).save("refinement_data/srep_pre_refinement.vtk")

srep.regularize_skeleton()
refiner = RefinerWrapper(pe)
srep.polydata = refiner.refine(srep)

pv.wrap(srep.polydata).save("refinement_data/srep_pre_medialization.vtk")

# pv.wrap(srep.polydata).save("/tmp/srep_1_refinement_step.vtk")
# pv.wrap(srep.polydata).save("/tmp/refined_srep_new.vtk")
# pv.wrap(pe).save("/tmp/pe_warped.vtk")

srep.medialize_skeleton()

pv.wrap(srep.polydata).save("refinement_data/srep_post_refinement.vtk")

quit()

plt = pv.Plotter()
plt.add_mesh(srep.polydata, color="blue")
plt.add_mesh(pe, opacity=0.5)
plt.add_points(srep.get_points_from_index_list(srep.skeleton_point_indices), render_points_as_spheres=True, color="red")
plt.add_axes()
plt.show()

# srep.refine_spoke_lengths(input_mesh=pe) # taken care of by zhiyuan's regularization
# srep.medialize_skeleton()
# pv.wrap(srep.polydata).save("/tmp/srep_after_regularization.vtk")

# srep2 = Srep()
# srep2.fit_ellipsoid_radii(radii, num_spine_points=13, num_steps=4, spine_distribution=lambda x: np.abs(x))
# srep2.regularize_skeleton()
# refiner = Refiner(pe)
# refiner.refine(srep2)
# srep2.medialize_skeleton()

plt = pv.Plotter()
plt.add_mesh(srep.polydata, color="blue")
# plt.add_mesh(srep2.polydata, color="green")
plt.add_mesh(pe, opacity=0.5)
plt.add_points(srep.get_points_from_index_list(srep.skeleton_point_indices), render_points_as_spheres=True, color="red")
plt.add_axes()
# plt.show()

