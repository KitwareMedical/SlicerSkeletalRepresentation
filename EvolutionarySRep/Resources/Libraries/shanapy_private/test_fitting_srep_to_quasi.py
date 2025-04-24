import pyvista as pv
import numpy as np
from fitting_srep_to_quasi import get_perfect_ellipsoid_and_srep_from_quasi
from generate_points_vtk import save_points_vtk

fname = "/home/nick/dev/deformetrica_hotdog/data/hotdogs/deformed_hotdog60_remesh.vtk"

quasi = pv.PolyData(fname)

pe, srep, coc, vertices, crest_indices = get_perfect_ellipsoid_and_srep_from_quasi(quasi)

p = pv.Plotter()
p.add_mesh(pe, color='white', opacity=1, show_edges=True)
crest_pts = [pe.points[i] for i in crest_indices]
print("num crest points:", len(crest_pts))
crest_pts = np.array(crest_pts)
print("------------------------------")
p.add_points(crest_pts, render_points_as_spheres=True, point_size=20., color="red")
# p.show()

# pe.save("data/hotdogs/deformed_hotdogPE_remesh.vtk")
# srep.save("data/sreps/PE_srep_regular.vtk")

# save_points_vtk("data/coc/PE_coc.vtk", coc)
# save_points_vtk("data/vertices/PE_vertices.vtk", vertices)