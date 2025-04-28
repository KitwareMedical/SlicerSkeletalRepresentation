from srep_fitter import *
from deform_ellipsoids import *

reader = vtk.vtkPolyDataReader()
reader.SetFileName('../data/final_mesh/std_ell_label_SPHARM.vtk')
reader.Update()
mesh = reader.GetOutput()
spoke_poly, rx, ry, rz = fit_ellipsoid(mesh)
refined_srep = refine_srep(spoke_poly, mesh)
overlay_polydata(refined_srep, mesh)
print('done')
