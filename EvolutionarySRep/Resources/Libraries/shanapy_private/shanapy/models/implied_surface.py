import os
import numpy as np
import vtk
from shanapy.models import srep_io, srep_fitter, functional_maps as fm
from shanapy.models.srep import SrepWrapper
from shanapy import Linker
from shanapy import Interpolater
from shanapy.utils import viz
def implied_boundary_surf(srep=None, tau=1):
    """
    srep: type SrepWrapper, a collection of spokes

    the ordering of spokes is circular, from inner to outer skeleton
    Input tau: time point from [0, 1]; tau=0 means the skeleton; tau=1 means the boundary
    """
    vtk_srep = srep.data
    nrows, ncols = srep.rows, srep.cols
    surface_mesh = vtk.vtkPolyData()
    radial_steps = nrows // 2 + 1
    num_crest_points = nrows * 2 + (ncols - 2) * 2

    interpolate_level = 1
    interp = Interpolater()
    interpolated = interp.interpolate(vtk_srep, interpolate_level, nrows, ncols)
    scaled_spokes = [s.scale(tau) for s in interpolated]
    surface_mesh_visual, skeleton_visual = interp.construct_interp_implied_surface(vtk_srep, scaled_spokes, nrows, ncols, interpolate_level)
    implied_landmarks = vtk.vtkPolyData()# [s.getB() for s in scaled_spokes]
    implied_landmarks_pts = vtk.vtkPoints()
    for s in scaled_spokes:
        implied_landmarks_pts.InsertNextPoint(s.getB())
    implied_landmarks.SetPoints(implied_landmarks_pts)
    #viz.overlay_polydata(surface_mesh, srep.data)

    return surface_mesh_visual, scaled_spokes, implied_landmarks, skeleton_visual

def compute_deformation(srep, target_surf, neibor_surf, link_struct, tau):
    """
    Use thin plate spline algorithm to estimate the deformation
    of the whole wave front according to a few linked spokes

    tau = [0, 1]; tau = 0 when it is the boundary; tau = 1 when it is linking axis
    """
    ### 1. compute movement of linked points from link_struct
    linker_source_pts = vtk.vtkPoints()
    linker_target_pts = vtk.vtkPoints()

    linkee_source_pts = vtk.vtkPoints()
    linkee_target_pts = vtk.vtkPoints()

    for link_tuple in link_struct:
        target_spoke, neibor_spoke, target_id, neibor_id = link_tuple
        if np.isinf(target_spoke.ext):
            continue
        linker_source_pts.InsertNextPoint(target_spoke.getB())
        linkee_source_pts.InsertNextPoint(neibor_spoke.getB())

        linker_front_pt = target_spoke.p + (target_spoke.r + tau * target_spoke.ext) * target_spoke.U
        linkee_front_pt = neibor_spoke.p + (neibor_spoke.r + tau * neibor_spoke.ext) * neibor_spoke.U
        linker_target_pts.InsertNextPoint(linker_front_pt)
        linkee_target_pts.InsertNextPoint(linkee_front_pt)
    linker_source_pts.Modified()
    linker_target_pts.Modified()

    ### Interpolate deformation with thin-plate-spline
    linker_tps = vtk.vtkThinPlateSplineTransform()
    linker_tps.SetSourceLandmarks(linker_source_pts)
    linker_tps.SetTargetLandmarks(linker_target_pts)
    linker_tps.SetBasisToR()
    linker_tps.Modified()

    linkee_source_pts.Modified()
    linkee_target_pts.Modified()

    ### Interpolate deformation with thin-plate-spline
    linkee_tps = vtk.vtkThinPlateSplineTransform()
    linkee_tps.SetSourceLandmarks(linkee_source_pts)
    linkee_tps.SetTargetLandmarks(linkee_target_pts)
    linkee_tps.SetBasisToR()
    linkee_tps.Modified()
    
    ### 2. Estimate the deformation via TPS
    ### 3. Apply the deformation to other points on the surface
    deformed_linker_surf = deform_surface_pts(target_surf, linker_tps)
    deformed_linkee_surf = deform_surface_pts(neibor_surf, linkee_tps)

    ### 4. (optional) Change the direction of s-rep for the new surface
    return deformed_linker_surf, deformed_linkee_surf
def deform_surface_pts(surf, tps):
    ### Apply the deformation onto the spokes
    deformed_surf_pts = vtk.vtkPoints()
    # refined_srep is a polydata that collects spokes
    for i in range(surf.GetNumberOfPoints()):
        pt = surf.GetPoint(i)

        new_pt = tps.TransformPoint(pt)

        id0 = deformed_surf_pts.InsertNextPoint(new_pt)

    deformed_surf_pts.Modified()
    surf.SetPoints(deformed_surf_pts)
    surf.Modified()
    return surf

def implied_crest_surf(crest_spokes=None, tau=1):
    """
    Given the crest spokes (interpolated),

    connect spokes' ends with triangles to form a triangle mesh
    """
    pass
def retile_surf(poly_implied_boundary_surf, poly_target_boundary_surf, match_target=False):
    """
    Implied boundary surface (triangle mesh) may have different number of landmarks.
    This method unifies the number of landmarks between implied and target boundaries.

    Moreover, it (optionally) deforms the retiled surface to match the target boundary
    """
    pass
def link_level_surf(target_srep, neibor_srep, link_struct, tau=1):
    """
    Compute level surfaces due to link flow
    
    """
    target_surf, target_srep_interp, implied_landmarks_target, skeleton_visual = implied_boundary_surf(target_srep)
    neibor_surf, neibor_srep_interp, implied_landmarks_neibor, skeleton_visual = implied_boundary_surf(neibor_srep)

    ## have redundent points for visualization
    target_deform_vis, neibor_deform_vis = compute_deformation(target_srep, target_surf, neibor_surf, link_struct, tau)

    ## no redundent points
    target_deform, neibor_deform = compute_deformation(target_srep, implied_landmarks_target, implied_landmarks_neibor, link_struct, tau)
    #viz.overlay_polydata(target_deform, highlight=neibor_deform)
    target_skeleton_pts = vtk.vtkPoints()
    target_skeleton_poly = vtk.vtkPolyData()
    for i in range(len(target_srep_interp)):
        target_skeleton_pts.InsertNextPoint(target_srep_interp[i].p)
    target_skeleton_poly.SetPoints(target_skeleton_pts)
    target_skeleton_poly.Modified()
    target_skeleton_2_linking_vectors = viz.form_spokes_from_strata(target_skeleton_poly, target_deform)
#    neibor_srep_new = viz.form_spokes_from_strata(corresponding_skeleton_neibor, neibor_deform)
    # target_srep_new = viz.form_spokes_poly(vectors_target)
    # neibor_srep_new = viz.form_spokes_poly(vectors_neibor)
    return target_deform, neibor_deform, target_skeleton_2_linking_vectors, target_srep_interp, target_deform_vis, skeleton_visual

# case_id = '138494'

# target_surf_file = '/playpen/workspace/my_paper/linking/data/hipp_init_srep/' + case_id + '/hipp.vtk'
# nbr_surf_file = '/playpen/workspace/my_paper/linking/data/caud_init_srep/' + case_id + '/caud.vtk'
# reader = vtk.vtkPolyDataReader()
# reader.SetFileName(target_surf_file)
# reader.Update()
# target_mesh = reader.GetOutput()

# caud_mesh_reader = vtk.vtkPolyDataReader()
# caud_mesh_reader.SetFileName(nbr_surf_file)
# caud_mesh_reader.Update()
# nbr_mesh = caud_mesh_reader.GetOutput()

# hipp_srep_file = '/playpen/workspace/my_paper/linking/data/hipp_init_srep/' + case_id + '/header.xml'
# caud_srep_file = '/playpen/workspace/my_paper/linking/data/caud_init_srep/' + case_id + '/header.xml'
# hipp_up_renderable, hipp_down_renderable, hipp_crest_renderable, num_rows, num_cols = \
#                                                                             srep_io.readSrepFromXML(hipp_srep_file)
# caud_up_renderable, caud_down_renderable, caud_crest_renderable, num_rows, num_cols =\
#                                                                                       srep_io.readSrepFromXML(caud_srep_file)
# refined_hipp_down_poly = srep_fitter.refine_srep(hipp_down_renderable.data, target_mesh)
# refined_hipp_down_renderable = SrepWrapper(refined_hipp_down_poly, hipp_down_renderable.rows, hipp_down_renderable.cols)
# refined_caud_up_poly = srep_fitter.refine_srep(caud_up_renderable.data, nbr_mesh)
# refined_caud_up_renderable = SrepWrapper(refined_caud_up_poly, caud_up_renderable.rows, caud_up_renderable.cols)
# linker = Linker(refined_hipp_down_renderable, refined_caud_up_renderable)
# link_struct = linker.link()
# target_deform, neibor_deform, skeleton_2_linking, target_interp_spokes = link_level_surf(refined_hipp_down_renderable, refined_caud_up_renderable, link_struct, tau=1)
# hipp_with_link_length, link_indicator, length_linked, length_unlinked = fm.heatmap_to_object_surface(target_mesh, target_interp_spokes, skeleton_2_linking)
# appender = vtk.vtkAppendPolyData()
# appender.AddInputData(target_mesh)
# appender.AddInputData(skeleton_2_linking)
# appender.AddInputData(nbr_mesh)
# appender.Update()
# import pyvista
# plotter = pyvista.Plotter()
# plotter.add_mesh(target_deform)
# plotter.add_mesh(appender.GetOutput())
# plotter.show()

#viz.overlay_polydata(target_deform, appender.GetOutput())

#     #implied_boundary_surf(srep=hipp_down_renderable, tau=1.5)
# for step in np.arange(0, 1.2, 0.2):
#     target_deform, neibor_deform, _, _ = link_level_surf(hipp_down_renderable, caud_up_renderable, link_struct, tau=step)

#     appender = vtk.vtkAppendPolyData()
#     appender.AddInputData(target_mesh)
#     appender.AddInputData(nbr_mesh)
#     appender.Update()

#     viz.overlay_polydata(target_deform, appender.GetOutput(), highlight=neibor_deform)
# print('done')