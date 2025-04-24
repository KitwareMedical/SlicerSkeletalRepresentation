"""
This file defines linking maps that maps every vertex on a shape to a real value denoting the radial distances
"""
import os
import numpy as np
import vtk
from shanapy.models import srep_io, implied_surface, srep_fitter
from shanapy.models.srep import SrepWrapper
from shanapy import Linker
from shanapy import Interpolater
from shanapy.utils import viz
def compute_distance_strata(skeleton, bdry):
    """
    The skeleton and bdry are polydata that have correspondence
    
    """
    distances = []
    for i in range(skeleton.GetNumberOfPoints()):
        tail = np.array(skeleton.GetPoint(i))
        head = np.array(bdry.GetPoint(i))

        distances.append(np.linalg.norm(head-tail))
    return np.array(distances)
def compute_length(connection):
    """
    Input connection (vtkPolyData) connects skeletal points and wavefront points
    This function compute the Euclidean distances between pairs of such points
    """
    lengths = []
    for i in range(connection.GetNumberOfCells()):
        cell = connection.GetCell(i)
        tail = np.array(connection.GetPoint(cell.GetPointId(0)))
        head = np.array(connection.GetPoint(cell.GetPointId(1)))

        diff = head - tail
        lengths.append(np.linalg.norm(diff))
    return np.array(lengths)
def link_map(target_srep, neibor_srep, target_srep_unlinked):
    """
    Map every vertex on mesh1 and mesh2 to a real number denoting the distances from
    (1) the vertex (on boundary) to the linking axis for linked point or
    (2) the vertex (on boundary) to the skeletal structure to the linking axis for unlinked point

    Return a link map for target_srep at different regions
    """
    linker = Linker(target_srep, neibor_srep)
    link_struct = linker.link()
    target_deform, neibor_deform, corr_target_connection, corr_neibor_connection = implied_surface.link_level_surf(target_srep, neibor_srep, link_struct, tau=1)
    link_map_linked_region_target = compute_length(corr_target_connection)

    vtk_link_map = vtk.vtkDoubleArray()
    vtk_link_map.SetNumberOfComponents(1)
    for i in range(len(link_map_linked_region_target)):
        vtk_link_map.InsertNextTuple([link_map_linked_region_target[i]])
    target_deform.GetPointData().SetScalars(vtk_link_map)
    unlinked_map = link_map_unlinked_region(target_srep_unlinked, corr_target_connection)
    return link_map_linked_region_target, target_deform, vtk_link_map
def link_map_unlinked_region(target_srep, vec_skeleton_2_linking_target):
    """
    The link distance of unlinked points is defined as the distance from the boundary
    to the skeletal then to the linking axis

    The input s-reps are spokes of non-linked regions.
    """
    ## 1. interpolate sreps
    interpolate_level = 2
    interp = Interpolater()
    nrows, ncols = target_srep.rows, target_srep.cols
    interpolated = interp.interpolate(target_srep.data, interpolate_level, nrows, ncols)

    ## 2. add up linking lengths
    poly_appender = vtk.vtkAppendPolyData()

    for i in range(vec_skeleton_2_linking_target.GetNumberOfCells()):
        cell = vec_skeleton_2_linking_target.GetCell(i)
        tail = np.array(vec_skeleton_2_linking_target.GetPoint(cell.GetPointId(0)))
        head = np.array(vec_skeleton_2_linking_target.GetPoint(cell.GetPointId(1)))

        link_poly_data = viz.form_spoke_poly(tail, head)
        poly_appender.AddInputData(link_poly_data)
    poly_appender.Update()
def heatmap_to_object_surface(mesh, spokes, link_vectors, spokes_another_side):
    """
    Map the link length from implied boundary given by spokes to the object surface given by mesh

    Basic idea: a boundary point can be in the link region, then the closest point to this boundayr point should be a spoke head. Otherwise, a boundary point in the unlinked region is closer to a spoke tail.
    Input: mesh the object surface mesh
    Input spokes: interpolated linking spokes in the mesh for computing implied boundary. No redundent spokes
    Input link_vectors: the extended links of
    the spokes. Use this to compute the link length from the skeleton to the linking axis
    Input spokes_another_side: the unlinked spokes, if the point on mesh is not linked, find the spoke the spokes_another_side
    """
    spokes_ends_pt = vtk.vtkPoints()
    spokes_ends = vtk.vtkPolyData()
    for s in spokes:
        spokes_ends_pt.InsertNextPoint(s.p)
        spokes_ends_pt.InsertNextPoint(s.getB())
    spokes_ends.SetPoints(spokes_ends_pt)
    spokes_ends.Modified()

    point_locator = vtk.vtkPointLocator()
    result=vtk.vtkIdList()
    point_locator.SetDataSet(spokes_ends)
    link_map_on_mesh = vtk.vtkDoubleArray()
    link_map_on_mesh.SetNumberOfComponents(1)

    linked_or_not = []  ## 1: linked, ## 0: unlinked
    link_length_in_linked = []
    link_length_in_unlinked = []

    for i in range(mesh.GetNumberOfPoints()):
        curr_pt = mesh.GetPoint(i)
        id_closest = point_locator.FindClosestPoint(curr_pt)

        # link_tail = np.array(link_vectors.GetPoint(link_cell.GetPointId(0)))
        # link_head = np.array(link_vectors.GetPoint(link_cell.GetPointId(1)))
        if id_closest % 2 != 0:
            ## The closest point is on the boundary, meaning the current point is in the linked region
            link_tail = np.array(link_vectors.GetPoint(id_closest-1))
            link_head = np.array(link_vectors.GetPoint(id_closest))
            link_length = np.linalg.norm(link_head - link_tail)
            link_map_on_mesh.InsertNextTuple([link_length])

            ## statistics for linked/unlinked regions
            linked_or_not.append(1)
            link_length_in_linked.append(link_length)
        else:
            ## The closest point is on the skeleton, the link length is radius of spoke + link length
            unlinked_spokes_bdry_pt = vtk.vtkPoints()
            for s in spokes_another_side:
                unlinked_spokes_bdry_pt.InsertNextPoint(s.getB())
            unlinked_spokes_poly = vtk.vtkPolyData()
            unlinked_spokes_poly.SetPoints(unlinked_spokes_bdry_pt)
            unlinked_spokes_poly.Modified()
            point_locator_unlinked = vtk.vtkPointLocator()
            point_locator_unlinked.SetDataSet(unlinked_spokes_poly)

            id_closest_unlinked = point_locator_unlinked.FindClosestPoint(curr_pt)
            spoke_closest = spokes_another_side[id_closest_unlinked]
            
            link_tail = np.array(link_vectors.GetPoint(id_closest))
            link_head = np.array(link_vectors.GetPoint(id_closest+1))
            link_length = np.linalg.norm(link_head - link_tail)

            ## this is just an approximated radius, not strictly same with spoke radius
            radius = spoke_closest.r  #np.linalg.norm(np.array(curr_pt) - link_tail)
            link_map_on_mesh.InsertNextTuple([link_length + radius])

            ## statistics for linked/unlinked regions
            linked_or_not.append(0)
            link_length_in_unlinked.append(link_length + radius)
        #### For debug checking the closest point
        # line_segment = viz.form_spoke_poly(curr_pt, link_head)
        # line_segment2 = viz.form_spoke_poly(curr_pt, np.array(spokes_ends.GetPoint(id_closest)))
        # viz.overlay_polydata(line_segment, mesh, line_segment2)
        # if i > 5: break
    mesh.GetPointData().SetScalars(link_map_on_mesh)
    return mesh, linked_or_not, link_length_in_linked, link_length_in_unlinked

if __name__ == '__main__':
    case_id = '138494'

    target_surf_file = '/playpen/workspace/my_paper/linking/data/hipp_init_srep/' + case_id + '/hipp.vtk'
    nbr_surf_file = '/playpen/workspace/my_paper/linking/data/caud_init_srep/' + case_id + '/caud.vtk'
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(target_surf_file)
    reader.Update()
    target_mesh = reader.GetOutput()

    caud_mesh_reader = vtk.vtkPolyDataReader()
    caud_mesh_reader.SetFileName(nbr_surf_file)
    caud_mesh_reader.Update()
    nbr_mesh = caud_mesh_reader.GetOutput()

    hipp_srep_file = '/playpen/workspace/my_paper/linking/data/hipp_init_srep/' + case_id + '/header.xml'
    caud_srep_file = '/playpen/workspace/my_paper/linking/data/caud_init_srep/' + case_id + '/header.xml'
    hipp_up_renderable, hipp_down_renderable, hipp_crest_renderable, num_rows, num_cols = \
                                                                                srep_io.readSrepFromXML(hipp_srep_file)
    caud_up_renderable, caud_down_renderable, caud_crest_renderable, num_rows, num_cols =\
                                                                                          srep_io.readSrepFromXML(caud_srep_file)
    refined_hipp_down_poly = srep_fitter.refine_srep(hipp_down_renderable.data, target_mesh)
    refined_hipp_down_renderable = SrepWrapper(refined_hipp_down_poly, hipp_down_renderable.rows, hipp_down_renderable.cols)
    refined_caud_up_poly = srep_fitter.refine_srep(caud_up_renderable.data, nbr_mesh)
    refined_caud_up_renderable = SrepWrapper(refined_caud_up_poly, caud_up_renderable.rows, caud_up_renderable.cols)
    linker = Linker(refined_hipp_down_renderable, refined_caud_up_renderable)
    link_struct = linker.link()
    target_deform, neibor_deform, skeleton_2_linking, target_interp_spokes, \
        target_deform_vis, skeleton_vis = implied_surface.link_level_surf(refined_hipp_down_renderable, refined_caud_up_renderable, link_struct, tau=1)

    refined_hipp_up_poly = srep_fitter.refine_srep(hipp_up_renderable.data, target_mesh)
    refined_hipp_up_renderable = SrepWrapper(refined_caud_up_poly, hipp_up_renderable.rows, hipp_up_renderable.cols)
    interpolate_level = 3
    interp = Interpolater()
    interpolated_hipp_up = interp.interpolate(refined_hipp_up_poly, interpolate_level, hipp_up_renderable.rows, hipp_up_renderable.cols)

    hipp_with_link_length, link_indicator, length_linked, length_unlinked = heatmap_to_object_surface(target_mesh, target_interp_spokes, skeleton_2_linking, interpolated_hipp_up)

    # writer = vtk.vtkPolyDataWriter()
    # writer.SetFileName('test_hipp.vtk')
    # writer.SetInputData(hipp_with_link_length)
    # writer.Update()
    import pyvista
    plotter = pyvista.Plotter()
    plotter.add_mesh(hipp_with_link_length, scalars=hipp_with_link_length.GetPointData().GetScalars())
#    plotter.add_mesh(hipp_with_link_length, scalars=link_indicator)
    plotter.add_scalar_bar(title='Link distance (mm)', italic=True)
    opacity = 0.5
    #plotter.add_mesh(target_mesh, opacity=opacity)
#    plotter.add_mesh(nbr_mesh, opacity=opacity)
    #plotter.add_mesh(skeleton_2_linking, opacity=opacity)
#    plotter.add_mesh(target_deform_vis)#, scalars=compute_distance_strata(skeleton_vis, target_deform_vis))

    # plotter.add_mesh(pyvista.PolyData(caud_crest_renderable.data))
    # plotter.add_mesh(pyvista.PolyData(target_wavefront), scalars=target_link_map,
    #                  show_scalar_bar=True)

    plotter.show()

    print('done')