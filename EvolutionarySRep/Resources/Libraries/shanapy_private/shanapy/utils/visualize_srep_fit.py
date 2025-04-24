import vtk
from shanapy.utils import viz, obtain_class_labels as ocl
from shanapy.models import srep_io, srep_fitter as sf
import os
import numpy as np
import pyvista as pv
def remove_data_array(polydata):
    while polydata.GetPointData().GetNumberOfArrays() > 0:
        polydata.GetPointData().RemoveArray(0)
    return polydata
def compute_vtk_skeleton(spoke_vtk_file, srep_vtk= None):
    if srep_vtk is None:
        up_renderable, down_renderable, crest_renderable, num_rows, num_cols = \
                                                                                    srep_io.readSrepFromXML(spoke_vtk_file)
        up_poly = remove_data_array(up_renderable.data)
        down_poly = remove_data_array(down_renderable.data)
        crest_poly = remove_data_array(crest_renderable.data)
        appender = vtk.vtkAppendPolyData()
        appender.AddInputData(up_poly)
        appender.AddInputData(down_poly)
        appender.AddInputData(crest_poly)
        appender.Update()

        spokes_vtk = appender.GetOutput()
    else:
        spokes_vtk = srep_vtk
    num_rows = 5
    num_cols = 9

    skeleton_vtk = vtk.vtkPolyData()
    skeletal_pts = vtk.vtkPoints()
    skeletal_mesh = vtk.vtkCellArray()
    radii_arr = vtk.vtkDoubleArray()
    radii_arr.SetNumberOfComponents(1)
    radii_arr.SetName("Spoke radii (mm)")
    spokes_radii_arr = vtk.vtkDoubleArray()
    spokes_radii_arr.SetNumberOfComponents(1)
    spokes_radii_arr.SetName("Spoke radii (mm)")
    
    for i in range(0, spokes_vtk.GetNumberOfPoints(), 2):
        skel_pt = np.array(spokes_vtk.GetPoint(i))
        bdry_pt = np.array(spokes_vtk.GetPoint(i+1))
        id_pt = skeletal_pts.InsertNextPoint(skel_pt)
        length = [np.linalg.norm(bdry_pt - skel_pt)]
        radii_arr.InsertNextTuple(length)
        spokes_radii_arr.InsertNextTuple(length)
#        spokes_radii_arr.InsertNextTuple(length)
    skeleton_vtk.SetPoints(skeletal_pts)
    skeleton_vtk.GetPointData().AddArray(radii_arr)
    skeleton_vtk.GetPointData().SetActiveScalars("Spoke radii (mm)")

    while (spokes_vtk.GetPointData().GetNumberOfArrays()) > 0:
        arr_name = spokes_vtk.GetPointData().GetArrayName(0)
        spokes_vtk.GetPointData().RemoveArray(arr_name)
    spokes_vtk.GetPointData().AddArray(spokes_radii_arr)
    spokes_vtk.GetPointData().SetActiveScalars("Spoke radii")
    spokes_vtk.Modified()
    num_concen_circles = num_rows // 2
    num_crest_pt = (num_rows - 2) * 2 + num_cols * 2
    ### rows are radial lines
    for i in range(skeleton_vtk.GetNumberOfPoints()):
        quad = vtk.vtkQuad()
        curr_row = i // (num_concen_circles + 1)
        curr_col = i - curr_row * (num_concen_circles + 1)
        if curr_col >= 0 and curr_row >= 0 and curr_row < num_crest_pt - 1 and curr_col < (num_concen_circles):
            quad.GetPointIds().SetId(0, i)
            quad.GetPointIds().SetId(1, i + num_concen_circles + 1)
            quad.GetPointIds().SetId(2, i + num_concen_circles + 2)
            quad.GetPointIds().SetId(3, i + 1)
        elif curr_row == num_crest_pt - 1 and curr_col < num_concen_circles:
            ### connect the last radial line and 1st radial line
            quad.GetPointIds().SetId(0, i)
            quad.GetPointIds().SetId(1, curr_col)
            quad.GetPointIds().SetId(2, curr_col + 1)
            quad.GetPointIds().SetId(3, i + 1)
        skeletal_mesh.InsertNextCell(quad)
    skeleton_vtk.SetPolys(skeletal_mesh)

    return skeleton_vtk, spokes_vtk

def show_mesh_and_srep_data(case_id):
    target_surf_folder = '/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/'  + case_id

    for mesh_file in os.listdir(target_surf_folder):
        if mesh_file.endswith('.vtk') and "stx_noscale_" in mesh_file:
            target_surf_file = os.path.join(target_surf_folder, mesh_file)
            reader = vtk.vtkPolyDataReader()
            reader.SetFileName(target_surf_file)
            reader.Update()
            target_mesh = reader.GetOutput()
            hipp_srep_file = os.path.join(target_surf_folder, 'header.xml')

            skeleton_vtk, spokes_vtk, up_poly = compute_vtk_skeleton(hipp_srep_file)
            # frame_reader = vtk.vtkPolyDataReader()
            # frame_reader.SetFileName('/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/107524/frames.vtk')
            # frame_reader.Update()
            # frames = frame_reader.GetOutput()
            p = pv.Plotter()
            p.add_mesh(skeleton_vtk)
            p.add_mesh(up_poly)
#            p.add_mesh(frames, line_width=5)
        #    p.add_mesh(spokes_vtk)
            p.add_mesh(target_mesh, opacity=0.3, color="white")
            p.show()

            # viz.overlay_polydata(appender.GetOutput(), target_mesh)
            # viz.compare_first_spokes(srep_poly, srep_poly, target_mesh)

if __name__ == "__main__":
    show_mesh_and_srep_data('841812')
    # pos_ids, neg_ids = ocl.load_class_labels()

    # for i, p_id in enumerate(pos_ids):
    #     print(i, p_id)
    #     show_mesh_and_srep_data(p_id)


