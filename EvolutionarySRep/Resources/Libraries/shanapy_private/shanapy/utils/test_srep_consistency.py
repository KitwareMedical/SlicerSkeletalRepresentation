import vtk
import numpy as np
from shape_model import *
import os
from viz import *
from spoke import *

#spharm_mesh('../data/label_img/top4_label.nrrd', output_file_folder='../data/spharm_results/')
def form_spoke_poly(np_base_pt, np_bdry_pt):
    target_spoke_poly = vtk.vtkPolyData()
    ts_pts = vtk.vtkPoints()
    vs_line = vtk.vtkCellArray()
    id0 = ts_pts.InsertNextPoint(np_base_pt)
    id1 = ts_pts.InsertNextPoint(np_bdry_pt)

    line = vtk.vtkLine()
    line.GetPointIds().SetId(0, id0)
    line.GetPointIds().SetId(1, id1)
    vs_line.InsertNextCell(line)

    target_spoke_poly.SetPoints(ts_pts)
    target_spoke_poly.SetLines(vs_line)
    target_spoke_poly.Modified()
    return target_spoke_poly

def show_boundary_pts_connections():
#    os.system('cp ../data/spharm_results/top4_label_SPHARM.vtk ../data/final_mesh/')

    reader = vtk.vtkPolyDataReader()
    reader.SetFileName('../data/final_mesh/std_ell_label_SPHARM.vtk')
    reader.Update()
    ell_mesh = reader.GetOutput()
    transformer = vtk.vtkTransform()
    transformer.Translate((0, 0, 25))
    transformer.Update()
    transform_filter = vtk.vtkTransformPolyDataFilter()
    transform_filter.SetInputData(ell_mesh)
    transform_filter.SetTransform(transformer)
    transform_filter.Update()
    ell_mesh = transform_filter.GetOutput()

    reader = vtk.vtkPolyDataReader()
    reader.SetFileName('../data/final_mesh/top4_label_SPHARM.vtk')
    reader.Update()
    target_mesh = reader.GetOutput()
    assert ell_mesh.GetNumberOfPoints() == target_mesh.GetNumberOfPoints(), \
                "Different number of surface points between ellipsoid and target"
    num_pts = ell_mesh.GetNumberOfPoints()

    i = np.random.randint(num_pts)
    pt_ell = np.array(ell_mesh.GetPoint(i))

    print(i)
    pt_4 = np.array(target_mesh.GetPoint(i))

    reader.SetFileName('../data/final_mesh/top0_label_SPHARM.vtk')
    reader.Update()
    mesh0 = reader.GetOutput()
    pt_0 = np.array(mesh0.GetPoint(i))

    background_appender = vtk.vtkAppendPolyData()
    background_appender.AddInputData(target_mesh)
    background_appender.AddInputData(ell_mesh)
    background_appender.AddInputData(mesh0)
    background_appender.Update()
    overlay_polydata(form_spoke_poly(pt_ell, pt_0), background_appender.GetOutput(), form_spoke_poly(pt_ell, pt_4))


def get_spoke_dir(srep, spoke_id):
    base_pt_id = spoke_id * 2
    bdry_pt_id = base_pt_id + 1
    base_pt = np.array(srep.GetPoint(base_pt_id))
    bdry_pt = np.array(srep.GetPoint(bdry_pt_id))
    spoke = Spoke(base_pt=base_pt, bdry_pt=bdry_pt)
    return np.array(spoke.U)
def is_close_dir(dir1, dir2):
    """
    Input unit directions dir1 and dir2

    check if inner product is positive
    """
    return np.inner(dir1, dir2) > 0

def check_directions_spokes():
    data_folder = '../data_diff_mean/s_reps/'
    spoke_id = 8
    template_reader = vtk.vtkPolyDataReader()
    template_reader.SetFileName(data_folder + 'top_srep0.vtk')
    template_reader.Update()
    top_template_srep = template_reader.GetOutput()
    top_temp_dir = get_spoke_dir(top_template_srep, spoke_id)

    top_bend_reader = vtk.vtkPolyDataReader()
    top_bend_reader.SetFileName(data_folder + 'top_joint_srep0.vtk')
    top_bend_reader.Update()
    top_bend_srep = top_bend_reader.GetOutput()
    top_bend_temp_dir = get_spoke_dir(top_bend_srep, spoke_id)

    bot_template_reader = vtk.vtkPolyDataReader()
    bot_template_reader.SetFileName(data_folder + 'bot_srep0.vtk')
    bot_template_reader.Update()
    bot_template_srep = bot_template_reader.GetOutput()
    bot_temp_dir = get_spoke_dir(bot_template_srep, spoke_id)

    for i in range(60):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(data_folder + 'top_srep' + str(i) + '.vtk')
        reader.Update()
        top_srep = reader.GetOutput()
        top_dir = get_spoke_dir(top_srep, spoke_id)
        if not is_close_dir(top_dir, top_temp_dir):
            print("Top spoke is pointing to the opposite for " + str(i) + "th case")

        reader.SetFileName(data_folder + 'bot_srep' + str(i) + '.vtk')
        reader.Update()
        bot_srep = reader.GetOutput()
        bot_dir = get_spoke_dir(bot_srep, spoke_id)
        if not is_close_dir(bot_dir, top_temp_dir):
            print("Bot spoke is pointing to the opposite for " + str(i) + "th case")

        reader.SetFileName(data_folder + 'top_joint_srep' + str(i) + '.vtk')
        reader.Update()
        top_bend_srep = reader.GetOutput()
        top_bend_dir = get_spoke_dir(top_bend_srep, spoke_id)
        if not is_close_dir(top_bend_dir, top_bend_temp_dir):
            print("Top bend spoke is pointing to the opposite for " + str(i) + "th case")
    print('Done')
#check_directions_spokes()