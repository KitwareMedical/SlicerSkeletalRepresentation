import os
import vtk
data_folder = '../../data_diff_mean/'
#data_folder = 'test_data'

top_srep_file = os.path.join(data_folder, 's_reps/top_srep0.vtk')

## target, neighbor s-reps and target spoke
reader = vtk.vtkPolyDataReader()
reader.SetFileName(top_srep_file)
reader.Update()
target_srep = reader.GetOutput()

bot_srep_file_name = os.path.join(data_folder, 's_reps/bot_srep0.vtk')
bot_srep_reader = vtk.vtkPolyDataReader()
bot_srep_reader.SetFileName(bot_srep_file_name)
bot_srep_reader.Update()
neighbor_srep = bot_srep_reader.GetOutput()

top_mesh_file_name = os.path.join(data_folder, 'final_mesh/top0_label_SPHARM.vtk')
mesh_reader = vtk.vtkPolyDataReader()
mesh_reader.SetFileName(top_mesh_file_name)
mesh_reader.Update()
target_mesh = mesh_reader.GetOutput()
bot_mesh_file_name = os.path.join(data_folder, 'final_mesh/bot0_label_SPHARM.vtk')
mesh_reader = vtk.vtkPolyDataReader()
mesh_reader.SetFileName(bot_mesh_file_name)
mesh_reader.Update()
nbr_mesh = mesh_reader.GetOutput()

from shanapy import Linker
linker = Linker(target_srep, neighbor_srep)
linker.link()
linker.show_links(target_mesh, nbr_mesh)
