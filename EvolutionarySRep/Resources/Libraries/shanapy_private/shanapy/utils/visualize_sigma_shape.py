import os
import shape_model as sm
import ellipsoid_simulator as es
import vtk
def compute_heatmap(file_path, template_vtk, output_path):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(file_path)
    reader.Update()
    source_vtk = reader.GetOutput()

    dist_filter = vtk.vtkDistancePolyDataFilter()
    dist_filter.SetInputData(0, template_vtk)
    dist_filter.SetInputData(1, source_vtk)
    dist_filter.Update()
    heatmap = dist_filter.GetOutput()

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(output_path)
    writer.SetInputData(heatmap)
    writer.Update()
    
std = 0.1 / 3.0
bend = std
twist = 0.1
output_folder = '../../diff_mean_heatmaps/'
top_center = (0, 0, 2.5)
bot_center = (0, 0, -2.5)
top_reader = vtk.vtkPolyDataReader()
top_reader.SetFileName('../diff_mean_heatmaps/top_mean_label_SPHARM.vtk')
top_reader.Update()
top_vtk = top_reader.GetOutput()

bot_reader = vtk.vtkPolyDataReader()
bot_reader.SetFileName('../diff_mean_heatmaps/bot_mean_label_SPHARM.vtk')
bot_reader.Update()
bot_vtk = bot_reader.GetOutput()

## deform top ellipsoid
img_file_name = os.path.join(output_folder, 'top_sigma_label.nrrd')
es.generate_ellipsoid( \
                        bend=bend,
                        twist=std,
                        center=top_center,
                        output_file_name=img_file_name)

sm.spharm_mesh(img_file_name, output_file_folder=output_folder)
## move the final mesh to a folder
selected_mesh_file_name = os.path.join(output_folder,
                               'top_sigma_label_SPHARM.vtk')
heatmap_vtk_file_name = os.path.join(output_folder, 'top_heatmap_minus_sigma.vtk')
compute_heatmap(selected_mesh_file_name, top_vtk, heatmap_vtk_file_name)
## deform bot ellipsoi
img_file_name = os.path.join(output_folder, 'bot_sigma_label.nrrd')
es.generate_ellipsoid( \
                        bend=bend,
                        twist=twist,
                        center=bot_center,
                        output_file_name=img_file_name)
sm.spharm_mesh(img_file_name, output_file_folder=output_folder)
## move the final mesh to a folder
selected_mesh_file_name = os.path.join(output_folder,
                               'bot_sigma_label_SPHARM.vtk')
heatmap_vtk_file_name = os.path.join(output_folder, 'bot_heatmap_minus_sigma.vtk')
compute_heatmap(selected_mesh_file_name, bot_vtk, heatmap_vtk_file_name)
print('Done')
