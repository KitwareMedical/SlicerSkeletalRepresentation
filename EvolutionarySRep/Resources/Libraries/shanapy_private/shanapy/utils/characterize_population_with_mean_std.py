import vtk
from shanapy.simulations import ellipsoid_simulator as es
from shanapy.models import shape_model as sm, srep_fitter as sf
import os
import numpy as np
import logging
import shutil
top_ell_center = (0, 0, 2.5)
bot_ell_center = (0, 0, -2.5)

neighbor_object_names = ['top', 'bot']
objects_center = [top_ell_center, bot_ell_center]

def construct_two_standard_ellipsoids(final_std_mesh_name, output_folder=None):

    ### transform standard_ellipsoid to top and bottom positions
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(final_std_mesh_name)
    reader.Update()
    std_ell_mesh = reader.GetOutput()

    ell_meshes = []
    ell_sreps = []

    for i, center in enumerate(objects_center):
        transformer = vtk.vtkTransform()
        transformer.Translate((0, 0, center[2] * 10))
        transformer.Update()
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetInputData(std_ell_mesh)
        transform_filter.SetTransform(transformer)
        transform_filter.Update()
        ell_mesh = transform_filter.GetOutput()
        ell_meshes.append(ell_mesh)

        spoke_poly, rx, ry, rz = sf.fit_ellipsoid(ell_mesh)
        refined_srep = sf.refine_srep(spoke_poly, ell_mesh)
        ell_sreps.append(refined_srep)
        if output_folder is not None:
            writer = vtk.vtkPolyDataWriter()
            writer.SetFileName(os.path.join(output_folder, 'population_variation', neighbor_object_names[i] + '_ell.vtk'))
            writer.SetInputData(ell_mesh)
            writer.Update()

    return ell_meshes, ell_sreps
## deal with mean shape of a population_name
def compute_moment_shape(output_folder, bend_mu, twist_mu, moment_name='mean'):
    std_mesh_path = os.path.join(output_folder, 'final_mesh', 'std_ell_label_SPHARM.vtk')
    ell_meshes, ell_sreps = construct_two_standard_ellipsoids(std_mesh_path)

    variation_modes_folder = 'population_variation'
    output_heatmap_folder = os.path.join(output_folder, variation_modes_folder)
    if not os.path.exists(output_heatmap_folder):
        os.system('mkdir ' + output_heatmap_folder)

    final_mesh_names = []
    ## 1. generate images (responding A's and B's) containing one deformed ellipsoid
    for j, obj_name in enumerate(neighbor_object_names):
        img_file_name = os.path.join(output_heatmap_folder, obj_name \
                                     + '_' + moment_name + '_label.nrrd')
        es.generate_ellipsoid( \
                                bend=bend_mu[j],
                                twist=twist_mu[j],
                                center=objects_center[j],
                                output_file_name=img_file_name)

        ## 2. generate surface mesh via spharm-pdm
        sm.spharm_mesh(img_file_name, output_file_folder=output_heatmap_folder)
        ## move the final mesh to a folder
        selected_mesh_file_name = os.path.join(output_heatmap_folder,
                                       obj_name + '_' + moment_name + '_label_SPHARM.vtk')
        final_mesh_names.append(selected_mesh_file_name)

    ## 3. fit s-reps to those ellipsoids
    top_srep = sm.fit_srep(ell_meshes[0], ell_sreps[0], final_mesh_names[0])
    bot_srep = sm.fit_srep(ell_meshes[1], ell_sreps[1], final_mesh_names[1])
    top_srep_file_name = os.path.join(output_heatmap_folder, 'top_srep_' + moment_name + '.vtk')
    srep_writer = vtk.vtkPolyDataWriter()
    srep_writer.SetInputData(top_srep)
    srep_writer.SetFileName(top_srep_file_name)
    srep_writer.Update()

    bot_srep_file_name = os.path.join(output_heatmap_folder, 'bot_srep_' + moment_name + '.vtk')
    srep_writer.SetInputData(bot_srep)
    srep_writer.SetFileName(bot_srep_file_name)
    srep_writer.Update()
def compute_heatmap(output_folder, moment_name):
    std_mesh_path = os.path.join(output_folder, 'final_mesh', 'std_ell_label_SPHARM.vtk')
    ell_meshes, ell_sreps = construct_two_standard_ellipsoids(std_mesh_path)
    variation_modes_folder = 'population_variation'
    output_heatmap_folder = os.path.join(output_folder, variation_modes_folder)
    bot_reader = vtk.vtkPolyDataReader()
    bot_reader.SetFileName(os.path.join(output_heatmap_folder, 'bot_' + moment_name + '_label_SPHARM.vtk'))
    bot_reader.Update()
    bot_vtk = bot_reader.GetOutput()
    dist_filter = vtk.vtkDistancePolyDataFilter()
    dist_filter.SetSignedDistance(False)
    dist_filter.SetInputData(0, bot_vtk)
    dist_filter.SetInputData(1, ell_meshes[1])
    dist_filter.Update()

    top_reader = vtk.vtkPolyDataReader()
    top_reader.SetFileName(os.path.join(output_heatmap_folder, 'top_' + moment_name + '_label_SPHARM.vtk'))
    top_reader.Update()
    top_vtk = top_reader.GetOutput()
    top_dist_filter = vtk.vtkDistancePolyDataFilter()
    top_dist_filter.SetSignedDistance(False)
    top_dist_filter.SetInputData(0, top_vtk)
    top_dist_filter.SetInputData(1, ell_meshes[0])
    top_dist_filter.Update()

    # append_filter = vtk.vtkAppendPolyData()
    # append_filter.AddInputData(dist_filter.GetOutput())
    # append_filter.AddInputData(top_dist_filter.GetOutput())
    # append_filter.Update()
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(os.path.join(output_heatmap_folder, 'top_'+ moment_name + '_heatmap.vtk'))
    writer.SetInputData(top_dist_filter.GetOutput())
    writer.Update()

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(os.path.join(output_heatmap_folder, 'bot_' + moment_name + '_heatmap.vtk'))
    writer.SetInputData(dist_filter.GetOutput())
    writer.Update()

if __name__ == '__main__':
    population_folder = '/playpen/workspace/Simulate_Shapes/simulation_for_links/corr_population_a/'
    m_name = 'sigma'
    bend_top = 0.1/3.0 + 0.11
    compute_moment_shape(population_folder, [bend_top, 0.9 * bend_top], [0, 0.04], moment_name=m_name)
    compute_heatmap(population_folder, moment_name=m_name)

    m_name = 'minus_sigma'
    bend_top = -0.1 / 3.0 + 0.11
    compute_moment_shape(population_folder, [bend_top, 0.9 * bend_top], [0, 0.04], moment_name=m_name)
    compute_heatmap(population_folder, moment_name=m_name)

    m_name = 'mean'
    compute_moment_shape(population_folder, [0.11, 0.9 * 0.11], [0.0, 0.04], moment_name=m_name)
    compute_heatmap(population_folder, moment_name=m_name)
    population_folder = '/playpen/workspace/Simulate_Shapes/simulation_for_links/corr_population_b/'
    m_name = 'sigma'
    bend = [0.1, 0.09]
    twist = [0.1 / 3.0, 0.1 / 3.0]
    compute_moment_shape(population_folder, bend, twist, moment_name=m_name)
    compute_heatmap(population_folder, moment_name=m_name)

    m_name = 'minus_sigma'
    twist = [-0.1 / 3.0, -0.1 / 3.0]
    compute_moment_shape(population_folder, bend, twist, moment_name=m_name)
    compute_heatmap(population_folder, moment_name=m_name)

    m_name = 'mean'
    twist = [0.0, 0.0]
    compute_moment_shape(population_folder, bend, twist, moment_name=m_name)
    compute_heatmap(population_folder, moment_name=m_name)
    # population_folder = '/playpen/workspace/Simulate_Shapes/simulation_for_links/corr_population_a/'
    # std_mesh_path = os.path.join(population_folder, 'final_mesh', 'std_ell_label_SPHARM.vtk')
    # construct_two_standard_ellipsoids(std_mesh_path, population_folder)
    print('Done')