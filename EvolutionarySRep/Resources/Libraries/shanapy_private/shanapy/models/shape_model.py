import itk
import vtk
import numpy as np
import os
eps = np.finfo(float).eps
def spharm_mesh(image_path, label=255, output_file_folder=None, flip_template=''):
    """
    This function generate a moderately dense triangle
    mesh for the object (of label) in the image (in a nrrd file).

    Such a triangle mesh is smooth and has good correspondence with
    other samples in the population
    """
    file_name = str.split(os.path.split(image_path)[1], '.')[0]
    output_file_name = os.path.join(output_file_folder, file_name)
    cmd_gen_para_mesh = '../../third_party/spharm_bin/GenParaMeshCLP --label ' + str(label) + \
        ' ' + image_path + ' ' + output_file_name + '_para.vtk ' + output_file_name + '_surf.vtk'
    os.system(cmd_gen_para_mesh)

    cmd_gen_spharm_mesh = '../../third_party/spharm_bin/ParaToSPHARMMeshCLP ' + output_file_name + '_para.vtk ' \
                          + output_file_name + '_surf.vtk ' +  output_file_name + '_ --subdivLevel 10 --flipTemplateOn --flipTemplate ' + flip_template
#    print(cmd_gen_spharm_mesh)
    os.system(cmd_gen_spharm_mesh)
def surface_mesh(image, output_file_name = None):
    """
    Generate surface mesh for the input image
    if output_file_name is not None, save two files
    for this surface mesh, one is vtkPolyData,
    the other is stl file for better use of meshio
    """

    ## threshold image into binary (assume object value in the range [20, 256])
    thresholded = itk.BinaryThresholdImageFilter[type(image), itk.Image[itk.SS, 3]].New()
    thresholded.SetInput(image)
    thresholded.SetLowerThreshold(20)
    thresholded.SetUpperThreshold(256)
    thresholded.SetOutsideValue(0)
    thresholded.SetInsideValue(255)
    thresholded.Update()

    ## convert to mesh
    mesh = itk.BinaryMask3DMeshSource[itk.Image[itk.SS, 3], itk.Mesh[itk.F, 3]].New()
    mesh.SetInput(thresholded.GetOutput())
    mesh.SetObjectValue(255)
    mesh.Update()
    mesh_out = mesh.GetOutput()

    if output_file_name is not None:
        ## write to file (*.stl for connection with pygalmesh)
        mesh_writer = itk.MeshFileWriter[type(mesh_out)].New()
        mesh_writer.SetInput(mesh_out)
        mesh_writer.SetFileName(output_file_name)
        mesh_writer.Update()

        vtk_reader = vtk.vtkPolyDataReader()
        vtk_reader.SetFileName(output_file_name)
        vtk_reader.Update()
        stl_writer = vtk.vtkSTLWriter()
        stl_file_name = os.path.splitext(output_file_name)[0] + '.stl'
        stl_writer.SetFileName(stl_file_name)
        stl_writer.SetInputData(vtk_reader.GetOutput())
        stl_writer.Write()

    return mesh_out
def postprocess_mesh(input_file_name,
                        output_file_name=None):
    """
    Smooth and coarsen mesh designated by input_file_name (*.stl)
    with default parameters in each function
    """
    coarsen_mesh(input_file_name, output_file_name=output_file_name)
    reader = vtk.vtkSTLReader()
    reader.SetFileName(output_file_name)
    reader.Update()
    return smooth_mesh(reader.GetOutput())
def coarsen_mesh(input_file_name='standard_ellipsoid.stl',\
                 edge_size = 0.4,  \
                 facet_angle=25,   \
                 facet_size=2,     \
                 facet_distance=2, \
                 output_file_name=None):
    """
    Coarsen the mesh
    """
    import pygalmesh

    mesh = pygalmesh.remesh_surface(
        input_file_name,
        edge_size=edge_size,
        facet_angle=facet_angle,
        facet_size=facet_size,
        facet_distance=facet_distance,
        verbose=False,
    )
    if output_file_name is None:
        output_file_name = 'coarsen_ellipsoid.stl'
    mesh.write(output_file_name)
    return mesh

def smooth_mesh(mesh, pass_band = 0.01, num_iter = 40):
    """
    Use vtk filters (vtkWindowedSincPolyDataFilter) to  smooth the surface
    Input mesh: vtkPolyData
    """
    smooth_filter = vtk.vtkWindowedSincPolyDataFilter()
    smooth_filter.SetInputData(mesh)
    smooth_filter.SetPassBand(pass_band)
    smooth_filter.SetNumberOfIterations(num_iter)
    smooth_filter.Update()
    out_mesh = smooth_filter.GetOutput()
    return out_mesh

def fit_srep(ell_mesh, ell_srep, target_mesh_file_name, output_folder='../data'):
    """
    Given the original ellipsoid and the deformed
    apply TPS to estimate the stratified deformation.
    Apply the deformation on spokes
    """
    # reader = vtk.vtkPolyDataReader()
    # reader.SetFileName(ell_mesh_file_name)
    # reader.Update()
    # ell_mesh = reader.GetOutput()
    # transformer = vtk.vtkTransform()
    # transformer.Translate(center)
    # transformer.Update()
    # transform_filter = vtk.vtkTransformPolyDataFilter()
    # transform_filter.SetInputData(ell_mesh)
    # transform_filter.SetTransform(transformer)
    # transform_filter.Update()
    # ell_mesh = transform_filter.GetOutput()

    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(target_mesh_file_name)
    reader.Update()
    target_mesh = reader.GetOutput()
    assert ell_mesh.GetNumberOfPoints() == target_mesh.GetNumberOfPoints(), \
                "Different number of surface points between ellipsoid and target"
    num_pts = ell_mesh.GetNumberOfPoints()

    source_pts = vtk.vtkPoints()
    target_pts = vtk.vtkPoints()
    for i in range(num_pts):
        pt = [0] * 3
        ell_mesh.GetPoint(i, pt)
        source_pts.InsertNextPoint(pt)

        target_mesh.GetPoint(i, pt)
        target_pts.InsertNextPoint(pt)
    source_pts.Modified()
    target_pts.Modified()

    ### Interpolate deformation with thin-plate-spline
    tps = vtk.vtkThinPlateSplineTransform()
    tps.SetSourceLandmarks(source_pts)
    tps.SetTargetLandmarks(target_pts)
    tps.SetBasisToR()
    tps.Modified()

    ### Apply the deformation onto the spokes
    deformed_srep = vtk.vtkPolyData()
    deformed_spokes_ends = vtk.vtkPoints()
    deformed_spoke_lines = vtk.vtkCellArray()
    # refined_srep is a polydata that collects spokes
    for i in range(ell_srep.GetNumberOfCells()):
        base_pt_id = i * 2
        bdry_pt_id = i * 2 + 1
        s_pt = ell_srep.GetPoint(base_pt_id)
        b_pt = ell_srep.GetPoint(bdry_pt_id)

        new_s_pt = tps.TransformPoint(s_pt)
        new_b_pt = tps.TransformPoint(b_pt)

        id0 = deformed_spokes_ends.InsertNextPoint(new_s_pt)
        id1 = deformed_spokes_ends.InsertNextPoint(new_b_pt)

        spoke_line = vtk.vtkLine()
        spoke_line.GetPointIds().SetId(0, id0)
        spoke_line.GetPointIds().SetId(1, id1)
        deformed_spoke_lines.InsertNextCell(spoke_line)
    deformed_srep.SetPoints(deformed_spokes_ends)
    deformed_srep.SetLines(deformed_spoke_lines)
    deformed_srep.Modified()
    return deformed_srep
