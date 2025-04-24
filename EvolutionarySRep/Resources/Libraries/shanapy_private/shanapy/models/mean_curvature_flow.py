"""
contains a bunch of steps from Zhiyuan's original try at the MCF method.
"""
from vtk.util import numpy_support
from shanapy.utils import viz
from shanapy.models import srep_fitter
import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
from glob import glob
import pyvista
import vtk
import scipy
import os
import re
from shanapy.utils import obtain_class_labels as ocl
iter_num = 500
def rotate_polydata_pts(rot_mat, trans_ellip):
    """
    Input rot_mat 3x3 rotation matrix
    Input trans_ellip vtkPolyData
    Return vtkPoints
    """
    np_ellipsoid_vertices = []
    for i in range(trans_ellip.GetNumberOfPoints()):
        np_ellipsoid_vertices.append(np.array(trans_ellip.GetPoint(i)))
    np_ellipsoid_vertices = np.array(np_ellipsoid_vertices)

    rot_ellipsoid = np.matmul(np_ellipsoid_vertices, rot_mat)
    ellip_pts = trans_ellip.GetPoints()
    for i in range(trans_ellip.GetNumberOfPoints()):
        ellip_pts.SetPoint(i, rot_ellipsoid[i, :])
    ellip_pts.Modified()
    return ellip_pts
def get_thin_plate_spline_deform(input_target_mesh, input_source_mesh):
    ## sparsen boundary points to speed up computation
    # clean_ratio = 0.1
    # clean_data_polydata = vtk.vtkCleanPolyData()
    # clean_data_polydata.SetInputData(input_target_mesh)
    # clean_data_polydata.SetTolerance(clean_ratio)
    # clean_data_polydata.Update()
    # target_mesh = clean_data_polydata.GetOutput()
    
    # source_clean_data_polydata = vtk.vtkCleanPolyData()
    # source_clean_data_polydata.SetTolerance(clean_ratio)
    # source_clean_data_polydata.SetInputData(input_source_mesh)
    # source_clean_data_polydata.Update()
    # source_mesh = source_clean_data_polydata.GetOutput()

    target_mesh = input_target_mesh
    source_mesh = input_source_mesh
    # compute_distance_between_poly(target_mesh, source_mesh)
    source_pts = vtk.vtkPoints()
    target_pts = vtk.vtkPoints()
    for i in range(0, target_mesh.GetNumberOfPoints(), 10):
        source_pts.InsertNextPoint(source_mesh.GetPoint(i))
        target_pts.InsertNextPoint(target_mesh.GetPoint(i))
    tps = vtk.vtkThinPlateSplineTransform()
    tps.SetSourceLandmarks(source_pts)
    tps.SetTargetLandmarks(target_pts)
    tps.SetBasisToR()
    tps.Modified()

    return tps
def compute_distance_between_poly(poly1, poly2):
    distance_filter = vtk.vtkDistancePolyDataFilter()
    distance_filter.SetInputData(0, poly1)
    distance_filter.SetInputData(1, poly2)
    distance_filter.Update()
    dist = np.array(distance_filter.GetOutput().GetPointData().GetScalars().GetRange())
    print(dist)
def flow(mesh_file, show_flow=False):
    """
    Computes mean curvature flow steps on mesh
    Uses taubin_smooth technique for smoothing at each step
    returns near-ellipsoid mesh and all intermediate meshes
    """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(mesh_file)
    reader.Update()

    mesh = reader.GetOutput()
    tol = 0.05
    q = 1.0
    dt = 0.001
    orig_mesh = mesh
    prev_mesh = mesh
    thin_plate_spline_list = []
    for i in range(iter_num + 1):
        deformed_surface_writer = vtk.vtkPolyDataWriter()
        deformed_surface_writer.SetFileName('../../data/forward/' + str(i) + '.vtk')
        deformed_surface_writer.SetInputData(mesh)
        deformed_surface_writer.Update()

        taubin_smooth = vtk.vtkWindowedSincPolyDataFilter()
        taubin_smooth.SetInputData(mesh)
        taubin_smooth.SetNumberOfIterations(20)
        taubin_smooth.BoundarySmoothingOff()
        taubin_smooth.FeatureEdgeSmoothingOff()
        taubin_smooth.SetPassBand(0.01)
        taubin_smooth.NonManifoldSmoothingOn()
        taubin_smooth.NormalizeCoordinatesOn()

        taubin_smooth.Update()
        mesh = taubin_smooth.GetOutput()

        normal_generator = vtk.vtkPolyDataNormals()
        normal_generator.SetInputData(mesh)
        normal_generator.SplittingOff()
        normal_generator.ComputePointNormalsOn()
        normal_generator.ComputeCellNormalsOff()
        normal_generator.Update()
        mesh = normal_generator.GetOutput()

        curvatures = vtk.vtkCurvatures()
        curvatures.SetCurvatureTypeToMean()
        curvatures.SetInputData(mesh)
        curvatures.Update()

        mean_curvatures = curvatures.GetOutput().GetPointData().GetArray("Mean_Curvature")
        normals = normal_generator.GetOutput().GetPointData().GetNormals()

        mesh_pts = mesh.GetPoints()
        deform_fields = []
        pv_mean_curvatures = []
        for j in range(mesh.GetNumberOfPoints()):
            current_point = mesh.GetPoint(j)
            current_normal = np.array(normals.GetTuple3(j))
            current_mean_curvature = mean_curvatures.GetValue(j)

            pt = np.array(mesh_pts.GetPoint(j))
            pt -= dt * current_mean_curvature * current_normal
            deform_fields.append(pyvista.Arrow(pt, -dt * current_mean_curvature * current_normal))
            pv_mean_curvatures.append(-current_mean_curvature)
            mesh_pts.SetPoint(j, pt)
        if show_flow:
            plotter = pyvista.Plotter()
            pv_mesh = pyvista.PolyData(mesh)
            pv_mesh.point_arrays['mean_curvatures'] = pv_mean_curvatures
            pv_mesh.set_active_scalars('mean_curvatures')
            plotter.add_mesh(pv_mesh, opacity=0.4)
            for arrow in deform_fields:
                plotter.add_mesh(arrow)
            plotter.show()
        mesh_pts.Modified()
        mesh.SetPoints(mesh_pts)
        mesh.Modified()

        tps_deform = get_thin_plate_spline_deform(prev_mesh, mesh)

        prev_mesh = mesh
        thin_plate_spline_list.append(tps_deform)
    # new_mesh = apply_tps_on_spokes_poly(mesh, thin_plate_spline_list)
    # plotter = pyvista.Plotter()
    # plotter.add_mesh(new_mesh, color='red', opacity=0.3)
    # plotter.add_mesh(orig_mesh, color='blue', opacity=0.3)
    # plotter.show()

    return mesh, thin_plate_spline_list
def fit_srep_to_quasi_ellipsoid(mesh, ax = None, auto_flip=True):
    """
    Fit an srep the quasi-ellipsoid resulting from MCF
    
    """
    np_vertices = []
    for i in range(mesh.GetNumberOfPoints()):
        np_vertices.append(np.array(mesh.GetPoint(i)))
    np_vertices = np.array(np_vertices)
    center = np.mean(np_vertices, axis=0)

    centered_vertices = np_vertices - center[None, :]
    w, v = np.linalg.eig(np.matmul(centered_vertices.T, centered_vertices))
    idx = w.argsort()[::-1]
    w = w[idx]
    v = v[:, idx]

    ## correct the ellipsoid based on volume
    r0, r1, r2 = np.sqrt(w[0]), np.sqrt(w[1]), np.sqrt(w[2])
    ellipsoid_volume = 4 / 3.0 * np.pi * r0 * r1 * r2
    mass_filter = vtk.vtkMassProperties()
    mass_filter.SetInputData(mesh)
    mass_filter.Update()

    volume_factor = pow(mass_filter.GetVolume() / ellipsoid_volume, 1.0 / 3.0)
    r0 *= volume_factor
    r1 *= volume_factor
    r2 *= volume_factor

   # print(r0, r1, r2)
    ellipsoid_param = vtk.vtkParametricEllipsoid()
    ellipsoid_param.SetXRadius(r0)
    ellipsoid_param.SetYRadius(r1)
    ellipsoid_param.SetZRadius(r2)
    param_funct = vtk.vtkParametricFunctionSource()
    param_funct.SetUResolution(30)
    param_funct.SetVResolution(30)
    param_funct.SetParametricFunction(ellipsoid_param)
    param_funct.Update()
    best_fitting_ellipsoid = param_funct.GetOutput()
    trans_ellip = best_fitting_ellipsoid
    ##### To keep the orientation of ellipsoids consistent in a population
    # first_ellipsoid_file = '../../data/best_fitting_ellipsoid.npy'
    # if not os.path.exists(first_ellipsoid_file):
    #     with open(first_ellipsoid_file, 'wb') as f:
    #         print("save reference principal directions")
    #         np.save(f, v)
    # elif auto_flip:
    #     with open(first_ellipsoid_file, 'rb') as f:
    #         v_ref = np.load(f)

    #         cos_ev1 = np.dot(v_ref[:, 0], v[:, 0])
    #         cos_ev2 = np.dot(v_ref[:, 1], v[:, 1])
    #         cos_ev3 = np.dot(v_ref[:, 2], v[:, 2])

    #         # sin_ev1 = np.linalg.norm(np.cross(v_ref[:, 0], v[:, 0]))
    #         # sin_ev2 = np.linalg.norm(np.cross(v_ref[:, 1], v[:, 1]))
    #         # sin_ev3 = np.linalg.norm(np.cross(v_ref[:, 2], v[:, 2]))
    #         # print(cos_ev1, sin_ev1)
    #         if cos_ev1 < 0:
    #             print('revert 0')
    #             v[:, 0] *= -1
    #         if cos_ev2 < 0:
    #             print('revert 1')
    #             v[:, 1] *= -1
    #         if cos_ev3 < 0:
    #             print('revert 2')
    #             v[:, 2] *= -1
    # ## rotation
    # rot_mat = v.T
    # # 1. rotate surface mesh
    # ellip_pts = rotate_polydata_pts(rot_mat, trans_ellip)
    # trans_ellip.SetPoints(ellip_pts)
    # trans_ellip.Modified()

    # ## translation
    # transform = vtk.vtkTransform()
    # transform.Translate(center[0], center[1], center[2])
    # trans_filter = vtk.vtkTransformPolyDataFilter()
    # trans_filter.SetInputData(trans_ellip)
    # trans_filter.SetTransform(transform)
    # trans_filter.Update()
    # trans_ellip = trans_filter.GetOutput()
    ##########################################################################
    # ## fit s-reps to the ellipsoid
    srep_poly, up_spokes_poly, down_spokes_poly, crest_spokes_poly = srep_fitter.fit_ellipsoid(trans_ellip)#, predefined_eig_vec=v)
    deformed_up_srep_writer = vtk.vtkPolyDataWriter()
    deformed_up_srep_writer.SetFileName('../../data/model/up' + str(iter_num) + '.vtk')
    deformed_up_srep_writer.SetInputData(up_spokes_poly)
    deformed_up_srep_writer.Update()

    deformed_down_srep_writer = vtk.vtkPolyDataWriter()
    deformed_down_srep_writer.SetFileName('../../data/model/down' + str(iter_num) + '.vtk')
    deformed_down_srep_writer.SetInputData(down_spokes_poly)
    deformed_down_srep_writer.Update()

    deformed_crest_srep_writer = vtk.vtkPolyDataWriter()
    deformed_crest_srep_writer.SetFileName('../../data/model/crest' + str(iter_num) + '.vtk')
    deformed_crest_srep_writer.SetInputData(crest_spokes_poly)
    deformed_crest_srep_writer.Update()
    # plotter = pyvista.Plotter()
    # plotter.add_mesh(trans_ellip, color='red', opacity=0.3)
    # plotter.add_mesh(mesh, color='blue', opacity=0.3)
    # plotter.add_mesh(srep_poly)
    # plotter.show()

    ##### Debug code for visualization of consistent orientation
    # plotter = pyvista.Plotter()
    # plotter.add_mesh(mesh, opacity=0.2, color='red')
    # plotter.add_mesh(trans_ellip)

    # top_pt = np.array(trans_ellip.GetPoint(0))
    # plotter.add_arrows(top_pt, v[:, 0])
    # plotter.add_arrows(top_pt, v[:, 1])
    # plotter.add_arrows(top_pt, v[:, 2])

    # plotter.show()

    ##visualize s-rep and mesh
    # plotter = pyvista.Plotter()
    # plotter.add_mesh(trans_ellip, opacity=0.2, color='blue')
    # plotter.add_mesh(srep_poly, color='red')
    # plotter.add_mesh(mesh, opacity=0.3)
    # plotter.show()
    return v, srep_poly
def apply_tps_on_spokes_poly(srep_poly, tps_list):
    """
    vtk thin plate spline method is too slow
    """
    # print(len(tps_list))

    deformed_srep = vtk.vtkPolyData()
    deformed_srep.DeepCopy(srep_poly)

    deformations = []

    new_pts = vtk.vtkPoints()
    for i, tps in enumerate(tps_list[::-1]):
        for j in range(deformed_srep.GetNumberOfPoints()):
            new_pt = tps.TransformPoint(np.array(deformed_srep.GetPoint(j)))
            new_pts.InsertNextPoint(new_pt)
        # bdry_poly, skeletal_poly = viz.form_strata_from_spokes(deformed_srep)
        # transform_filter = vtk.vtkTransformPolyDataFilter()
        # transform_filter.SetTransform(tps)
        # transform_filter.SetInputData(bdry_poly)
        # transform_filter.Update()
        # deformed_bdry = transform_filter.GetOutput()

        # skeletal_transform_filter = vtk.vtkTransformPolyDataFilter()
        # skeletal_transform_filter.SetTransform(tps)
        # skeletal_transform_filter.SetInputData(skeletal_poly)
        # skeletal_transform_filter.Update()
        # deformed_skeletal = skeletal_transform_filter.GetOutput()

        # deformed_srep = viz.form_spokes_from_strata(deformed_skeletal, deformed_bdry)

        # transform_filter = vtk.vtkTransformPolyDataFilter()
        # transform_filter.SetTransform(tps)
        # transform_filter.SetInputData(deformed_srep)
        # transform_filter.Update()
        # deformed_srep = transform_filter.GetOutput()
        # plotter = pyvista.Plotter()
        # plotter.add_mesh(deformed_srep)
        # plotter.show()

        deformed_srep.SetPoints(new_pts)
        deformed_srep.Modified()
    return deformed_srep
def copy_mesh_file(input_dir='/playpen/data/non-aligned/CaudL/', target_dir='/playpen/workspace/my_paper/linking/data/nonaligned_caud_sreps/'):
    pos_ids, neg_ids = ocl.load_class_labels()
    for patient_folder in os.listdir(target_dir):
        if patient_folder not in pos_ids and patient_folder not in neg_ids:
            cmd = 'rm -r ' + target_dir + patient_folder
            os.system(cmd)
            continue
        mesh_file_name = None
        for files in glob(target_dir + patient_folder + "/*.vtk"):
            mesh_file_name = files
        if mesh_file_name is None:
            cmd = 'cp ' + input_dir + 'stx_noscale_' + patient_folder + '*.vtk ' + target_dir + patient_folder
            os.system(cmd)
            
def batch_initialize_sreps(input_dir='/playpen/data/non-aligned/CaudL', output_dir='/playpen/workspace/my_paper/linking/data/nonaligned_caud_sreps/', group_ids=None):
    # dummy_pt = np.array([[0, 0, 1]])
    # ax = viz.viz_directions_distribution(dummy_pt, dummy_pt, show_now=False)
    ev1, ev2, ev3 = [], [], []
    process_num = 0
    spoke0 = []
    for mesh_file in os.listdir(input_dir):
        # if process_num < 95:
        #     process_num += 1
        #     continue

        if re.match(r"(.*)\.vtk", mesh_file) == None:
            continue

        patient_folder = mesh_file[12:18]

        is_pos = False
        if patient_folder not in group_ids:
            print('Not considered')
            continue

        # if os.path.exists(os.path.join(output_dir, patient_folder, 'crest.vtp')):
        #     print('Finished')
        #     continue
    
#         elif patient_folder in pos_ids:
#             print('Pos case')
#             #continue
#             is_pos = True
#         else:
# #            continue
#             print('Neg case')
        process_num += 1
        mesh_file_path = os.path.join(input_dir, mesh_file)
        ell_mesh, tps_list = flow(mesh_file_path)
        eigen_vectors, srep_poly = fit_srep_to_quasi_ellipsoid(ell_mesh)
        # pt0 = np.array(srep_poly.GetPoint(0))
        # pt1 = np.array(srep_poly.GetPoint(1))
        # dir0 = (pt1 - pt0) / np.linalg.norm(pt1-pt0)
        # spoke0.append(dir0)
        # obj_srep = apply_tps_on_spokes_poly(srep_poly, tps_list)

        # reader = vtk.vtkPolyDataReader()
        # reader.SetFileName(mesh_file_path)
        # reader.Update()

        # mesh = reader.GetOutput()
        # plotter = pyvista.Plotter()
        # plotter.add_mesh(mesh, opacity=0.4)
        # plotter.add_mesh(obj_srep)
        # plotter.show()

        ev1.append(eigen_vectors[:, 0])
        ev2.append(eigen_vectors[:, 1])
        ev3.append(eigen_vectors[:, 2])

        output_file_path = os.path.join(output_dir, patient_folder)
        if not os.path.exists(output_file_path):
            print('mkdir ' + output_file_path)
            os.system('mkdir ' + output_file_path)

        cmd = 'Slicer --no-main-window --python-script backward_mcf.py ' + output_file_path
        print(cmd)
        os.system(cmd)
        # if process_num > 10:
        #     break
        print('Finished ' + str(process_num) + ' cases')
    ev1 = np.array(ev1)
    ev2 = np.array(ev2)
    ev3 = np.array(ev3)
    # with open('../../data/eigen_vectors1.npy', 'wb') as ev1_f:
    #     np.save(ev1_f, ev1)
    # with open('../../data/eigen_vectors2.npy', 'wb') as ev2_f:
    #     np.save(ev2_f, ev2)
    # with open('../../data/eigen_vectors3.npy', 'wb') as ev3_f:
    #     np.save(ev3_f, ev3)

    # spoke0s = np.array(spoke0)
    # print(spoke0s.shape)
    # with open('../../data/spoke0.npy', 'wb') as spoke0_f:
    #     np.save(spoke0_f, spoke0s)

    # viz.viz_directions_distribution(ev1, spoke0s)
if __name__ == '__main__':
    mesh_file_path = '/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/103430/stx_noscale_103430_V06_t1w_RAI_Bias_label_pp_surfSPHARM.vtk'
    ell_mesh, tps_list = flow(mesh_file_path, True)
    # pos_ids, neg_ids = ocl.load_class_labels()
    # #batch_initialize_sreps(input_dir='/playpen/data/non-aligned/CaudL', output_dir='/playpen/workspace/my_paper/linking/data/nonaligned_caud_sreps/', group_ids=pos_ids)
    # #batch_initialize_sreps(input_dir='/playpen/data/non-aligned/CaudL', output_dir='/playpen/workspace/my_paper/linking/data/nonaligned_caud_sreps/', group_ids=neg_ids)

    # #batch_initialize_sreps(input_dir='/playpen/data/non-aligned/HippL', output_dir='/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/', group_ids=pos_ids)
    # batch_initialize_sreps(input_dir='/playpen/data/non-aligned/HippL', output_dir='/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/', group_ids=neg_ids)

    #copy_mesh_file()
    print('done')

# import numpy as np
# from shanapy.utils import viz
# with open('../../data/eigen_vectors1.npy', 'rb') as ev1_f:
#     ev1 = np.load(ev1_f)
# with open('../../data/eigen_vectors2.npy', 'rb') as ev2_f:
#     ev2 = np.load(ev2_f)
# with open('../../data/eigen_vectors3.npy', 'rb') as ev3_f:
#     ev3 = np.load(ev3_f)

# viz.viz_directions_distribution(ev3, ev1)
