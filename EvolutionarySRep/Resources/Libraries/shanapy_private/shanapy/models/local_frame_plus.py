import vtk
import os
import pyvista as pv
from shanapy.utils import visualize_srep_fit as vsf, viz
from shanapy.models import srep_fitter as sf, principal_nested_spheres as pns

import numpy as np
delta = 0.1#0.001
def form_basis_vectors(curr_pt, next_pt, pt_next_radial_line, flip_tau = False, normal_end=None, flip_theta=True):
    """
    next_pt and curr_pt form direction along tau
    pt_next_radial_line is on the direction of which theta increases
    """
    curr_pt = curr_pt.astype('float').squeeze()
    next_pt = next_pt.astype('float').squeeze()
    pt_next_radial_line = pt_next_radial_line.astype('float').squeeze()
    tau_vec = next_pt - curr_pt
    tau_vec /= np.linalg.norm(tau_vec)
    if flip_tau: tau_vec *= -1

    theta_vec = pt_next_radial_line - curr_pt
    theta_vec /= np.linalg.norm(theta_vec)
    if flip_theta: theta_vec *= -1

    if normal_end is None:
        n_vec = np.cross(tau_vec, theta_vec)
        n_vec /= np.linalg.norm(n_vec)
        new_theta_vec = np.cross(n_vec, tau_vec)
        new_theta_vec /= np.linalg.norm(new_theta_vec)
    else:
        new_theta_vec = theta_vec
        n_vec = normal_end - curr_pt
        n_vec /= np.linalg.norm(n_vec)
    return tau_vec, new_theta_vec, n_vec
def selected_frames_feats_within_object(feats, num_feats_per_spoke = 13):
    ## TODO: double check if the selected are the required
    #selected_spoke_ids = [0, 72//2, 72, 72+72//2, 72//4-3, 72 + 72//4 - 3, 143 + 24//4, 169 - 24//4]
    selected_spoke_ids = np.arange(168)

    n, d = feats.shape
    local_frame_features = np.array([]).reshape(n, 0)
    for spoke_id in selected_spoke_ids:
        start_idx = spoke_id * num_feats_per_spoke
        # skeletal_pt = feats[:, start_idx:start_idx+3]
        # bdry_pt = feats[:, start_idx+6:start_idx+9]
        # spoke_dir_skeletal_f = feats[:, start_idx+3:start_idx+6]
        # spoke_dir_bdry_f = feats[:, start_idx+9:start_idx+num_feats_per_spoke-1]
        # spoke_length = feats[:, start_idx+num_feats_per_spoke-1]
        local_frame_features = np.concatenate((local_frame_features, feats[:, start_idx:start_idx+num_feats_per_spoke]), axis=1)

    return local_frame_features, selected_spoke_ids
def frechet_mean(vecs):
    pns_model = pns.PNS(vecs)
    pns_model.fit()
    _, PNS_coords = pns_model.output
    mean_vec = PNS_coords['mean']
    return mean_vec
def selected_frames_feats_global(hipp_frames, caud_frames, selected_spoke_ids):
    """
    Frames features w.r.t. central coordinate system
    Return features include origin positions
    """
    up_center_spoke_id = 72//4-3#selected_spoke_ids[4]
    groups, case_num, spoke_num, frame_pts_num, _ = hipp_frames.shape
    all_origins_ccs = []
    all_normal_dirs = []
    all_ccs = []
    object_centers = []
    for group_frames in range(groups):
        for case_id in range(case_num):
            ## compute origin of CCS for each case
            curr_hipp_frame = hipp_frames[group_frames][case_id]
            curr_caud_frame = caud_frames[group_frames][case_id]

            hipp_center = curr_hipp_frame[up_center_spoke_id][0]
            caud_center = curr_caud_frame[up_center_spoke_id][0]
            object_centers.append(hipp_center)
            object_centers.append(caud_center)

            ccs_center = (hipp_center + caud_center) / 2
            ## compute frame vector directions by average corresponding directions with PNS
            tau_vecs, theta_vecs, normal_vecs = [], [], []
            hipp_normals = []
            caud_normals = []
            for spoke_id in selected_spoke_ids:
                ## add hippocampus frames
                center_pt = curr_hipp_frame[spoke_id][0]
                tau_end = curr_hipp_frame[spoke_id][1]
                theta_end = curr_hipp_frame[spoke_id][2]
                normal_end = curr_hipp_frame[spoke_id][3]
                hipp_normals.append(np.concatenate((tau_end, theta_end, normal_end)))
                ## add caudate frames
                center_pt = curr_caud_frame[spoke_id][0]
                tau_end = curr_caud_frame[spoke_id][1]
                theta_end = curr_caud_frame[spoke_id][2]
                normal_end = curr_caud_frame[spoke_id][3]
                caud_normals.append(np.concatenate((tau_end, theta_end, normal_end)))

            all_normal_dirs.append(np.array(hipp_normals).flatten())
            all_normal_dirs.append(np.array(caud_normals).flatten())

            all_ccs.append(ccs_center)

            ## compute relative coordinates of each origin for current hipp and caud
            curr_hipp_origins = []
            curr_caud_origins = []
            for spoke_id in selected_spoke_ids:
                curr_hipp_center_pt = curr_hipp_frame[spoke_id][0]
#                hipp_vect_to_global_center = coord_in_ccs(ccs, curr_hipp_center_pt)
                vect_frame_to_hipp_center = curr_hipp_center_pt# - hipp_center
                curr_hipp_origins.append(vect_frame_to_hipp_center)# np.concatenate((hipp_vect_to_global_center, vect_frame_to_hipp_center)))
                curr_caud_center_pt = curr_caud_frame[spoke_id][0]
#                caud_vect_to_global_center = coord_in_ccs(ccs, curr_caud_center_pt)
                vect_frame_to_caud_center = curr_caud_center_pt# - caud_center
                curr_caud_origins.append(vect_frame_to_caud_center)# np.concatenate((caud_vect_to_global_center, vect_frame_to_caud_center)))
            hipp_feats_css = np.array(curr_hipp_origins).flatten()
            caud_feats_css = np.array(curr_caud_origins).flatten()
            all_origins_ccs.append(hipp_feats_css)
            all_origins_ccs.append(caud_feats_css)
    return np.array(all_normal_dirs, dtype=float), np.array(all_origins_ccs, dtype=float), all_ccs, object_centers

def coord_in_ccs(ccs, pt):
    origin_ccs = ccs[0]
    return pt - origin_ccs # same with global c.s.
    x = np.dot(pt - origin_ccs, ccs[1])
    y = np.dot(pt - origin_ccs, ccs[2])
    z = np.dot(pt - origin_ccs, ccs[3])
    return np.array([x, y, z], dtype=float)
def local_frame_near_skeleton(obj_folder, input_id, alpha=1, input_skeleton_tps=None, initialize=False):
    target_srep_xml = '/playpen/workspace/my_paper/linking/data/' + obj_folder + '/' + input_id + '/header.xml'
    target_surf_file = '/playpen/workspace/my_paper/linking/data/' + obj_folder + '/' + input_id + '/stx_noscale_' + input_id + '_V06_t1w_RAI_Bias_label_pp_surfSPHARM.vtk'
    if not os.path.exists(target_surf_file):
        print("Surface not exists:", target_surf_file)
        return None, None, None, None, None, None
    ## initialize s-rep
    if initialize:
        output_path = '/playpen/workspace/my_paper/linking/data/' + obj_folder + '/' + input_id
        if not os.path.exists(output_path + '/header.xml'):
            print("Creating s-rep for ", input_id)
            os.system("Slicer --no-main-window --python-script /playpen/workspace/Simulate_Shapes/shanapy/models/cmd_initialize_srep.py " + target_surf_file + " " + output_path)

    ## special frames on ellipsoid
    ell_mesh_reader = vtk.vtkPolyDataReader()
    ell_mesh_reader.SetFileName('/playpen/workspace/my_paper/linking/data/'+ obj_folder + '/ellipsoid.vtk')
#    ell_mesh_reader.SetFileName('/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/ellipsoid.vtk')
    ell_mesh_reader.Update()
    ell_mesh = ell_mesh_reader.GetOutput()

    ell_srep_reader = vtk.vtkPolyDataReader()
    ell_srep_reader.SetFileName('/playpen/workspace/my_paper/linking/data/'+ obj_folder + '/ell_srep.vtk')
#    ell_srep_reader.SetFileName('/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/ell_srep.vtk')
    ell_srep_reader.Update()
    ell_srep = ell_srep_reader.GetOutput()

    target_mesh_reader = vtk.vtkPolyDataReader()
    target_mesh_reader.SetFileName(target_surf_file)
    target_mesh_reader.Update()
    target_mesh = target_mesh_reader.GetOutput()
    target_skeleton_vtk, initial_spokes_vtk = vsf.compute_vtk_skeleton(target_srep_xml)
#    target_spokes_vtk = sf.refine_srep(initial_spokes_vtk, target_mesh)
    ell_skeleton_vtk, ell_spoke_vtk = vsf.compute_vtk_skeleton("", ell_srep)
    
    p = pv.Plotter()
    p.add_mesh(ell_skeleton_vtk, opacity=0.2, color="grey")
    p.add_mesh(ell_srep, opacity=0.2, color="white")
    p.add_mesh(target_skeleton_vtk, opacity=0.2)
    p.add_mesh(target_mesh, opacity=0.2, color="white")
    obj_plt = pv.Plotter()
    obj_plt.add_mesh(target_mesh, color='white', opacity=0.2)
    obj_plt.add_mesh(initial_spokes_vtk)

    if input_skeleton_tps is None:
        skeleton_tps = vtk.vtkThinPlateSplineTransform()
        skeleton_tps.SetSourceLandmarks(ell_skeleton_vtk.GetPoints())
        skeleton_tps.SetTargetLandmarks(target_skeleton_vtk.GetPoints())
        skeleton_tps.SetBasisToR()
        skeleton_tps.Modified()
    else:
        skeleton_tps = input_skeleton_tps

    tangent_vector_taus = []
    tangent_vector_thetas = []
    normal_vectors = []
    selected_centers = []
    selected_spoke_ids = []
    num_steps = 2
    num_crest_pts = 24
    total_skeletal_pts = ell_skeleton_vtk.GetNumberOfPoints()
    num_interior_pts = (total_skeletal_pts - num_crest_pts) // 2
    selected_rows = [0, 12] #list(range(0, 24))#[0, 12]

    selected_cols = [0] #list(range(0, num_steps + 1)) # [0]
    for i in range(total_skeletal_pts):
        curr_row = i // (num_steps + 1)
        curr_col = i - curr_row * (num_steps + 1)
        curr_pt =  np.array(ell_skeleton_vtk.GetPoint(i))
        next_pt = np.array(ell_skeleton_vtk.GetPoint(i+1))
        pt_next_radial_line = np.array(ell_skeleton_vtk.GetPoint(i + num_steps + 1))
        pt_prev_radial_line = np.array(ell_skeleton_vtk.GetPoint(i - num_steps - 1))

        if i >= total_skeletal_pts - num_crest_pts:# and (i-total_skeletal_pts+num_crest_pts) in selected_crest_spokes:
            ### Frames on the fold curve, tau direction is decided by the interior neighbor point
            i_neighbor = (i - 2*num_interior_pts) * (num_steps + 1)
            interior_neighbor = np.array(ell_skeleton_vtk.GetPoint(i_neighbor))
            ## the last crest point connect to the first crest point
            if i == total_skeletal_pts - 1:
                next_pt = np.array(ell_skeleton_vtk.GetPoint(2 * num_interior_pts))
            tau_vec, new_theta_vec, n_vec = form_basis_vectors(curr_pt, interior_neighbor,next_pt, flip_tau=True)
        elif i < num_interior_pts:
            ## Frames on one side
            if curr_col == num_steps:
                ## Frames at outmost interior points
                i_neighbor = curr_row + 2 * num_interior_pts
                next_pt = np.array(ell_skeleton_vtk.GetPoint(i_neighbor))
            if i == num_interior_pts-1:
                ## the outmost point on the last radial line connect to the first radial line
                pt_next_radial_line = np.array(ell_skeleton_vtk.GetPoint(num_steps))
            elif curr_row in selected_rows and curr_col == 0:
                ## Because the selected rows (spines' ends) have colinear neighbors, choose different neighbors
                pt_next_radial_line = np.array(ell_skeleton_vtk.GetPoint(i + num_steps + 2))
            tau_vec, new_theta_vec, n_vec = form_basis_vectors(curr_pt, next_pt, pt_next_radial_line)
        else:
            ## the other side on the skeleton, the direction of increasing theta is also counter-clock-wise
            if curr_row == num_crest_pts:
                ## the first radial line point to the last radial line
                pt_prev_radial_line = np.array(ell_skeleton_vtk.GetPoint(total_skeletal_pts - 2))
            elif (curr_row - num_crest_pts) in selected_rows and curr_col == 0:
                ## colinear neighbors for points near the spine's ends
                pt_prev_radial_line = np.array(ell_skeleton_vtk.GetPoint(i - num_steps))

            if curr_col == num_steps:
                ## outmost point
                i_neighbor = curr_row - num_crest_pts + 2 * num_interior_pts
                next_pt = np.array(ell_skeleton_vtk.GetPoint(i_neighbor))
            tau_vec, new_theta_vec, n_vec = form_basis_vectors(curr_pt, next_pt, pt_prev_radial_line)

        selected_centers.append(curr_pt)
        tangent_vector_taus.append(tau_vec)
        tangent_vector_thetas.append(new_theta_vec)
        normal_vectors.append(n_vec)
        vec = pv.Arrow(curr_pt, tau_vec)
        normal_arrow = pv.Arrow(curr_pt, n_vec)
        new_theta_arr = pv.Arrow(curr_pt, new_theta_vec)

        selected_spoke_ids.append(i)
        # p.add_mesh(vec, color="red")
        # p.add_mesh(new_theta_arr, color="white")
        # p.add_mesh(normal_arrow, color="blue")
        # base_pt = np.array(ell_spoke_vtk.GetPoint(i * 2))
        # bdry_pt = np.array(ell_spoke_vtk.GetPoint(i * 2 + 1))
        # p.add_mesh(viz.form_spoke_poly(base_pt, bdry_pt), line_width=2, color='red')

#        if i > num_interior_pts: break

    tangent_vector_taus = np.array(tangent_vector_taus)
    tangent_vector_thetas = np.array(tangent_vector_thetas)
    normal_vectors = np.array(normal_vectors)
    centers = np.array(selected_centers)

    num_frames, _ = centers.shape
    ## frames on a onion skin
    target_spokes_vtk = initial_spokes_vtk

    vtk_base_pts1 = vtk.vtkPoints()
    vtk_base_pts2 = vtk.vtkPoints()
    vtk_bdry_pts1 = vtk.vtkPoints()
    vtk_bdry_pts2 = vtk.vtkPoints()
    level_pts = []
    for i_spoke in range(num_frames):
        base_pt_id = i_spoke * 2
        bdry_pt_id = base_pt_id + 1

        base_pt = np.array(target_spokes_vtk.GetPoint(base_pt_id))
        bdry_pt = np.array(target_spokes_vtk.GetPoint(bdry_pt_id))
        spoke_dir = bdry_pt - base_pt
        radius = np.linalg.norm(spoke_dir)
        spoke_dir /= radius
        # if i_spoke < 2:
        #     p.add_mesh(viz.form_spoke_poly(base_pt, bdry_pt), line_width=2, color='red')

        tau_spoke = 0.2 * alpha
        bdry_pt = base_pt + tau_spoke * radius * spoke_dir
        level_pts.append(bdry_pt)
        level_pts.append(base_pt)
        if i_spoke < num_interior_pts or i_spoke >= 2 * num_interior_pts:
            ## one side of object
            vtk_base_pts1.InsertNextPoint(base_pt)
            vtk_bdry_pts1.InsertNextPoint(bdry_pt)
        else:
            ## the other side of object
            vtk_base_pts2.InsertNextPoint(base_pt)
            vtk_bdry_pts2.InsertNextPoint(bdry_pt)

    level_poly = pv.PolyData(np.array(level_pts))
    obj_plt.add_mesh(level_poly)

    skeleton_bdry_tps1 = vtk.vtkThinPlateSplineTransform()
    skeleton_bdry_tps1.SetSourceLandmarks(vtk_base_pts1)
    skeleton_bdry_tps1.SetTargetLandmarks(vtk_bdry_pts1)
    skeleton_bdry_tps1.Update()

    skeleton_bdry_tps2 = vtk.vtkThinPlateSplineTransform()
    skeleton_bdry_tps2.SetSourceLandmarks(vtk_base_pts2)
    skeleton_bdry_tps2.SetTargetLandmarks(vtk_bdry_pts2)
    skeleton_bdry_tps2.Update()

    ell_onion_skin_tau_ends = []
    ell_onion_skin_theta_ends = []
    ell_onion_skin_centers = []
    obj_onion_skin_tau_ends = []
    obj_onion_skin_theta_ends = []
    obj_onion_skin_centers = []
    obj_skeletal_frames = []
    obj_onion_frames = []
    deformed_centers_via_tps = []
    ### the top side that normals point to
    for i in range(num_frames):
        if i < num_interior_pts or i >= 2 * num_interior_pts:
            ## one side of object
            skeleton_bdry_tps = skeleton_bdry_tps1
        else:
            skeleton_bdry_tps = skeleton_bdry_tps2

        ## Deform skeletal frames to the object
        # obj_skeletal_center     = np.array(skeleton_tps.TransformPoint(centers[i, :]))
        # obj_skeletal_tau_end    = obj_skeletal_center + np.array(skeleton_tps.TransformVectorAtPoint(centers[i, :], alpha * delta*tangent_vector_taus[i, :]))
        # obj_skeletal_theta_end  = obj_skeletal_center + np.array(skeleton_tps.TransformVectorAtPoint(centers[i, :], alpha * delta*tangent_vector_thetas[i, :]))
        # obj_skeletal_normal_end = obj_skeletal_center + np.array(skeleton_tps.TransformVectorAtPoint(centers[i, :], alpha * delta*normal_vectors[i, :]))

        obj_skeletal_tau_end    = np.array(skeleton_tps.TransformPoint(centers[i, :] + delta*tangent_vector_taus[i, :]))
        obj_skeletal_theta_end  = np.array(skeleton_tps.TransformPoint(centers[i, :] + delta*tangent_vector_thetas[i, :]))
        obj_skeletal_normal_end = np.array(skeleton_tps.TransformPoint(centers[i, :] + delta*normal_vectors[i, :]))
        obj_skeletal_center     = np.array(skeleton_tps.TransformPoint(centers[i, :]))
        obj_skeletal_frames.append((obj_skeletal_center, obj_skeletal_tau_end, obj_skeletal_theta_end, obj_skeletal_normal_end))


        ### Deform ellipsoidal skeleton to object's skeleton, Use this deformation to deform the onion skin near the skeleton
        obj_onion_tau_end      = np.array(skeleton_bdry_tps.TransformPoint(obj_skeletal_tau_end))
        obj_onion_theta_end    = np.array(skeleton_bdry_tps.TransformPoint(obj_skeletal_theta_end))
        obj_onion_normal_end   = np.array(skeleton_bdry_tps.TransformPoint(obj_skeletal_normal_end))
        obj_onion_frame_center = np.array(skeleton_bdry_tps.TransformPoint(obj_skeletal_center))

        obj_onion_frames.append((obj_onion_frame_center, obj_onion_tau_end, obj_onion_theta_end, obj_onion_normal_end))
        deformed_centers_via_tps.append(obj_onion_frame_center)
        # tau_vec, new_theta_vec, n_vec = form_basis_vectors(obj_skeletal_center, obj_skeletal_tau_end, obj_skeletal_theta_end)
        # vec = pv.Arrow(obj_skeletal_center, tau_vec)
        # normal_arrow = pv.Arrow(obj_skeletal_center, n_vec)
        # new_theta_arr = pv.Arrow(obj_skeletal_center, new_theta_vec)
        # obj_plt.add_mesh(vec, color='red')
        # obj_plt.add_mesh(new_theta_arr, color="white")
        # obj_plt.add_mesh(normal_arrow, color='blue')
        # obj_plt.add_mesh(pv.PolyData(obj_skeletal_normal_end))


        # tau_vec, new_theta_vec, n_vec = form_basis_vectors(obj_onion_frame_center, obj_onion_tau_end, obj_onion_theta_end)
        # vec = pv.Arrow(obj_onion_frame_center, tau_vec)
        # normal_arrow = pv.Arrow(obj_onion_frame_center, n_vec)
        # new_theta_arr = pv.Arrow(obj_onion_frame_center, new_theta_vec)

        # obj_plt.add_mesh(vec, color='red')
        # obj_plt.add_mesh(new_theta_arr, color="white")
        # obj_plt.add_mesh(normal_arrow, color='blue')        

#    obj_plt.show()
    p.add_mesh(pv.PolyData(np.array(deformed_centers_via_tps)))
#    p.show()

    lines = [] # connections between skeletal points and corresponding boundary points

    within_feats = []

    for i in range(num_frames):
        ## add top boundary points with skeletal points
        skeletal_pts = obj_skeletal_frames[i][0]
        bdry_pts1 = obj_onion_frames[i][0]
        line1 = viz.form_spoke_poly(skeletal_pts, bdry_pts1)
        lines.append(line1)

        np_line1 = bdry_pts1 - skeletal_pts
        line1_length = np.linalg.norm(np_line1)
        dir_feats = map_dir_to_local_frame(np_line1, np.array(obj_skeletal_frames[i]))
        dir_feats2 = map_dir_to_local_frame(np_line1, np.array(obj_onion_frames[i]))
        within_feats.append(np.concatenate((skeletal_pts, dir_feats, bdry_pts1, dir_feats2, [line1_length])))

    within_feats = np.array(within_feats)

    bdry_pts_vec = []
    for i in range(target_mesh.GetNumberOfPoints()):
        np_pt = np.array(target_mesh.GetPoint(i))
        bdry_pts_vec.append(np_pt)
    bdry_pts_vec = np.array(bdry_pts_vec).flatten()
    return np.array(lines), np.array(obj_skeletal_frames), within_feats, np.array(obj_onion_frames), bdry_pts_vec, skeleton_tps

# def special_lines_local_frames(obj_folder, input_id):
#     target_srep_xml = '/playpen/workspace/my_paper/linking/data/' + obj_folder + '/' + input_id + '/header.xml'
#     target_surf_file = '/playpen/workspace/my_paper/linking/data/' + obj_folder + '/' + input_id + '/stx_noscale_' + input_id + '_V06_t1w_RAI_Bias_label_pp_surfSPHARM.vtk'
#     ## initialize s-rep
#     # print("Creating s-rep for ", input_id)
#     # output_path = '/playpen/workspace/my_paper/linking/data/' + obj_folder + '/' + input_id
#     # os.system("Slicer --no-main-window --python-script /playpen/workspace/Simulate_Shapes/shanapy/models/cmd_initialize_srep.py " + target_surf_file + " " + output_path)

#     ## special frames on ellipsoid
#     ell_mesh_reader = vtk.vtkPolyDataReader()
#     ell_mesh_reader.SetFileName('/playpen/workspace/my_paper/linking/data/'+ obj_folder + '/ellipsoid.vtk')
#     ell_mesh_reader.Update()
#     ell_mesh = ell_mesh_reader.GetOutput()

#     ell_srep_reader = vtk.vtkPolyDataReader()
#     ell_srep_reader.SetFileName('/playpen/workspace/my_paper/linking/data/'+ obj_folder + '/ell_srep.vtk')
#     ell_srep_reader.Update()
#     ell_srep = ell_srep_reader.GetOutput()

#     target_mesh_reader = vtk.vtkPolyDataReader()
#     target_mesh_reader.SetFileName(target_surf_file)
#     target_mesh_reader.Update()
#     target_mesh = target_mesh_reader.GetOutput()
#     target_skeleton_vtk, initial_spokes_vtk = vsf.compute_vtk_skeleton(target_srep_xml)
#     target_spokes_vtk = sf.refine_srep(initial_spokes_vtk, target_mesh)
#     ell_skeleton_vtk, ell_spoke_vtk = vsf.compute_vtk_skeleton("", ell_srep)

#     # srep_plotter = pv.Plotter()
#     # srep_plotter.add_mesh(target_mesh, opacity=0.3, color='white')
#     # srep_plotter.add_mesh(target_spokes_vtk)
#     # srep_plotter.show()
#     ## deformation between two skeletons
#     skeleton_tps = vtk.vtkThinPlateSplineTransform()
#     skeleton_tps.SetSourceLandmarks(ell_skeleton_vtk.GetPoints())
#     skeleton_tps.SetTargetLandmarks(target_skeleton_vtk.GetPoints())
#     skeleton_tps.SetBasisToR()
#     skeleton_tps.Modified()

#     tangent_vector_taus = []
#     tangent_vector_thetas = []
#     normal_vectors = []
#     selected_centers = []
#     selected_spoke_ids = []
#     num_steps = 2
#     num_crest_pts = 24
#     total_skeletal_pts = ell_skeleton_vtk.GetNumberOfPoints()
#     num_interior_pts = (total_skeletal_pts - num_crest_pts) // 2
#     selected_rows = [0, 12] #list(range(0, 24))#[0, 12]
#     selected_crest_spokes = [24//4, 24//4 * 3] #list(range(0, 24)) #[24//4, 24//4 * 3]
#     selected_cols = [0] #list(range(0, num_steps + 1)) # [0]

#     p = pv.Plotter()
#     p.add_mesh(ell_skeleton_vtk, color="grey")
#     p.add_mesh(ell_srep)
#     ### Select special frames from ellipsoid
#     ### Tau, theta and normal are three basis vectors that follow right hand rule
#     for i in range(total_skeletal_pts):
#         curr_row = i // (num_steps + 1)
#         curr_col = i - curr_row * (num_steps + 1)
#         curr_pt =  np.array(ell_skeleton_vtk.GetPoint(i))
#         next_pt = np.array(ell_skeleton_vtk.GetPoint(i+1))
#         pt_next_radial_line = np.array(ell_skeleton_vtk.GetPoint(i + num_steps + 1))
#         pt_prev_radial_line = np.array(ell_skeleton_vtk.GetPoint(i - num_steps - 1))

#         if i >= total_skeletal_pts - num_crest_pts:# and (i-total_skeletal_pts+num_crest_pts) in selected_crest_spokes:
#             ### Frames on the fold curve, tau direction is decided by the interior neighbor point
#             i_neighbor = (i - 2*num_interior_pts) * (num_steps + 1)
#             interior_neighbor = np.array(ell_skeleton_vtk.GetPoint(i_neighbor))

#             ## the last crest point connect to the first crest point
#             if i == total_skeletal_pts - 1:
#                 next_pt = np.array(ell_skeleton_vtk.GetPoint(2 * num_interior_pts))
#             tau_vec, new_theta_vec, n_vec = form_basis_vectors(curr_pt, interior_neighbor,next_pt, flip_tau=True)
#         elif i < num_interior_pts:
#             ## Frames on one side
#             if curr_col == num_steps:
#                 ## Frames at outmost interior points
#                 i_neighbor = curr_row + 2 * num_interior_pts
#                 next_pt = np.array(ell_skeleton_vtk.GetPoint(i_neighbor))
#             if i == num_interior_pts-1:
#                 ## the outmost point on the last radial line connect to the first radial line
#                 pt_next_radial_line = np.array(ell_skeleton_vtk.GetPoint(num_steps))
#             elif curr_row in selected_rows and curr_col == 0:
#                 ## Because the selected rows (spines' ends) have colinear neighbors, choose different neighbors
#                 pt_next_radial_line = np.array(ell_skeleton_vtk.GetPoint(i + num_steps + 2))
#             tau_vec, new_theta_vec, n_vec = form_basis_vectors(curr_pt, next_pt, pt_next_radial_line)
#         else:
#             ## the other side on the skeleton, the direction of increasing theta is also counter-clock-wise
#             if curr_row == num_crest_pts:
#                 ## the first radial line point to the last radial line
#                 pt_prev_radial_line = np.array(ell_skeleton_vtk.GetPoint(total_skeletal_pts - 2))
#             elif (curr_row - num_crest_pts) in selected_rows and curr_col == 0:
#                 ## colinear neighbors for points near the spine's ends
#                 pt_prev_radial_line = np.array(ell_skeleton_vtk.GetPoint(i - num_steps))

#             if curr_col == num_steps:
#                 ## outmost point
#                 i_neighbor = curr_row - num_crest_pts + 2 * num_interior_pts
#                 next_pt = np.array(ell_skeleton_vtk.GetPoint(i_neighbor))
#             tau_vec, new_theta_vec, n_vec = form_basis_vectors(curr_pt, next_pt, pt_prev_radial_line)

#         selected_centers.append(curr_pt)
#         tangent_vector_taus.append(tau_vec)
#         tangent_vector_thetas.append(new_theta_vec)
#         normal_vectors.append(n_vec)
#         vec = pv.Arrow(curr_pt, tau_vec)
#         normal_arrow = pv.Arrow(curr_pt, n_vec)
#         new_theta_arr = pv.Arrow(curr_pt, new_theta_vec)

#         selected_spoke_ids.append(i)
#         p.add_mesh(vec, color="red")
#         p.add_mesh(new_theta_arr, color="white")
#         p.add_mesh(normal_arrow, color="blue")

#     #p.show()
#     tangent_vector_taus = np.array(tangent_vector_taus)
#     tangent_vector_thetas = np.array(tangent_vector_thetas)
#     normal_vectors = np.array(normal_vectors)
#     centers = np.array(selected_centers)

#     ### Deform frames on skeleton, given the deformation via TPS that smoothly deform tangent vectors
#     ### The deformed tangent vectors may not be orthogonal. However, they are still tangent vectors on
#     ### the deformed skeleton.
#     ### Compute the new normals via cross product according to the tangent vectors.
#     ### Then fix one tangent vector and compute the other othogonal to it via cross product.
#     num_pts, _ = centers.shape
    
#     tau_ends = centers + delta * tangent_vector_taus
#     theta_ends = centers + delta * tangent_vector_thetas
#     normal_ends = centers + delta * normal_vectors
#     skeletal_frames_plot = pv.Plotter()
#     bdry_plot = pv.Plotter()
#     skeletal_frames_plot.add_mesh(target_skeleton_vtk)
#     skeletal_frames_plot.add_mesh(target_mesh, opacity=0.2, color='white')

#     #skeletal_frames_plot.add_mesh(target_spokes_vtk)
#     skeletal_tau_end = []
#     skeletal_theta_end = []
#     skeletal_center_pt = []
#     skeletal_frames = []
#     for i in range(num_pts):
#         curr_tau_end = np.array(skeleton_tps.TransformPoint(tau_ends[i, :]))
#         curr_theta_end = np.array(skeleton_tps.TransformPoint(theta_ends[i, :]))
#         curr_center = np.array(skeleton_tps.TransformPoint(centers[i, :]))

#         # if True:# i == 2:
#         #     pv_sphere = pv.Sphere(0.1, center=curr_tau_end)
#         #     skeletal_frames_plot.add_mesh(pv_sphere)
#         #     pv_center = pv.Sphere(0.1, center=curr_center)
#         #     skeletal_frames_plot.add_mesh(pv_center, color='red')

#         #     pv_sphere = pv.Sphere(0.1, center=tau_ends[i, :])
#         #     pv_center = pv.Sphere(0.1, center=centers[i, :])
#         #     p.add_mesh(pv_sphere)
#         #     p.add_mesh(pv_center, color='red')
#         #     p.show()
#         tau_vec, new_theta_vec, n_vec = form_basis_vectors(curr_center, curr_tau_end, curr_theta_end)
#         ## only display the first frame in every object to debug
#         if True:#i == 0:
#             vec = pv.Arrow(curr_center, tau_vec)
#             normal_arrow = pv.Arrow(curr_center, n_vec)
#             new_theta_arr = pv.Arrow(curr_center, new_theta_vec)
#             skeletal_frames_plot.add_mesh(vec, color="red")
#             skeletal_frames_plot.add_mesh(new_theta_arr, color="white")
#             skeletal_frames_plot.add_mesh(normal_arrow, color="blue")
#         skeletal_frames.append((curr_center, curr_tau_end, curr_center + delta * new_theta_vec, curr_center + delta * n_vec))

#         skeletal_tau_end.append(curr_center + delta * tau_vec)
#         skeletal_theta_end.append(curr_center + delta * new_theta_vec)
#         skeletal_center_pt.append(curr_center)

#         # bdry_plot.add_mesh(vec, color="red")
#         # bdry_plot.add_mesh(new_theta_arr, color="white")
#         # bdry_plot.add_mesh(normal_arrow, color="blue")

#     skeletal_tau_end = np.array(skeletal_tau_end)
#     skeletal_theta_end = np.array(skeletal_theta_end)
#     skeletal_center_pt = np.array(skeletal_center_pt)
# #    skeletal_frames_plot.show()

#     ### Local frames on boundary
#     ### Deform frames on the skeleton to the boundary of an object. 

#     bdry_plot.add_mesh(target_mesh, opacity=0.2, color="white")
#     vtk_base_pts1 = vtk.vtkPoints()
#     vtk_base_pts2 = vtk.vtkPoints()
#     vtk_bdry_pts1 = vtk.vtkPoints()
#     vtk_bdry_pts2 = vtk.vtkPoints()
#     for i_spoke in range(num_pts):
#         base_pt_id = i_spoke * 2
#         bdry_pt_id = base_pt_id + 1

#         if i_spoke < num_interior_pts or i_spoke >= 2 * num_interior_pts:
#             ## one side of object
#             vtk_base_pts1.InsertNextPoint(target_spokes_vtk.GetPoint(base_pt_id))
#             vtk_bdry_pts1.InsertNextPoint(target_spokes_vtk.GetPoint(bdry_pt_id))
#         else:
#             ## the other side of object
#             vtk_base_pts2.InsertNextPoint(target_spokes_vtk.GetPoint(base_pt_id))
#             vtk_bdry_pts2.InsertNextPoint(target_spokes_vtk.GetPoint(bdry_pt_id))

#     skeleton_bdry_tps1 = vtk.vtkThinPlateSplineTransform()
#     skeleton_bdry_tps1.SetSourceLandmarks(vtk_base_pts1)
#     skeleton_bdry_tps1.SetTargetLandmarks(vtk_bdry_pts1)
#     skeleton_bdry_tps1.Update()

#     skeleton_bdry_tps2 = vtk.vtkThinPlateSplineTransform()
#     skeleton_bdry_tps2.SetSourceLandmarks(vtk_base_pts2)
#     skeleton_bdry_tps2.SetTargetLandmarks(vtk_bdry_pts2)
#     skeleton_bdry_tps2.Update()

#     bdry_frames = []
#     ### the top side that normals point to
#     for i in range(num_pts):
#         if i < num_interior_pts or i >= 2 * num_interior_pts:
#             ## one side of object
#             skeleton_bdry_tps = skeleton_bdry_tps1
#         else:
#             skeleton_bdry_tps = skeleton_bdry_tps2

#         deformed_tau_end = np.array(skeleton_bdry_tps.TransformPoint(skeletal_tau_end[i, :]))
#         deformed_theta_end = np.array(skeleton_bdry_tps.TransformPoint(skeletal_theta_end[i, :]))
#         deformed_center = np.array(target_spokes_vtk.GetPoint(i * 2 + 1))

#         tau_vec, new_theta_vec, n_vec = form_basis_vectors(deformed_center, deformed_tau_end, deformed_theta_end)
#         vec = pv.Arrow(deformed_center, tau_vec)
#         normal_arrow = pv.Arrow(deformed_center, n_vec)
#         new_theta_arr = pv.Arrow(deformed_center, new_theta_vec)
#         bdry_plot.add_mesh(vec, color='red')
#         bdry_plot.add_mesh(new_theta_arr, color="white")
#         bdry_plot.add_mesh(normal_arrow, color='blue')
#         bdry_frames.append((deformed_center, deformed_tau_end, deformed_center + delta * new_theta_vec, deformed_center + delta * n_vec))

#     # bdry_plot.show()

#     ### Generate within-object features w.r.t. the local frames
#     line_seg_plot = pv.Plotter()
#     line_seg_plot.add_mesh(target_mesh, opacity=0.3, color="white")
#     lines = []
#     within_feats = []
#     extended_within_feats = []
#     for i in range(num_pts):
#         ## add top boundary points with skeletal points
#         skeletal_pts = skeletal_frames[i][0]
#         bdry_pts1 = bdry_frames[i][0]
#         line1 = viz.form_spoke_poly(skeletal_pts, bdry_pts1)
#         line_seg_plot.add_mesh(line1)
#         lines.append(line1)

#         np_line1 = bdry_pts1 - skeletal_pts
#         line1_length = np.linalg.norm(np_line1)
#         dir_feats = map_dir_to_local_frame(np_line1, np.array(skeletal_frames[i]))
#         dir_feats2 = map_dir_to_local_frame(np_line1, np.array(bdry_frames[i]))
#         within_feats.append(np.concatenate((skeletal_pts, dir_feats, bdry_pts1, dir_feats2, [line1_length])))

# #    line_seg_plot.show()
#     within_feats = np.array(within_feats)

#     return lines, skeletal_frames, within_feats, bdry_frames
def map_dir_to_local_frame(u_dir, frame_pts):
    """
    Input u_dir: a unit vector in global c.s.
    Input frame_pts 4x3 matrix: 4 frame points that form a 3 basis vectors of a frame
    Output 3-tuple
    """
    vec_u = frame_pts[1, :] - frame_pts[0, :]
    vec_u = vec_u / np.linalg.norm(vec_u)

    vec_t = frame_pts[2, :] - frame_pts[0, :]
    vec_t = vec_t / np.linalg.norm(vec_t)

    vec_n = frame_pts[3, :] - frame_pts[0, :]
    vec_n = vec_n / np.linalg.norm(vec_n)

    u_dir /= np.linalg.norm(u_dir)
    return np.array([np.dot(u_dir, vec_u), np.dot(u_dir, vec_t), np.dot(u_dir, vec_n)])

def between_features_local_frames(hipp_skeletal_frames, caud_skeletal_frames):
    """
    Each skeletal frames consists of 5 selected skeletal frames: 3 on spine, 2 on the orthogonal direction
    """
    link_plotter = pv.Plotter()
    ret_feats = []
    for i in range(len(hipp_skeletal_frames)):
        hipp_frame_origin = hipp_skeletal_frames[i][0]
        caud_frame_origin = caud_skeletal_frames[i][0]

        line = caud_frame_origin - hipp_frame_origin
        line_length = np.linalg.norm(line)
        line_dir = line / line_length

        hipp_proj = map_dir_to_local_frame(line_dir, np.array(hipp_skeletal_frames[i]))
        caud_proj = map_dir_to_local_frame(line_dir, np.array(caud_skeletal_frames[i]))
        dir_feats = np.concatenate((hipp_frame_origin, hipp_proj, caud_frame_origin, caud_proj))
        all_feats = np.append(dir_feats, [line_length])
        ret_feats.append(all_feats)

    return np.array(ret_feats)
def commensurate_length(feats, num_feats_per_vec):
    print("Commensurate lengths...")
    n, d = feats.shape
    num_vectors = d // num_feats_per_vec
    lengths = []
    for i in range(num_vectors):
        id_curr_len = (i + 1) * num_feats_per_vec - 1
        log_feats = np.log(feats[:, id_curr_len])
        geo_mean = np.exp(np.mean(log_feats))
        feats[:, id_curr_len] = (log_feats - np.mean(log_feats))
        lengths.append(log_feats - np.mean(log_feats))
    return feats, np.array(lengths).T
def commensurate_length_at_id(feats, id_curr_len):
    print("Commensurate lengths...")
    n, d = feats.shape

    log_feats = np.log(feats[:, id_curr_len])
    geo_mean = np.exp(np.mean(log_feats))
    feats[:, id_curr_len] = (log_feats - np.mean(log_feats))

    return feats

### More recent implementation is in local_frames_spike_n_slab.py
# def load_features(regenerate_feats=False, align_pts=False):
#     print("Load features associated with the special local frames...")
#     intra_feats = []
#     pos_inter_feats = []
#     neg_inter_feats = []

#     pos_hipp_feats = []
#     pos_caud_feats = []
#     neg_hipp_feats = []
#     neg_caud_feats = []
#     pos_hipp_line_tmp = None
#     pos_caud_line_tmp = None
#     all_pos_hipp_skeletal_frames = []
#     all_pos_caud_skeletal_frames = []
#     all_neg_hipp_skeletal_frames = []
#     all_neg_caud_skeletal_frames = []

#     all_pos_hipp_bdry_frames = []
#     all_pos_caud_bdry_frames = []
#     all_neg_hipp_bdry_frames = []
#     all_neg_caud_bdry_frames = []

#     if regenerate_feats:
#         for i, pos_id in enumerate(pos_ids):
#             print(i, pos_id)
#             pos_hipp_line, pos_hipp_skeletal_frames, pos_hipp_within_feats, pos_hipp_bdry_frames = lfp.special_lines_local_frames('nonaligned_hipp_sreps', pos_id)
#             pos_caud_line, pos_caud_skeletal_frames, pos_caud_within_feats, pos_caud_bdry_frames = lfp.special_lines_local_frames('nonaligned_caud_sreps', pos_id)
#             all_pos_hipp_skeletal_frames.append(pos_hipp_skeletal_frames)
#             all_pos_caud_skeletal_frames.append(pos_caud_skeletal_frames)
#             all_pos_hipp_bdry_frames.append(pos_hipp_bdry_frames)
#             all_pos_caud_bdry_frames.append(pos_caud_bdry_frames)

#             # show_pts_on_sphere(pos_hipp_within_feats[:, :3])
#             # show_pts_on_sphere(pos_caud_within_feats[:, :3])
#             # if pos_hipp_line_tmp is None:
#             #     pos_hipp_line_tmp = hipp_within_frames[0][1] - hipp_within_frames[0][0]
#             #     pos_caud_line_tmp = caud_within_frames[0][1] - caud_within_frames[0][0]
#             #     pos_hipp_line_tmp /= np.linalg.norm(pos_hipp_line_tmp)
#             #     pos_caud_line_tmp /= np.linalg.norm(pos_caud_line_tmp)
#             # else:
#             #     new_theta_dir = hipp_within_frames[0][1] - hipp_within_frames[0][0]
#             #     new_theta_dir /= np.linalg.norm(new_theta_dir)
#             #     angle = np.dot(pos_hipp_line_tmp, new_theta_dir)
#             #     print('Positive angle: ', angle)

#             pos_inter = lfp.between_features_local_frames(pos_hipp_skeletal_frames, pos_caud_skeletal_frames)
#             pos_hipp_feats.append(pos_hipp_within_feats.flatten())
#             pos_caud_feats.append(pos_caud_within_feats.flatten())
#             # pos_hipp_feats.append(pos_hipp_within_ext.flatten())
#             # pos_caud_feats.append(pos_caud_within_ext.flatten())

#             pos_inter_feats.append(pos_inter.flatten())

#         for i, neg_id in enumerate(neg_ids):
#             neg_id = neg_ids[i]
# #            if neg_id == "995004": neg_id = "974562" # replace outlier with another case
#             print(neg_id)
#             neg_hipp_lines, neg_hipp_skeletal_frames, neg_hipp_within_feats, neg_hipp_bdry_frames = lfp.special_lines_local_frames('nonaligned_hipp_sreps', neg_id)
#             neg_caud_lines, neg_caud_skeletal_frames, neg_caud_within_feats, neg_caud_bdry_frames = lfp.special_lines_local_frames('nonaligned_caud_sreps', neg_id)
#             all_neg_hipp_skeletal_frames.append(neg_hipp_skeletal_frames)
#             all_neg_caud_skeletal_frames.append(neg_caud_skeletal_frames)
#             all_neg_hipp_bdry_frames.append(neg_hipp_bdry_frames)
#             all_neg_caud_bdry_frames.append(neg_caud_bdry_frames)

#             # new_theta_dir = hipp_within_frames[0][1] - hipp_within_frames[0][0]
#             # new_theta_dir /= np.linalg.norm(new_theta_dir)
#             # angle = np.dot(pos_hipp_line_tmp, new_theta_dir)
#             # print('Negative angle: ', angle)

#             neg_inter = lfp.between_features_local_frames(neg_hipp_skeletal_frames, neg_caud_skeletal_frames)
#             #neg_intra_feats = np.concatenate((neg_hipp_within_feats, neg_caud_within_feats))

#             neg_hipp_feats.append(neg_hipp_within_feats.flatten())
#             neg_caud_feats.append(neg_caud_within_feats.flatten())
#             # neg_hipp_feats.append(neg_hipp_within_ext.flatten())
#             # neg_caud_feats.append(neg_caud_within_ext.flatten())

#             neg_inter_feats.append(neg_inter.flatten())

#         inter_feats = np.concatenate((np.array(pos_inter_feats), np.array(neg_inter_feats)))
#         hipp_feats  = np.concatenate((np.array(pos_hipp_feats),   np.array(neg_hipp_feats)))
#         caud_feats  = np.concatenate((np.array(pos_caud_feats),   np.array(neg_caud_feats)))
#         hipp_frames = np.array([all_pos_hipp_skeletal_frames, all_neg_hipp_skeletal_frames, all_pos_hipp_bdry_frames, all_neg_hipp_bdry_frames])
#         caud_frames = np.array([all_pos_caud_skeletal_frames, all_neg_caud_skeletal_frames, all_pos_caud_bdry_frames, all_neg_caud_bdry_frames])
#         # if os.path.exists('/playpen/workspace/Simulate_Shapes/data/local_frames_feats.npy'):
#         #     os.system('rm /playpen/workspace/Simulate_Shapes/data/local_frames_feats.npy')
#         # with open('/playpen/workspace/Simulate_Shapes/data/local_frames_feats.npy', 'wb') as f:
#         #     np.save(f, intra_feats)
#         #     np.save(f, inter_feats)
#         #     np.save(f, hipp_feats)
#         #     np.save(f, caud_feats)
#         #     np.save(f, hipp_frames)
#         #     np.save(f, caud_frames)
#     else:
#         with open('/playpen/workspace/Simulate_Shapes/data/local_frames_feats.npy', 'rb') as f:
#             intra_feats = np.load(f)
#             inter_feats = np.load(f)
#             hipp_feats = np.load(f)
#             caud_feats = np.load(f)
#             hipp_frames = np.load(f)
#             caud_frames = np.load(f)
#     ## Euclideanize and commensurate data
#     # euc_hipp_feats, aligned_pts_hipp = euclideanize_feats(hipp_feats, 13, align_pts)
#     # euc_caud_feats, aligned_pts_caud = euclideanize_feats(caud_feats, 13, align_pts)

#     # # dir_hipp_feats, len_hipp_feats = separate_feats(euc_hipp_feats)
#     # # dir_caud_feats, len_caud_feats = separate_feats(euc_caud_feats)

#     # euc_intra_feats = np.concatenate((euc_hipp_feats, euc_caud_feats), axis=1)
#     # euc_inter_feats, aligned_pts_inter = euclideanize_feats(inter_feats, 13, align_pts)

#     # intra_feats = commensurate_length(euc_intra_feats, 13)
#     # inter_feats = commensurate_length(euc_inter_feats, 13)

#     # ## extend features with aligned pts
#     # intra_feats = np.concatenate((intra_feats, aligned_pts_hipp, aligned_pts_caud), axis=1)
#     # inter_feats = np.concatenate((inter_feats, aligned_pts_inter), axis=1)

#     return intra_feats, inter_feats, hipp_feats, caud_feats, hipp_frames, caud_frames

if __name__ == '__main__':
    hipp_within_lines, hipp_within_frames, hipp_within_feats, _ = special_lines_local_frames('nonaligned_hipp_sreps', '553295')
    caud_within_lines, caud_within_frames, caud_within_feats, _ = special_lines_local_frames('nonaligned_caud_sreps', '553295')
    between_feats = between_features_local_frames(hipp_within_frames, caud_within_frames)
    print(between_feats.shape)