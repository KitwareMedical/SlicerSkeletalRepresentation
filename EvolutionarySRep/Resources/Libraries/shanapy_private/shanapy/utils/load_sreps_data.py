from shanapy.utils import obtain_class_labels as ocl, viz, procrustes
from shanapy.models import srep_fitter
import vtk
import os, glob
import numpy as np
from scipy import spatial

pos_ids, neg_ids = ocl.load_class_labels()

def connect_spokes(spoke_poly):
    lines = vtk.vtkCellArray()
    for i in range(0, spoke_poly.GetNumberOfPoints(), 2):
        id0 = i
        id1 = i+1

        spoke_line = vtk.vtkLine()
        spoke_line.GetPointIds().SetId(0, id0)
        spoke_line.GetPointIds().SetId(1, id1)
        lines.InsertNextCell(spoke_line)
    spoke_poly.SetLines(lines)
    return spoke_poly

def load_srep(input_folder_name):
    up_file_name = os.path.join(input_folder_name, 'up.vtp')
    up_reader = vtk.vtkXMLPolyDataReader()
    up_reader.SetFileName(up_file_name)
    up_reader.Update()
    up_spokes = up_reader.GetOutput()
#    up_spokes = connect_spokes(up_spokes)

    down_file_name = os.path.join(input_folder_name, 'down.vtp')
    down_reader = vtk.vtkXMLPolyDataReader()
    down_reader.SetFileName(down_file_name)
    down_reader.Update()
    down_spokes = down_reader.GetOutput()
#    down_spokes = connect_spokes(down_spokes)
    
    crest_file_name = os.path.join(input_folder_name, 'crest.vtp')
    crest_reader = vtk.vtkXMLPolyDataReader()
    crest_reader.SetFileName(crest_file_name)
    crest_reader.Update()
    crest_spokes = crest_reader.GetOutput()
#    crest_spokes = connect_spokes(crest_spokes)

    return up_spokes, down_spokes, crest_spokes

def obtain_features(spokes_poly):
    directions = [] # k x 3
    radii = []
    bdry_pts = [] # k x 3
    skeletal_pts = [] # k x 3
    for i in range(0, spokes_poly.GetNumberOfPoints(), 2):
        id0 = i
        id1 = i+1

        base_pt = np.array(spokes_poly.GetPoint(id0))
        bdry_pt = np.array(spokes_poly.GetPoint(id1))

        diff = bdry_pt - base_pt
        directions.append((diff) / np.linalg.norm(diff))
        radii.append(np.linalg.norm(diff))
        bdry_pts.append(bdry_pt)
        skeletal_pts.append(base_pt)
    return np.array(directions), np.array(radii), np.array(bdry_pts), np.array(skeletal_pts)
def obtain_bdry_skeletal_pts(up_spokes, down_spokes, crest_spokes, surface_mesh):
    """
    Input spokes are skeletal points and spoke lengths and dirs
    """
    appender = vtk.vtkAppendPolyData()
    appender.AddInputData(up_spokes)
    appender.AddInputData(down_spokes)
    appender.AddInputData(crest_spokes)
    appender.Update()
    new_spokes = appender.GetOutput()

    spoke_poly = vtk.vtkPolyData()
    spoke_pts = vtk.vtkPoints()
    spoke_lines = vtk.vtkCellArray()

    spoke_poly_appender = vtk.vtkAppendPolyData()

    skeletal_pts = []
    bdry_pts = []
    for i in range(new_spokes.GetNumberOfPoints()):
        base_pt = np.array(new_spokes.GetPoint(i))
        arr = new_spokes.GetPointData()
        radii = arr.GetArray("spokeLength").GetValue(i)

        spoke_dir = np.array(arr.GetArray("spokeDirection").GetTuple3(i))
        bdry_pt = base_pt + radii * spoke_dir
        spoke_poly = viz.form_spoke_poly(base_pt, bdry_pt)
        spoke_poly_appender.AddInputData(spoke_poly)

        skeletal_pts.append(base_pt)
        bdry_pts.append(bdry_pt)
    spoke_poly_appender.Update()
    
    initial_spokes = spoke_poly_appender.GetOutput()
#    viz.overlay_polydata(initial_spokes, surface_mesh)
    refined_spokes = srep_fitter.refine_srep(initial_spokes, surface_mesh)
#    viz.overlay_polydata(refined_spokes, surface_mesh)
    _, _, all_bdry_pts, all_skeletal_pts = obtain_features(refined_spokes)
    return np.array(all_skeletal_pts), np.array(all_bdry_pts), refined_spokes

def obtain_all_bdry_pts(up_spokes, down_spokes, crest_spokes, surface_mesh):
    # _, _, up_bdry, up_skeletal = obtain_features(up_spokes)
    # _, _, down_bdry, down_skeletal = obtain_features(down_spokes)
    # _, _, crest_bdry, crest_skeletal = obtain_features(crest_spokes)

    # all_bdry_pts = np.concatenate((up_bdry, down_bdry, crest_bdry))
    # all_skeletal_pts = np.concatenate((up_skeletal, down_skeletal, crest_skeletal))
    appender = vtk.vtkAppendPolyData()
    appender.AddInputData(up_spokes)
    appender.AddInputData(down_spokes)
    appender.AddInputData(crest_spokes)
    appender.Update()
    new_spokes = appender.GetOutput()
    refined_spokes = srep_fitter.refine_srep(new_spokes, surface_mesh)
    _, _, all_bdry_pts, all_skeletal_pts = obtain_features(refined_spokes)
    return all_bdry_pts, all_skeletal_pts
def procrustes_alignment(template, moving_shape):
    mtx1, mtx2, disparity = spatial.procrustes(template, moving_shape)
    return mtx2
def update_spokes(skeletal_pts, bdry_pts):
    k, _ = skeletal_pts.shape
    appender = vtk.vtkAppendPolyData()
    for i in range(k):
        spoke_poly = viz.form_spoke_poly(skeletal_pts[i, :], bdry_pts[i, :])
        appender.AddInputData(spoke_poly)
    appender.Update()
    new_spokes = appender.GetOutput()

    return new_spokes
def select_key_skeletal_pts(skeleton, bdry_pts, num_up_spokes, num_crest_spokes, surface_mesh, spokes_poly):
    
    spoke_1st_spine_end_up = viz.form_spoke_poly(skeleton[0, :], bdry_pts[0, :])
    spoke_1st_spine_end_down = viz.form_spoke_poly(skeleton[num_up_spokes, :], bdry_pts[num_up_spokes, :])
    spoke_2nd_spine_end_up = viz.form_spoke_poly(skeleton[num_up_spokes//2, :], bdry_pts[num_up_spokes//2, :])

    spoke_2nd_spine_end_down = viz.form_spoke_poly(skeleton[num_up_spokes//2+num_up_spokes, :], bdry_pts[num_up_spokes//2 + num_up_spokes, :])

    spoke_center_up = viz.form_spoke_poly(skeleton[num_up_spokes//4, :], bdry_pts[num_up_spokes//4, :])
    spoke_center_down = viz.form_spoke_poly(skeleton[num_up_spokes//4+num_up_spokes, :], bdry_pts[num_up_spokes//4+num_up_spokes, :])

    spoke_fold_1st = viz.form_spoke_poly(skeleton[-num_crest_spokes//4, :], bdry_pts[-num_crest_spokes//4, :])
    spoke_fold_2nd = viz.form_spoke_poly(skeleton[-3*num_crest_spokes//4, :], bdry_pts[-3*num_crest_spokes//4, :])
    selected_appender = vtk.vtkAppendPolyData()
    selected_appender.AddInputData(spoke_1st_spine_end_up)
    selected_appender.AddInputData(spoke_1st_spine_end_down)
    selected_appender.AddInputData(spoke_2nd_spine_end_up)
    selected_appender.AddInputData(spoke_2nd_spine_end_down)
    selected_appender.AddInputData(spoke_center_up)
    selected_appender.AddInputData(spoke_center_down)
    selected_appender.AddInputData(spoke_fold_1st)
    selected_appender.AddInputData(spoke_fold_2nd)
    selected_appender.Update()
#    viz.overlay_polydata(spokes_poly, surface_mesh, selected_appender.GetOutput())

    np_selection = np.concatenate((skeleton[0, :], bdry_pts[0, :], skeleton[num_up_spokes, :], bdry_pts[num_up_spokes, :], \
                                   skeleton[num_up_spokes//2, :], bdry_pts[num_up_spokes//2, :], skeleton[num_up_spokes//2+num_up_spokes, :], bdry_pts[num_up_spokes//2 + num_up_spokes, :],\
                                   skeleton[num_up_spokes//4, :], bdry_pts[num_up_spokes//4, :], skeleton[num_up_spokes//4+num_up_spokes, :], bdry_pts[num_up_spokes//4+num_up_spokes, :],\
                                   skeleton[-num_crest_spokes//4, :], bdry_pts[-num_crest_spokes//4, :], skeleton[-3*num_crest_spokes//4, :], bdry_pts[-3*num_crest_spokes//4, :]))

    return selected_appender.GetOutput(), np_selection
def load_group(case_ids, input_dir, template=None, vtk_dir=None, obj=None):
    all_dirs = []
    all_radii = []
    all_bdry = []
    all_skeleton = []
    surfaces = []

    # first_ellipsoid_file = '../../data/best_fitting_ellipsoid.npy'
    # v_ref = []
    # if os.path.exists(first_ellipsoid_file):
    #     with open(first_ellipsoid_file, 'rb') as f:
    #         print("load reference principal directions")
    #         v_ref = np.load(f)

    temp_bdry = template
    processed = 0
    for case_id in case_ids:
#        if case_id == "751794": continue
        up_spokes, down_spokes, crest_spokes = load_srep(input_dir + case_id)

#        surface_file = vtk_dir + '/stx_noscale_' + case_id + '_V06_t1w_RAI_Bias_label_pp_surfSPHARM.vtk'
        #surface_file = os.path.join(input_dir, case_id, obj)
        for files in glob.glob(input_dir + case_id + "/*.vtk"):
            surface_file = files
            break
        # if not os.path.exists(surface_file):
        #     surface_file = '/playpen/data/non-aligned/CaudL/stx_noscale_' + case_id + '_V06_t1w_RAI_Bias_label_pp_surfSPHARM.vtk'
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(surface_file)
        reader.Update()
        surface_mesh = reader.GetOutput()

        #bdry_pts, skeletal_pts = obtain_all_bdry_pts(up_spokes, down_spokes, crest_spokes, surface_mesh)
        skeletal_pts, bdry_pts, spokes_poly = obtain_bdry_skeletal_pts(up_spokes, down_spokes, crest_spokes, surface_mesh)
        _, np_selected_spokes = select_key_skeletal_pts(skeletal_pts, bdry_pts, up_spokes.GetNumberOfPoints(), crest_spokes.GetNumberOfPoints(), surface_mesh, spokes_poly)
        #viz.overlay_polydata(spokes_poly, surface_mesh)
        if temp_bdry is None:
            temp_bdry = bdry_pts
        d, rotated_bdry_pts, tform = procrustes.procrustes(temp_bdry, bdry_pts)
        dist_before = procrustes.procrustes_distance(temp_bdry, bdry_pts)
        # if dist_before > 100:
        #     print(case_id)
        #     viz.compare_first_spokes(spokes_poly, spokes_poly, surface_mesh)
        #invert_rotation_bdry = procrustes.transform_points(tform, rotated_bdry_pts)
        #dist_invert = procrustes.procrustes_distance(temp_bdry, invert_rotation_bdry)
        #print(dist_before, dist_invert)
#        viz.scatter_pts(temp_bdry, rotated_bdry_pts)
        #rotated_bdry_pts = procrustes_alignment(temp_bdry, bdry_pts)

        rotated_skeletal_pts = procrustes.transform_points(tform, skeletal_pts)
        new_spokes = update_spokes(rotated_skeletal_pts, rotated_bdry_pts)

        dirs, radii, bdry, skeleton = obtain_features(new_spokes)
        all_dirs.append(dirs)
        all_radii.append(radii)#down_dirs[0, :])
        all_bdry.append(bdry)
        all_skeleton.append(skeleton)
        surfaces.append(surface_mesh)
        processed += 1
        # from matplotlib import pyplot as plt
        # fig = plt.figure()
        # ax = fig.add_subplot(projection='3d')
        # ax.scatter(skeleton[:, 0], skeleton[:, 1], skeleton[:, 2])

        # plt.show()

        print('Finished: ' + str(processed))

    return {"dirs": all_dirs, "radii": all_radii, "bdry_pts": all_bdry, "skeletal_pts": all_skeleton, "template":temp_bdry, "ids": case_ids, "key_spokes": np_selected_spokes}

def load_sreps(input_dir = '/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/', data_file = '../../data/asd_hipp_data.npy', vtk_dir='/playpen/data/non-aligned/HippL/', obj_name="hipp.vtk", redo=False):

    if os.path.exists(data_file) and not redo:
        with open(data_file, 'rb') as f:
            pos_feats = np.load(f, allow_pickle=True)
            neg_feats = np.load(f, allow_pickle=True)
        return pos_feats, neg_feats
    pos_feats = load_group(pos_ids, input_dir, vtk_dir=vtk_dir, obj=obj_name)
    neg_feats = load_group(neg_ids, input_dir, pos_feats['template'], vtk_dir= vtk_dir, obj=obj_name)

    with open(data_file, 'wb') as f:
        print("writing file")
        np.save(f, pos_feats)
        np.save(f, neg_feats)

    return np.array(pos_feats), np.array(neg_feats)
if __name__ == '__main__':
    pos_feats, _ = load_sreps()
    y = np.array(['pos'] * 34 + ['neg'] * 143)
    viz.pca_scatter_plot(pos_feats['bdry_pts'], y)

# print('Done')
#viz.viz_directions_distribution(np.array(first_up_dirs), np.array(first_down_dirs))