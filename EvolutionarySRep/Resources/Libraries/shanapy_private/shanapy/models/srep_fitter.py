"""
Zhiyuans original code for fitting sreps to ellipsoids.
"""
from cmath import tau
from re import I
import vtk
import numpy as np
import logging
import pyvista as pv
from scipy.interpolate import interp1d
try:
    from .ellipsoid_srep_coords import ellipse_theta_tau_to_x_y, ellipsoid_xm_ym_tau_to_x_y_z
except ImportError:
    from ellipsoid_srep_coords import ellipse_theta_tau_to_x_y, ellipsoid_xm_ym_tau_to_x_y_z

eps = np.finfo(float).eps

def fit_ellipsoid_nick(radii, num_spine_points = 10, spine_distribution: callable=None, spine_distribution_is_from_kde: bool=False, num_steps = 4, invcdf: callable=None):
    """
    Derive s-rep for analytic ellipsoid
    with a desired distribution of spine sampling

    Assume you are fitting ellipsoid with principal radii in the e1, e2, e3 direction
    Assume raddi is a 3-tuple (a,b,c) with a > b > c

    num_steps (minimum 2) refers to the number of sample points between the spine and the fold curve of the skeleton.
    If num_steps is 2, spokes will only come from the spine or fold curve

    spine_distribution is a callable

    from_kde means that the distribution has domain (-inf, inf) and area under the curve is 1. 
    Therefore it needs boundary reflection to make it have domain (-1, 1)
    
    The srep representation we have here contains lines from the skeleton to the boundary
    """

    vertex_indices = 0, 0

    assert(len(radii) == 3)
    a, b, c = radii
    assert(a > b)
    assert(b > c)

    # r1, r2 will be the principal radii of the medial axis (ellipse)
    # r1 in the same direction as a, r2 in the same direction as b

    r1 = (a**2 - c**2) / a
    r2 = (b**2 - c**2) / b

    spine_points = None
    # if invcdf is not None:

    if spine_distribution is None:
        spine_points = np.linspace(-1., 1., num_spine_points)
    else:
        # create sampling based on distribution
        pos = np.linspace(-0.99999, 0.99999, 1000)
        cdf_pos = None
        if spine_distribution_is_from_kde:
            # border reflection for kde
            spine_distribution2 = lambda x: (spine_distribution(x) + spine_distribution(2 - x) + spine_distribution(-2 - x)) * (x <= 1) * (x >= -1)
            cdf_pos = np.cumsum(spine_distribution2(pos)) / pos.shape[0] * 2
            cdf_pos = cdf_pos / max(cdf_pos)
        else:
            cdf_pos = np.cumsum(spine_distribution(pos)) / pos.shape[0] * 2
            print(max(cdf_pos))
            cdf_pos = cdf_pos / max(cdf_pos)
        f = interp1d(cdf_pos, pos)
        # print(min(cdf_pos), max(cdf_pos))
        x_new = np.linspace(min(cdf_pos), max(cdf_pos), num_spine_points)
        # x_new = np.linspace(-1, 1, num_spine_points)
        spine_points = f(x_new)
    
    skeleton_thetas = np.arccos(spine_points) # all in interval [0, pi]
    skeleton_thetas = np.concatenate([skeleton_thetas, -skeleton_thetas[1:-1]]) # get the "underside" of the ellipse as well
    
    fold_x_arr, fold_y_arr = ellipse_theta_tau_to_x_y(skeleton_thetas, np.ones_like(skeleton_thetas), a=r1, b=r2)
    # vertices are at fold_{x , y}_arr[{0 , num_steps - 1}]

    # need to tile spine_thetas according to num_steps
    # spine and fold are handled separately

    tau1_arr = np.concatenate([
        np.ones_like(skeleton_thetas) * t for t in np.linspace(0,1,num_steps)[1:-1]
    ])

    skeleton_thetas = np.tile(skeleton_thetas, num_steps - 2)
    
    # excludes skeleton points on spine and fold
    skeleton_x_arr, skeleton_y_arr = ellipse_theta_tau_to_x_y(theta_arr=skeleton_thetas, tau_arr=tau1_arr, a=r1, b=r2)

    spine_x_arr = spine_points * (r1**2 - r2**2) / r1
    spine_y_arr = np.zeros_like(spine_x_arr)

    nonfold_skeleton_x_arr = np.concatenate([spine_x_arr, skeleton_x_arr])
    nonfold_skeleton_y_arr = np.concatenate([spine_y_arr, skeleton_y_arr])

    nonfold_tau2_arr = np.ones_like(nonfold_skeleton_x_arr) 
    fold_tau2_arr = np.ones_like(fold_x_arr) # tau2 = 1 because we want boundary point

    # computing points on the skelton now

    noncrest_x_arr, noncrest_y_arr, noncrest_z_arr \
        = ellipsoid_xm_ym_tau_to_x_y_z(nonfold_skeleton_x_arr, nonfold_skeleton_y_arr, nonfold_tau2_arr, a, b, c)
    
    # include the "underside" of the ellipsoid too
    # right now all the noncrest_z_arr are positive-- need to take their negations for the underside

    nonfold_skeleton_x_arr = np.tile(nonfold_skeleton_x_arr, 2)
    nonfold_skeleton_y_arr = np.tile(nonfold_skeleton_y_arr, 2)
    noncrest_x_arr = np.tile(noncrest_x_arr, 2)
    noncrest_y_arr = np.tile(noncrest_y_arr, 2)
    noncrest_z_arr = np.concatenate([noncrest_z_arr, -noncrest_z_arr])

    crest_x_arr, crest_y_arr, crest_z_arr \
        = ellipsoid_xm_ym_tau_to_x_y_z(fold_x_arr, fold_y_arr, fold_tau2_arr, a, b, c)
    # vertices are at crest_{x | y}_arr[{0 | num_steps - 1}]
    
    all_medial_xyz = np.concatenate([ 
        np.concatenate([ # non fold
            nonfold_skeleton_x_arr[:, np.newaxis],
            nonfold_skeleton_y_arr[:, np.newaxis],
            np.zeros_like(nonfold_skeleton_y_arr[:, np.newaxis]),
        ], axis=1),
        np.concatenate([ # fold
            fold_x_arr[:, np.newaxis],
            fold_y_arr[:, np.newaxis],
            np.zeros_like(fold_y_arr[:, np.newaxis]),
        ], axis=1)
    ], axis=0)

    all_boundary_xyz = np.concatenate([ 
        np.concatenate([ # non crest
            noncrest_x_arr[:, np.newaxis],
            noncrest_y_arr[:, np.newaxis],
            noncrest_z_arr[:, np.newaxis],
        ], axis=1),
        np.concatenate([ # crest
            crest_x_arr[:, np.newaxis],
            crest_y_arr[:, np.newaxis],
            crest_z_arr[:, np.newaxis],
        ], axis=1)
    ], axis=0)

    vertex_indices = nonfold_skeleton_x_arr.shape[0], nonfold_skeleton_x_arr.shape[0] + num_spine_points - 1
    coc_indices = [2 * vertex_indices[0], 2 * vertex_indices[1]]
    coc_indices = np.array(coc_indices)
    vertex_indices = [2 * vertex_indices[0] + 1, 2 * vertex_indices[1] + 1]
    vertex_indices = np.array(vertex_indices)

    # take care of rotations
    ### Convert spokes to visualizable elements

    spokes_poly = vtk.vtkPolyData()
    spokes_pts = vtk.vtkPoints()
    spokes_cells = vtk.vtkCellArray()

    for i in range(all_medial_xyz.shape[0]):
        id_s = spokes_pts.InsertNextPoint(all_medial_xyz[i, :])
        id_b = spokes_pts.InsertNextPoint(all_boundary_xyz[i, :])

        spoke = vtk.vtkLine()
        spoke.GetPointIds().SetId(0, id_s)
        spoke.GetPointIds().SetId(1, id_b)
        spokes_cells.InsertNextCell(spoke)

    spokes_poly.SetPoints(spokes_pts)
    spokes_poly.SetLines(spokes_cells)

    return spokes_poly, coc_indices, vertex_indices

def fit_ellipsoid(standard_ellipsoid, voxel_size=0.1, num_crest_points = 24, predefined_eig_vec=None, ell_axis=None):
    """
    Input a standard ellipsoid surface mesh
    Derive the s-rep for the ellipsoid
    """
    num_pts = standard_ellipsoid.GetNumberOfPoints()

    coords_mat = np.zeros((num_pts, 3))
    # import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d import Axes3D
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    for i in range(num_pts):
        pt = standard_ellipsoid.GetPoint(i)
        coords_mat[i, :] = pt
        # if i % 5 == 0:
        #     ax.scatter(pt[0], pt[1], pt[2])
    input_center = np.mean(coords_mat, axis=0)[np.newaxis, :]
    centered_coords = coords_mat - input_center

    # covariance = np.cov(centered_coords.T) / num_pts
    # _, s, vh = np.linalg.svd(covariance)
    # rx, ry, rz = 2 * np.sqrt(s * num_pts).T
#    plt.show()

    s, vh = np.linalg.eig(np.matmul(centered_coords.T, centered_coords))
    mass_filter = vtk.vtkMassProperties()
    mass_filter.SetInputData(standard_ellipsoid)
    mass_filter.Update()
    volume = mass_filter.GetVolume()

    idx = s.argsort()[::-1]
    s = s[idx]
    vh = vh[:, idx]
    if predefined_eig_vec is not None:
        # print("Replace principal directions")
        # print(vh)
        # print(predefined_eig_vec)
        vh = predefined_eig_vec
    rx, ry, rz = np.sqrt(s)

    ellipsoid_volume = 4 / 3.0 * np.pi * rx * ry * rz
    volume_factor = pow(volume/ ellipsoid_volume, 1.0 / 3.0)

#    volume_factor = 0.8
    ### notations consistent with wenqi's presentation
    rx *= volume_factor
    ry *= volume_factor
    rz *= volume_factor

    # if not ell_axis:
    #     rx, ry, rz = ell_axis
    mrx_o = (rx*rx-rz*rz)/rx
    mry_o = (ry*ry-rz*rz)/ry

    # ELLIPSE_SCALE = 0.9
    ELLIPSE_SCALE = 1.0
    mrb = mry_o * ELLIPSE_SCALE
    mra = mrx_o * ELLIPSE_SCALE

    delta_theta = 2 * np.pi / num_crest_points
    num_steps = 3
    skeletal_pts_x = np.zeros((num_crest_points, num_steps))
    skeletal_pts_y = np.zeros((num_crest_points, num_steps))
    skeletal_pts_z = np.zeros((num_crest_points, num_steps))

    bdry_up_x = np.zeros((num_crest_points, num_steps))
    bdry_up_y = np.zeros((num_crest_points, num_steps))
    bdry_up_z = np.zeros((num_crest_points, num_steps))

    bdry_down_x = np.zeros((num_crest_points, num_steps))
    bdry_down_y = np.zeros((num_crest_points, num_steps))
    bdry_down_z = np.zeros((num_crest_points, num_steps))

    crest_bdry_pts = np.zeros((num_crest_points, 3))
    crest_skeletal_pts = np.zeros((num_crest_points, 3))
    for i in range(num_crest_points):
        theta = np.pi - delta_theta * i
        x = mra * np.cos(theta)
        y = mrb * np.sin(theta)

        mx_ = (mra * mra - mrb * mrb) * np.cos(theta) / mra
        my_ = .0
        dx_ = x - mx_
        dy_ = y - my_

        step_size = 1.0 / float(num_steps-1)

        for j in range(num_steps):
            sp_x = mx_ + step_size * j * dx_
            sp_y = my_ + step_size * j * dy_

            # if i == 0 and j == 0:
            #     sp_x = mx_ + step_size * j * dx_ - 1
            # elif i == 0 and j == 1:
            #     sp_x = mx_ + step_size * j * dx_-0.75
            # elif i == num_crest_points / 2 and j == 1:
            #     sp_x = mx_ + step_size * j * dx_+0.75
            # elif i == num_crest_points / 2 and j == 0:
            #     sp_x = mx_ + step_size * j * dx_+1
            skeletal_pts_x[i, j] = sp_x
            skeletal_pts_y[i, j] = sp_y
            sin_spoke_angle = sp_y * mrx_o
            cos_spoke_angle = sp_x * mry_o

            # normalize to [-1, 1]
            l = np.sqrt(sin_spoke_angle ** 2 + cos_spoke_angle ** 2)
            if l >  eps:
                sin_spoke_angle /= l
                cos_spoke_angle /= l
            cos_phi = l / (mrx_o * mry_o)
            sin_phi = np.sqrt(1 - cos_phi ** 2)
            bdry_x = rx * cos_phi * cos_spoke_angle
            bdry_y = ry * cos_phi * sin_spoke_angle
            bdry_z = rz * sin_phi
            bdry_up_x[i, j] = bdry_x
            bdry_up_y[i, j] = bdry_y
            bdry_up_z[i, j] = bdry_z

            bdry_down_x[i, j] = bdry_x
            bdry_down_y[i, j] = bdry_y
            bdry_down_z[i, j] = -bdry_z

            ## if at the boundary of the ellipse, add crest spokes
            if j == num_steps - 1:
                cx = rx * cos_spoke_angle - sp_x
                cy = ry * sin_spoke_angle - sp_y
                cz = 0
                vec_c = np.asarray([cx, cy, cz])
                norm_c = np.linalg.norm(vec_c)
                dir_c = np.asarray([bdry_x - sp_x, bdry_y - sp_y, 0.0])
                dir_c = dir_c / np.linalg.norm(dir_c)

                crest_spoke = norm_c * dir_c
                crest_bdry_x = crest_spoke[0] + sp_x
                crest_bdry_y = crest_spoke[1] + sp_y
                crest_bdry_z = 0.0

                crest_bdry_pts[i] = np.asarray([crest_bdry_x, crest_bdry_y, crest_bdry_z])
                crest_skeletal_pts[i] = np.asarray([sp_x, sp_y, 0.0])
    ### Rotate skeletal/implied boundary points as boundary points of the ellipsoid
    rot_obj = np.flipud(vh.T)
    ## make this rotation matrix same with c++ computation with Eigen3
    # rot_obj[0, :] *= -1
    # rot_obj[-1, :] *= -1

    concate_skeletal_pts = np.concatenate((skeletal_pts_x.flatten()[:, np.newaxis], \
                                           skeletal_pts_y.flatten()[:, np.newaxis], \
                                           skeletal_pts_z.flatten()[:, np.newaxis]), \
                                                axis=1)
    concate_bdry_up_pts = np.concatenate((bdry_up_x.flatten()[:, np.newaxis], \
                                      bdry_up_y.flatten()[:, np.newaxis], \
                                      bdry_up_z.flatten()[:, np.newaxis]), axis=1)
    concate_bdry_down_pts = np.concatenate((bdry_down_x.flatten()[:, np.newaxis], \
                                            bdry_down_y.flatten()[:, np.newaxis], \
                                            bdry_down_z.flatten()[:, np.newaxis]), axis=1)

    second_moment_srep = np.matmul(concate_skeletal_pts.T, concate_skeletal_pts)
    s_srep, v_srep = np.linalg.eig(second_moment_srep)

    rot_srep = np.flipud(v_srep.T)

    rotation = np.flipud(np.flipud(np.matmul(rot_obj, rot_srep)).T).T

    transformed_concate_skeletal_pts = np.matmul(concate_skeletal_pts, rotation) + input_center
    transformed_concate_bdry_up_pts = np.matmul(concate_bdry_up_pts, rotation) + input_center
    transformed_concate_bdry_down_pts = np.matmul(concate_bdry_down_pts, rotation) + input_center
    transformed_crest_bdry_pts = np.matmul(crest_bdry_pts, rotation) + input_center
    transformed_crest_skeletal_pts = np.matmul(crest_skeletal_pts, rotation) + input_center

    # p = pv.Plotter()
    # points = pv.PolyData(concate_skeletal_pts)
    # trans_pts = pv.PolyData(transformed_concate_skeletal_pts)
    # p.add_mesh(points, color='red', point_size=10)
    # p.add_mesh(trans_pts, color='blue', point_size=10)
    # p.show()

    ### Convert spokes to visualizable elements
    up_spokes_poly = vtk.vtkPolyData()
    up_spokes_pts = vtk.vtkPoints()
    up_spokes_cells = vtk.vtkCellArray()
    down_spokes_poly = vtk.vtkPolyData()
    down_spokes_pts = vtk.vtkPoints()
    down_spokes_cells = vtk.vtkCellArray()
    crest_spokes_poly = vtk.vtkPolyData()
    crest_spokes_pts = vtk.vtkPoints()
    crest_spokes_cells = vtk.vtkCellArray()

    for i in range(concate_skeletal_pts.shape[0]):
        id_s = up_spokes_pts.InsertNextPoint(transformed_concate_skeletal_pts[i, :])
        id_b = up_spokes_pts.InsertNextPoint(transformed_concate_bdry_up_pts[i, :])

        id_sdwn = down_spokes_pts.InsertNextPoint(transformed_concate_skeletal_pts[i, :])
        id_down = down_spokes_pts.InsertNextPoint(transformed_concate_bdry_down_pts[i, :])

        up_spoke = vtk.vtkLine()
        up_spoke.GetPointIds().SetId(0, id_s)
        up_spoke.GetPointIds().SetId(1, id_b)
        up_spokes_cells.InsertNextCell(up_spoke)

        down_spoke = vtk.vtkLine()
        down_spoke.GetPointIds().SetId(0, id_sdwn)
        down_spoke.GetPointIds().SetId(1, id_down)
        down_spokes_cells.InsertNextCell(down_spoke)


    up_spokes_poly.SetPoints(up_spokes_pts)
    up_spokes_poly.SetLines(up_spokes_cells)
    down_spokes_poly.SetPoints(down_spokes_pts)
    down_spokes_poly.SetLines(down_spokes_cells)

    for i in range(num_crest_points):
        id_crest_s = crest_spokes_pts.InsertNextPoint(transformed_crest_skeletal_pts[i, :])
        id_crest_b = crest_spokes_pts.InsertNextPoint(transformed_crest_bdry_pts[i, :])
        crest_spoke = vtk.vtkLine()
        crest_spoke.GetPointIds().SetId(0, id_crest_s)
        crest_spoke.GetPointIds().SetId(1, id_crest_b)
        crest_spokes_cells.InsertNextCell(crest_spoke)
    crest_spokes_poly.SetPoints(crest_spokes_pts)
    crest_spokes_poly.SetLines(crest_spokes_cells)

 
    append_filter = vtk.vtkAppendPolyData()
    append_filter.AddInputData(up_spokes_poly)
    append_filter.AddInputData(down_spokes_poly)
    append_filter.AddInputData(crest_spokes_poly)
    append_filter.Update()
    return append_filter.GetOutput(), up_spokes_poly, down_spokes_poly, crest_spokes_poly

def refine_srep(srep_poly, input_mesh):
    """
    Optimize spokes in srep_poly, such that
    the tips are touching the input_mesh more closely
    """

#    num_spokes = srep_poly.GetNumberOfCells()
    num_pts = srep_poly.GetNumberOfPoints()
    num_spokes = num_pts // 2
    radii_array = np.zeros(num_spokes)
    dir_array = np.zeros((num_spokes, 3))
    base_array = np.zeros((num_spokes,3))

    ### read the parameters from s-rep
    for i in range(num_spokes):
        id_base_pt = i * 2
        id_bdry_pt = id_base_pt + 1
        base_pt = np.array(srep_poly.GetPoint(id_base_pt))
        bdry_pt = np.array(srep_poly.GetPoint(id_bdry_pt))

        radius = np.linalg.norm(bdry_pt - base_pt)
        direction = (bdry_pt - base_pt) / radius

        radii_array[i] = radius
        dir_array[i, :] = direction
        base_array[i, :] = base_pt
    def obj_func(radii, grad=None):
        """
        Square of signed distance from tips
        of spokes to the input_mesh
        """
        implicit_distance = vtk.vtkImplicitPolyDataDistance()
        implicit_distance.SetInput(input_mesh)
        total_loss = 0
        for i, radius in enumerate(radii):
            direction = dir_array[i, :]
            base_pt   = base_array[i, :]
            bdry_pt   = base_pt + radius * direction

            dist = implicit_distance.FunctionValue(bdry_pt)
            total_loss += dist ** 2
        return total_loss
    # from scipy import optimize as opt
    # minimum = opt.fmin(obj_func, radii_array)

    # minimizer = minimum[0]

    ### optimize the variables (i.e., radii)
    import nlopt
    opt = nlopt.opt(nlopt.LN_NEWUOA, len(radii_array))
    opt.set_min_objective(obj_func)
    opt.set_maxeval(2000)
    minimizer = opt.optimize(radii_array)

    ## update radii of s-rep and return the updated
    num_diff_spokes = 0
    arr_length = vtk.vtkDoubleArray()
    arr_length.SetNumberOfComponents(1)
    arr_length.SetName("spokeLength")

    arr_dirs = vtk.vtkDoubleArray()
    arr_dirs.SetNumberOfComponents(3)
    arr_dirs.SetName("spokeDirection")
    for i in range(num_spokes):
        id_base_pt = i * 2
        id_bdry_pt = id_base_pt + 1
        base_pt = base_array[i, :]
        radius = minimizer[i]
        direction = dir_array[i, :]

        new_bdry_pt = base_pt + radius * direction
        arr_length.InsertNextValue(radius)
        arr_dirs.InsertNextTuple(direction)
        srep_poly.GetPoints().SetPoint(id_bdry_pt, new_bdry_pt)

        if np.abs(np.linalg.norm(new_bdry_pt - base_pt) - radii_array[i]) > eps:
            num_diff_spokes += 1
#        srep_poly.SetPoint

    logging.info('{} spokes have been changed in the optimization.'.format(num_diff_spokes))
    srep_poly.GetPointData().AddArray(arr_length)
    srep_poly.GetPointData().AddArray(arr_dirs)
    srep_poly.Modified()
    return srep_poly
# from shanapy.utils.viz import *
# coarsen_mesh_reader = vtk.vtkSTLReader()
# coarsen_mesh_reader.SetFileName('coarsen_ellipsoid.stl')
# coarsen_mesh_reader.Update()
# mesh = coarsen_mesh_reader.GetOutput()

