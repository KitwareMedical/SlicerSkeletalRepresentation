### This file is a demo script for proposing a new method regarding to non-Euclidean AJIVE
### The method rotates data points to a common mean on the sphere
### Then these data points are approximated in linear space centered at the north pole
### Then apply AJIVE to decompose
### Finally rotate the decomposed component back

import vtk
import numpy as np
import os
import logging
from shanapy import models as sm
from shanapy import utils

def collect_more_dirs_around_ends(total_num_sreps=60, data_folder='../../data_diff_mean/'):
    """
    Collect more dirs from top and bot in the same order
    These spokes are around ends of objects, to avoid the bending axis

    Return top_dirs (n x k x 3)
    bot_dirs (n x k x 3)
    """
    logging.warning('The data folder is ' + data_folder)
    srep_folder = os.path.join(data_folder, 's_reps')

    top_dirs = []
    bot_dirs = []
    top_bend_dirs = []
    ## select spokes near the ends of objects
    ### 111 spokes, too many
    # range_starts = [0, 8, 20, 32, 44]
    # range_ends = [5, 17, 29, 41, 49]

    ### much less spokes in each ellipsoid
    spoke_ids = [8, 32, 44, 68, 80, 104, 116, 140]
    for i in range(total_num_sreps):
        if i == 17 or i == 47:
#        if i == 6 or i == 8:
            ## inconsistent ordering; exclude these cases
            continue
        top_srep_file_name = os.path.join(srep_folder, 'top_srep' + str(i) + '.vtk')
        bot_srep_file_name = os.path.join(srep_folder, 'bot_srep' + str(i) + '.vtk')
        top_joint_srep_file = os.path.join(srep_folder, 'top_joint_srep' + str(i) + '.vtk')

        top_reader = vtk.vtkPolyDataReader()
        top_reader.SetFileName(top_srep_file_name)
        top_reader.Update()
        top_srep = top_reader.GetOutput()

        bot_reader = vtk.vtkPolyDataReader()
        bot_reader.SetFileName(bot_srep_file_name)
        bot_reader.Update()
        bot_srep = bot_reader.GetOutput()

        gt_reader = vtk.vtkPolyDataReader()
        gt_reader.SetFileName(top_joint_srep_file)
        gt_reader.Update()
        top_joint_srep = gt_reader.GetOutput()

        append_filter = vtk.vtkAppendPolyData()
        top_dirs_this_obj = []
        bot_dirs_this_obj = []
        top_dirs_bend_only = []
        # for iter_num, start in enumerate(range_starts):
        #     id_start = start * 3
        #     id_end = range_ends[iter_num] * 3

        #     spoke_ids = np.arange(id_start, id_end)
        for spoke_id in spoke_ids:
            base_pt_id = spoke_id * 2
            bdry_pt_id = base_pt_id + 1
            base_pt = np.array(top_srep.GetPoint(base_pt_id))
            bdry_pt = np.array(top_srep.GetPoint(bdry_pt_id))
            top_spoke = sm.Spoke(base_pt=base_pt, bdry_pt=bdry_pt)

            bot_base_pt = np.array(bot_srep.GetPoint(base_pt_id))
            bot_bdry_pt = np.array(bot_srep.GetPoint(bdry_pt_id))
            bot_spoke = sm.Spoke(base_pt=bot_base_pt, bdry_pt=bot_bdry_pt)

            top_bend_base_pt = np.array(top_joint_srep.GetPoint(base_pt_id))
            top_bend_bdry_pt = np.array(top_joint_srep.GetPoint(bdry_pt_id))
            top_bend_spoke = sm.Spoke(base_pt=top_bend_base_pt, bdry_pt=top_bend_bdry_pt)

            ## convert to spherical coordinates
            # theta_top, phi_top, _ = cart2sph(top_spoke.U)
            # theta_bot, phi_bot, _ = cart2sph(bot_spoke.U)
            # theta_gt, phi_gt, _ = cart2sph(top_bend_spoke.U)
            # top_dirs_this_obj.append([theta_top, phi_top])
            # bot_dirs_this_obj.append([theta_bot, phi_bot])
            # top_dirs_bend_only.append([theta_gt, phi_gt])

            top_dirs_this_obj.append(top_spoke.U)
            bot_dirs_this_obj.append(bot_spoke.U)
            top_dirs_bend_only.append(top_bend_spoke.U)
        top_dirs_vec = np.array(top_dirs_this_obj)
        bot_dirs_vec = np.array(bot_dirs_this_obj)
        top_bend_vec = np.array(top_dirs_bend_only)
        top_dirs.append(top_dirs_vec)
        bot_dirs.append(bot_dirs_vec)
        top_bend_dirs.append(top_bend_vec)
    top_dirs = np.array(top_dirs)
    bot_dirs = np.array(bot_dirs)
    top_bend_dirs = np.array(top_bend_dirs)

    return top_dirs, bot_dirs, top_bend_dirs

def collect_dirs(num_dirs=1, data_folder = '../data_back'):
    """
    Select pairs of directions in #total_num_sreps configurations.
    For each configuration, select #num_dirs spokes' directions from
    top and bottom ellipsoid respectively.
    The collection of top directions form the first returned block,
    while the bottom forms the second block.


    The pairs consists of a top spoke and the bottom spoke found by the
    algorithm similar to compute linking structure.

    In preliminary research only one top spoke and one bot spoke are
    selected from each configuration. Thus, the first returned block
    is of shape (total_num_sreps, 3), because the direction is 3-tuple.

    """
#    srep_folder = '../data/s_reps/'
    srep_folder = os.path.join(data_folder, 's_reps')
    linker = Linker()
    selected_top_dirs = []
    linked_bot_dirs = []
    linked_spoke_ids = []
    top_bend_dirs = []
    num_changed_linked_spoke = 0
    num_north_target = 0
    total_num_sreps = 60#len([name for name in os.listdir(srep_folder) if os.path.isfile(os.path.join(srep_folder, name))]) // 2 - 1

    for i in range(total_num_sreps):
        top_srep_file_name = os.path.join(srep_folder, 'top_srep' + str(i) + '.vtk')
        bot_srep_file_name = os.path.join(srep_folder, 'bot_srep' + str(i) + '.vtk')
        top_joint_srep_file = os.path.join(srep_folder, 'top_joint_srep' + str(i) + '.vtk')

        top_reader = vtk.vtkPolyDataReader()
        top_reader.SetFileName(top_srep_file_name)
        top_reader.Update()
        top_srep = top_reader.GetOutput()

        bot_reader = vtk.vtkPolyDataReader()
        bot_reader.SetFileName(bot_srep_file_name)
        bot_reader.Update()
        bot_srep = bot_reader.GetOutput()

        top_joint_reader = vtk.vtkPolyDataReader()
        top_joint_reader.SetFileName(top_joint_srep_file)
        top_joint_reader.Update()
        top_joint_srep = top_joint_reader.GetOutput()

        num_spokes = top_srep.GetNumberOfCells()
        spoke_ids = np.random.randint(0, num_spokes, size=num_dirs)
        if num_dirs is None:
            # select all spokes
            spoke_ids = np.arange(num_spokes)

        ## temporaily select one spoke
        spoke_ids = np.array([82])
#        spoke_ids = np.array([10])
        for spoke_id in spoke_ids:
            base_pt_id = spoke_id * 2
            bdry_pt_id = base_pt_id + 1
            base_pt = np.array(top_srep.GetPoint(base_pt_id))
            bdry_pt = np.array(top_srep.GetPoint(bdry_pt_id))

            s = bdry_pt - base_pt
            r = np.linalg.norm(s)
            u = s / r
            target_spoke = Spoke(r, u, base_pt, None)
            ### due to inconsistently correspondence in orientation of surface mesh generated by SPHARM-PDM
            ### here exclude those spokes in the opposite hemisphere
            if geodesic_dist(target_spoke.U, np.array([0, 0, 1])) < np.pi / 2:
                logging.warn('The target spoke is on the north hemisphere: ' + str(i))
                num_north_target += 1
                continue

            ### ground truth of joint variation of top spokes
            top_joint_base_pt = np.array(top_joint_srep.GetPoint(base_pt_id))
            top_joint_bdry_pt = np.array(top_joint_srep.GetPoint(bdry_pt_id))
            bend_dir = (top_joint_bdry_pt - top_joint_base_pt) / np.linalg.norm(top_joint_bdry_pt - top_joint_base_pt)
            top_bend_dirs.append(bend_dir.tolist())

            ### find linked spoke that links in original ellipsoids
            # target_spoke, linked_spoke, linked_spoke_id =\
            #                                         linker.link_spoke(target_spoke, bot_srep)

            # if linked_spoke_id != 10:
            #     continue
            ### find link correspondence with same spoke id
            spoke_id = 10
            base_pt_id = spoke_id * 2
            bdry_pt_id = base_pt_id + 1

            bot_base_pt = np.array(bot_srep.GetPoint(base_pt_id))
            bot_bdry_pt = np.array(bot_srep.GetPoint(bdry_pt_id))
            bot_s = bot_bdry_pt - bot_base_pt
            bot_r = np.linalg.norm(bot_s)
            bot_u = bot_s / bot_r
            linked_spoke = Spoke(bot_r, bot_u, bot_base_pt, None)
            linked_spoke_id = spoke_id

            # if linked_spoke_id > 71:
            #     logging.warn('The linked spoke is > 71, case id ' + str(i) + ', its linked id is ' + str(linked_spoke_id))
            #     num_changed_linked_spoke += 1
            #     continue
            ## test code for visualizing variations' directions of bending and twisting

            selected_top_dirs.append(target_spoke.U.tolist())
            linked_bot_dirs.append(linked_spoke.U.tolist())
            linked_spoke_ids.append(linked_spoke_id)

            # target_spoke_poly = form_spoke_poly(base_pt, bdry_pt)
            # gt_spoke_poly = form_spoke_poly(bot_base_pt, bot_bdry_pt)
            # linked_spoke_poly = form_spoke_poly(linked_spoke.p, linked_spoke.getB())
            # append_filter = vtk.vtkAppendPolyData()
            # append_filter.AddInputData(target_spoke_poly)
            # append_filter.AddInputData(gt_spoke_poly)

            # append_filter.Update()
            # overlay_polydata(target_spoke_poly, linked_spoke_poly, gt_spoke_poly)
    np_selected_top_dirs = np.asarray(selected_top_dirs)
    np_linked_bot_dirs = np.asarray(linked_bot_dirs)
    np_linked_spoke_ids = np.asarray(linked_spoke_ids)
    np_bend_dirs = np.asarray(top_bend_dirs)

    logging.warn('In total there are ' + str(num_north_target) + ' cases have south spoke')
    return np_selected_top_dirs, np_linked_bot_dirs, np_linked_spoke_ids, np_bend_dirs


def jive_analysis(X, Y, init_rank_x, init_rank_y):
    """
    Decompose the joint and individual variations of X and Y
    These two blocks are supposed to be zero centered.

    Return the join and individual variations
    """
    from jive.AJIVE import AJIVE
    ajive = AJIVE(init_signal_ranks={'x': init_rank_x, 'y': init_rank_y})
    ajive.fit(blocks={'x': X, 'y': Y})
    return ajive
## previous name: test_ajive_spherical_coords
def get_ajive_joint(top_block, bot_block, rank_top=2, rank_bot=2):
    """
    Input two blocks of dim n x d

    Return reconstructed joint components
    """
    ajive_result = jive_analysis(top_block, bot_block, rank_top, rank_bot)
    x_joint = ajive_result.get_full_block_estimates()['x']['joint'].to_numpy()

    y_joint = ajive_result.get_full_block_estimates()['y']['joint'].to_numpy()

    x_recon = x_joint#.reshape((n, -1, 2))
    y_recon = y_joint#.reshape((n, -1, 2))

    return x_recon, y_recon

def euclidean_joint_analysis():
    top, bot, top_bend = collect_more_dirs_around_ends()
    utils.viz_directions_distribution(top.reshape(-1, 3), bot.reshape(-1, 3), pts_joint=top_bend.reshape(-1, 3))
    n, k, _ = top.shape
    mu_top = np.mean(top.reshape(n, -1), axis=0)[np.newaxis, :]
    mu_bot = np.mean(bot.reshape(n, -1), axis=0)[np.newaxis, :]

    centered_top = top.reshape(n, -1) - mu_top
    centered_bot = bot.reshape(n, -1) - mu_bot
    #grid_search_ranks(24, centered_top, centered_bot, top_bend=top_bend)
    top_joint_component, bot_joint_component = get_ajive_joint(centered_top, centered_bot, 2, 13)
    top_joint_component = top_joint_component + mu_top
    bot_joint_component = bot_joint_component + mu_bot
    ## compute the mean square error
    print("============== Results from AJIVE Euclidean coords ===================")
    top_mse = ((top_joint_component - top_bend.reshape(n, -1)) ** 2).mean()
    print("Mean square error of top joint component is %f" % top_mse)

    bot_mse = ((bot_joint_component - bot.reshape(n, -1)) ** 2).mean()
    print("Mean square error of bot joint component is %f" % bot_mse)
    
    top_dist_from_unit_vectors = utils.avg_distance_from_sphere(top_joint_component.reshape(-1, 3))
    bot_dist_from_unit_vectors = utils.avg_distance_from_sphere(bot_joint_component.reshape(-1, 3))
    print("Average distance from top unit sphere: %f" % top_dist_from_unit_vectors)
    print("Average distance from bot unit sphere: %f" % bot_dist_from_unit_vectors)
    print("======================================================================")

    ## visualize the result joint component
    utils.viz_joint_variation_on_sphere(top.reshape(-1, 3), bot.reshape(-1, 3),
                                  top_joint_component.reshape(-1, 3),
                                  bot_joint_component.reshape(-1, 3),title="Euclidean AJIVE",
                                  x_joint_gt=top_bend.reshape(-1, 3))
    # viz_directions_distribution(top_joint_component.reshape(-1, 3),
    #                             bot_joint_component.reshape(-1, 3), title="Euclidean AJIVE")
def grid_search_ranks(max_rank, centered_top, centered_bot, sph_top=None, top_bend=None):
    """
    max_rank: integer
    centered_top / centered_bot: two centered blocks (X, Y) of shape n x k x d, where k is #spokes
    sph_top: if centered_bot/top is represented by spherical coords, to reconstruct the
               cartesian coords for top and bot, need sph_top of shape n x k x 3, the (az, elev,r)
    top_bend: ground truth represented in cartesian coords
    """
    ## grid search the best ranks
    min_diff = []
    indices = []

    n = centered_top.shape[0]
    for i in range(2, max_rank):
        for j in range(2, max_rank):
            top_joint_component, _ = get_ajive_joint(centered_top, centered_bot, i, j)
            ## map back to sphere
            if sph_top is not None:
                top_joint_component = top_joint_component.reshape(n, -1, 2)
                top_joint_component = np.concatenate((top_joint_component, sph_top[:, :, 2][:, :, np.newaxis]), axis=2)
                top_joint_component = utils.sph2cart(top_joint_component)
            else:
                ## input cartesian coords
                top_joint_component = top_joint_component.reshape(n, -1, 3)
            ## compute the distance from the ground truth
            min_diff.append(((top_joint_component - top_bend) ** 2).mean())
            indices.append((i, j))
    # plt.plot(min_diff)
    # plt.show()
    # print("Minimal difference from ground truth is %f" % np.nanmin(min_diff))
    # print("Corresponding ranks are {}".format(indices[np.nanargmin(min_diff)]))
def uncentered_sph_coords_analysis():
    """
    Convert cartesian coordinates to spherical coordinates and perform joint analysis with AJIVE
    """
    top, bot, top_bend = collect_more_dirs_around_ends()
    n, k, _ = top.shape
    sph_top = utils.cart2sph(top)
    sph_bot = utils.cart2sph(bot)
    log_top = sph_top[:, :, :2]
    log_bot = sph_bot[:, :, :2]

    ## use a grid search of best ranks (3, 13)
    top_joint_component, bot_joint_component = get_ajive_joint(log_top.reshape(n, -1), log_bot.reshape(n, -1), 3, 13)
    ## map back to cartesian coords
    top_joint_component = top_joint_component.reshape(n, -1, 2)
    top_joint_component = np.concatenate((top_joint_component, sph_top[:, :, 2][:, :, np.newaxis]), axis=2)
    top_joint_component = utils.sph2cart(top_joint_component)

    ## compute the mean square error
    print("============== Results from AJIVE spherical coords ===================")
    top_mse = ((top_joint_component - top_bend) ** 2).mean()
    print("Mean square error of top joint component is %f" % top_mse)

    bot_joint_component_sph = np.concatenate((bot_joint_component.reshape(n, -1, 2), sph_bot[:, :, 2][:, :, np.newaxis]), axis=2)
    bot_joint_cart = utils.sph2cart(bot_joint_component_sph)
    bot_mse = ((bot_joint_cart - bot) ** 2).mean()
    print("Mean square error of bot joint component is %f" % bot_mse)
    print("======================================================================")
def centered_sph_coords_analysis():
    top, bot, top_bend = collect_more_dirs_around_ends()
    n, k, _ = top.shape
    sph_top = utils.cart2sph(top)
    sph_bot = utils.cart2sph(bot)
    az_elev_top = sph_top[:, :, :2].reshape(n, -1)
    az_elev_bot = sph_bot[:, :, :2].reshape(n, -1)

    mu_top = np.mean(az_elev_top, axis=0)[np.newaxis, :]
    mu_bot = np.mean(az_elev_bot, axis=0)[np.newaxis, :]
    ## mu_theta = mu_phi = 0
    centered_top = az_elev_top - mu_top
    centered_bot = az_elev_bot - mu_bot

    top_joint_component, bot_joint_component = get_ajive_joint(centered_top.reshape(n, -1), centered_bot.reshape(n, -1), 2, 13)# 6, 2)
    ## map back to sphere
    # exp_top_joint = exp_maps(top_joint_component)
    top_joint_component = (top_joint_component + mu_top).reshape(n, -1, 2)
    top_joint_component = np.concatenate((top_joint_component, sph_top[:, :, 2][:, :, np.newaxis]), axis=2)
    top_joint_component = utils.sph2cart(top_joint_component)
    ## compute the mean square error
    print("============== Results from AJIVE aligned spherical coords ===================")
    top_mse = ((top_joint_component - top_bend) ** 2).mean()
    print("Mean square error of top joint component is %f" % top_mse)

    bot_joint_component = bot_joint_component + mu_bot
    bot_joint_component_sph = np.concatenate((bot_joint_component.reshape(n, -1, 2), sph_bot[:, :, 2][:, :, np.newaxis]), axis=2)
    bot_joint_cart = utils.sph2cart(bot_joint_component_sph)
    bot_mse = ((bot_joint_cart - bot) ** 2).mean()
    print("Mean square error of bot joint component is %f" % bot_mse)

    top_dist_from_unit_vectors = utils.avg_distance_from_sphere(top_joint_component.reshape(-1, 3))
    bot_dist_from_unit_vectors = utils.avg_distance_from_sphere(bot_joint_cart.reshape(-1, 3))
    print("Average distance from top unit sphere: %f" % top_dist_from_unit_vectors)
    print("Average distance from bot unit sphere: %f" % bot_dist_from_unit_vectors)

    print("======================================================================")

    ## visualize the resulting joint component
    utils.viz_joint_variation_on_sphere(top.reshape(-1, 3), bot.reshape(-1, 3),
                              top_joint_component.reshape(-1, 3),
                                  bot_joint_cart.reshape(-1, 3),title="Non-Euclidean AJIVE",
                                  x_joint_gt=top_bend.reshape(-1, 3))
    # viz_directions_distribution(top_joint_component.reshape(-1, 3),
    #                             bot_joint_cart.reshape(-1, 3), title="Non-Euclidean AJIVE")
