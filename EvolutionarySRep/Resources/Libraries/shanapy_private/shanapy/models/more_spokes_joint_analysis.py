#from multi_object_shape_analysis import *
#from multi_view_geodesics import *
#from s2_geometry import *
from shanapy import models as sm
from shanapy import utils
from matplotlib import pyplot as plt
import numpy as np
def avg_distance_from_sphere(x):
    """
    Input x: n x 3, where n is number of points on 2-sphere
    Return the average distance of data points to the sphere
    """
    lengths = np.sqrt(np.sum(x ** 2, axis=1))
    ones = np.ones_like(lengths)
    l1_dist = np.average(np.abs(ones - lengths))
    return l1_dist
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
            top_joint_component, _ = test_ajive_spherical_coords(centered_top, centered_bot, i, j)
            ## map back to sphere
            if sph_top is not None:
                top_joint_component = top_joint_component.reshape(n, -1, 2)
                top_joint_component = np.concatenate((top_joint_component, sph_top[:, :, 2][:, :, np.newaxis]), axis=2)
                top_joint_component = sph2cart(top_joint_component)
            else:
                ## input cartesian coords
                top_joint_component = top_joint_component.reshape(n, -1, 3)
            ## compute the distance from the ground truth
            min_diff.append(((top_joint_component - top_bend) ** 2).mean())
            indices.append((i, j))
    plt.plot(min_diff)
    plt.show()
    print("Minimal difference from ground truth is %f" % np.nanmin(min_diff))
    print("Corresponding ranks are {}".format(indices[np.nanargmin(min_diff)]))
def test_ajive_spherical_coords(top_block, bot_block, rank_top=2, rank_bot=2):
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

def joint_analysis_spherical_coords():
    """
    Convert cartesian coordinates to spherical coordinates and perform joint analysis with AJIVE
    """
    top, bot, top_bend = collect_more_dirs_around_ends(74)
    n, k, _ = top.shape
    sph_top = cart2sph(top)
    sph_bot = cart2sph(bot)
    log_top = sph_top[:, :, :2]
    log_bot = sph_bot[:, :, :2]

    ## use a grid search of best ranks (3, 13)
    top_joint_component, bot_joint_component = test_ajive_spherical_coords(log_top.reshape(n, -1), log_bot.reshape(n, -1), 3, 13)
    ## map back to cartesian coords
    top_joint_component = top_joint_component.reshape(n, -1, 2)
    top_joint_component = np.concatenate((top_joint_component, sph_top[:, :, 2][:, :, np.newaxis]), axis=2)
    top_joint_component = sph2cart(top_joint_component)

    ## compute the mean square error
    print("============== Results from AJIVE spherical coords ===================")
    top_mse = ((top_joint_component - top_bend) ** 2).mean()
    print("Mean square error of top joint component is %f" % top_mse)

    bot_joint_component_sph = np.concatenate((bot_joint_component.reshape(n, -1, 2), sph_bot[:, :, 2][:, :, np.newaxis]), axis=2)
    bot_joint_cart = sph2cart(bot_joint_component_sph)
    bot_mse = ((bot_joint_cart - bot) ** 2).mean()
    print("Mean square error of bot joint component is %f" % bot_mse)
    print("======================================================================")

def joint_analysis_aligned_sph_coords():
    top, bot, top_bend = collect_more_dirs_around_ends()
    n, k, _ = top.shape
    sph_top = cart2sph(top)
    sph_bot = cart2sph(bot)
    az_elev_top = sph_top[:, :, :2].reshape(n, -1)
    az_elev_bot = sph_bot[:, :, :2].reshape(n, -1)

    mu_top = np.mean(az_elev_top, axis=0)[np.newaxis, :]
    mu_bot = np.mean(az_elev_bot, axis=0)[np.newaxis, :]
    ## mu_theta = mu_phi = 0
    centered_top = az_elev_top - mu_top
    centered_bot = az_elev_bot - mu_bot

    top_joint_component, bot_joint_component = test_ajive_spherical_coords(centered_top.reshape(n, -1), centered_bot.reshape(n, -1), 2, 13)# 6, 2)
    ## map back to sphere
    # exp_top_joint = exp_maps(top_joint_component)
    top_joint_component = (top_joint_component + mu_top).reshape(n, -1, 2)
    top_joint_component = np.concatenate((top_joint_component, sph_top[:, :, 2][:, :, np.newaxis]), axis=2)
    top_joint_component = sph2cart(top_joint_component)
    ## compute the mean square error
    print("============== Results from AJIVE aligned spherical coords ===================")
    top_mse = ((top_joint_component - top_bend) ** 2).mean()
    print("Mean square error of top joint component is %f" % top_mse)

    bot_joint_component = bot_joint_component + mu_bot
    bot_joint_component_sph = np.concatenate((bot_joint_component.reshape(n, -1, 2), sph_bot[:, :, 2][:, :, np.newaxis]), axis=2)
    bot_joint_cart = sph2cart(bot_joint_component_sph)
    bot_mse = ((bot_joint_cart - bot) ** 2).mean()
    print("Mean square error of bot joint component is %f" % bot_mse)

    top_dist_from_unit_vectors = avg_distance_from_sphere(top_joint_component.reshape(-1, 3))
    bot_dist_from_unit_vectors = avg_distance_from_sphere(bot_joint_cart.reshape(-1, 3))
    print("Average distance from top unit sphere: %f" % top_dist_from_unit_vectors)
    print("Average distance from bot unit sphere: %f" % bot_dist_from_unit_vectors)

    print("======================================================================")

    ## visualize the resulting joint component
    viz_joint_variation_on_sphere(top.reshape(-1, 3), bot.reshape(-1, 3),
                              top_joint_component.reshape(-1, 3),
                                  bot_joint_cart.reshape(-1, 3),title="Non-Euclidean AJIVE",
                                  x_joint_gt=top_bend.reshape(-1, 3))
    # viz_directions_distribution(top_joint_component.reshape(-1, 3),
    #                             bot_joint_cart.reshape(-1, 3), title="Non-Euclidean AJIVE")
def non_euclidean_jive():
    """
    General non-Euclidean jive not only suitable for directional data
    """
    joint_analysis_aligned_sph_coords()

# joint_analysis_aligned_sph_coords()