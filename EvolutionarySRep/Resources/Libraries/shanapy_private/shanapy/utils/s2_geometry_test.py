from multi_view_geodesics import *
from viz import *
from s2_geometry import *
import numpy as np
from multi_object_shape_analysis import *

def test_cart2sph():
    data = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    data_sph = cart2sph(data)
    print(data_sph)
    data_cart = sph2cart(data_sph)
    print(data_cart)

def test_geod_center(x):
    return geodesic_center(x)
def test_geod_dist():
    ## single directional data
    # top_features, bot_features, linked_ids, bend_dirs = collect_dirs(data_folder='../data')
    ## fake directional data
    top_features = np.array([[1, 0, 0], [1, 0, 0]])
    bot_features = np.array([[0, 1, 0], [0, 1, 0]])

    ## multi-spokes data
    # top_features, bot_features, top_bend = collect_more_dirs_around_ends()
    #    top_bend_mu = geodesic_center(top_bend.reshape(n, -1))
    n = top_features.shape[0]
    top_mu = geodesic_center(top_features.reshape(n, -1))
    bot_mu = geodesic_center(bot_features.reshape(n, -1))


    dist = geodesic_dist(top_mu, bot_mu)
    print(dist * 180/ np.pi)
    # dist2 = geodesic_dist(top_mu, top_bend_mu)
    # print(dist2 * 180 / np.pi)
    # dist3 = geodesic_dist(bot_mu, top_bend_mu)
    # print(dist3 * 180 / np.pi)
def test_exp_maps():
    pt = np.array([0, 1, 0])[np.newaxis, :]
    tan_vect = log_maps(x=pt)
    print(tan_vect)
    print(exp_maps(x=tan_vect))
test_exp_maps()
# test_geod_dist()

# test_cart2sph()
print('done')