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

def parallel_transport(vec_a, vec_b):
    """

    """
    pass
def normalize_pts(x):
    """
    Map instances in x onto a unit sphere

    x: n x d matrix
    Return normalized matrix
    """
    n, d = x.shape
    ret_mat = []
    for i in range(n):
        ret_mat.append(x[i, :] / np.linalg.norm(x[i, :]))
    return np.array(ret_mat)

def rotate_to_north_pole(v, angle=None):
    """
    Rotate a unit vector v to the north pole of a unit sphere

    Return the rotation matrix
    See Rodrigues' rotation formula: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation
    Input: 1D array (i.e., a point on a unit sphere)
    """
    d = len(v) # dimension of the feature space

    ## normalize vector v
    v = v / np.linalg.norm(v)
    ## north pole coordinates d-dimension [0, ..., 0, 1]
    north_pole = [0] * (d - 1)
    north_pole.append(1)
    north_pole = np.asarray(north_pole)

    inner_prod = np.inner(north_pole, v)
    if not angle:
        angle = np.arccos(np.clip(inner_prod, -1, 1))
    if np.abs(inner_prod - 1) < 1e-15:
        return np.eye(d)
    elif np.abs(inner_prod + 1) < 1e-15:
        return -np.eye(d)
    c = v - north_pole * inner_prod
    c = c / np.linalg.norm(c)
    A = np.outer(north_pole, c) - np.outer(c, north_pole)

    rot = np.eye(d) + np.sin(angle)*A + (np.cos(angle) - 1)*(np.outer(north_pole, north_pole)\
                                                             + np.outer(c, c))
    return rot

def sph2cart(rtp):
    """
    Transform spherical to Cartesian coordinates.
    [X,Y,Z] = sph2cart(rthetaphi) transforms corresponding elements of
    data stored in spherical coordinates (azimuth TH, elevation PHI,
    radius R) to Cartesian coordinates X,Y,Z.  The arrays TH, PHI, and
    R must be the same size (or any of them can be scalar).  TH and
    PHI must be in radians.
 
    TH is the counterclockwise angle in the xy plane measured from the
    positive x axis.  PHI is the elevation angle from the xy plane.

    Input rthetaphi:  phi, theta
    Return matrix: n x 3
    """
    if len(rtp.shape) == 2:
        az, elev = rtp[:, 0], rtp[:, 1]
        r = np.ones_like(az)

        z = np.multiply(r, np.sin(elev))[:, np.newaxis]
        rcoselev = np.multiply(r, np.cos(elev))
        x = np.multiply(rcoselev, np.cos(az))[:, np.newaxis]
        y = np.multiply(rcoselev, np.sin(az))[:, np.newaxis]
        return np.hstack((x, y, z))
    else:
        ## input n x k x 2
        n = rtp.shape[0]
        ret = []
        for ni in range(n):
            feat_slice = rtp[ni, :, :]
            feat_cart = sph2cart(feat_slice)
            ret.append(feat_cart)
        return np.array(ret)
def cart2sph(xyz):
    """
    Transform Cartesian to spherical coordinates.
    [TH,PHI,R] = cart2sph(X,Y,Z) transforms corresponding elements of
    data stored in Cartesian coordinates X,Y,Z to spherical
    coordinates (azimuth TH, elevation PHI, and radius R).  The arrays
    X,Y, and Z must be the same size (or any of them can be scalar).
    TH and PHI are returned in radians.
 
    TH is the counterclockwise angle in the xy plane measured from the
    positive x axis.  PHI is the elevation angle from the xy plane.

    Input xyz: n x 3 or n x k x 3
    Return n x 3 or n x k x 3
    """
    # x, y, z = xyz
    # XsqPlusYsq = x**2 + y**2
    # r = m.sqrt(XsqPlusYsq + z**2)               # r
    # elev = m.atan2(z,m.sqrt(XsqPlusYsq))
    # az = m.atan2(y,x)

    ## vectorization to speedup
    if len(xyz.shape) == 2:
        xy = xyz[:,0]**2 + xyz[:,1]**2
        r = np.sqrt(xy + xyz[:,2]**2)  #r
        elev = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from Z-axis down
        az = np.arctan2(xyz[:,1], xyz[:,0]) # az (i.e., theta)
        return np.hstack((az[:, np.newaxis], elev[:, np.newaxis], r[:, np.newaxis]))
    else:
        ## input n x k x 3
        n = xyz.shape[0]
        ret = []
        for ni in range(n):
            feat_slice = xyz[ni, :, :]
            feat_sph = cart2sph(feat_slice)
            ret.append(feat_sph)
        return np.array(ret)
def geodesic_center(x):
    """
    Compute geodesic center of x
    input x: n x d matrix, in which each sample (d-tuple) is on the unit sphere
    """
    TOL = 1e-9
    err = 1
    mu = x[0, :]
    tau = 1 # similar to learning rate
    while err > TOL:
        projected_data = []

        ## rotate to north pole, assuming they are already on the unit sphere
        ## this is analogous to center the data
        ## while centering data doesn't make direct sense on sphere
        rot_mat = rotate_to_north_pole(mu)
        rot_x = np.matmul(x, rot_mat.T)

        ### map to tangent plane centered on the north pole
        for i in range(rot_x.shape[0]):
            log_pt = log_map(x=rot_x[i, :])
            projected_data.append(log_pt)

        ## TODO: fit a circle instead of a direct averaging
        ## mu on the tangent plane
        new_mu = np.mean(projected_data, 0)

        ## mu may be off the sphere
        new_mu_exp = exp_map(x = new_mu, t = tau)

        ## mu that is now on the unit sphere
        new_mu_exp = new_mu_exp / np.linalg.norm(new_mu_exp)

        ## rotate mu back
        inv_rot = np.linalg.inv(rot_mat)
        new_mu_exp = np.matmul(inv_rot, new_mu_exp)

        ## update iterate variables
        err = np.abs(geodesic_dist(new_mu_exp, mu))
        mu = new_mu_exp
        ## reduce the learning rate
        tau = tau / 2


    return mu
def center_at_north_pole(x):
    """
    Parallel transport the input points to the north pole
    of the unit sphere.

    Input x: n x 3 # n unit vectors
    Return centered matrix and rotation matrix
    """
    ### 1. compute intrinsic mean of data points
    mu = geodesic_center(x)

    ### 2. rotate center to north pole
    rot_mat = rotate_to_north_pole(mu)

    ### 3. rotate other data to around north pole
    return np.matmul(x, rot_mat.T), rot_mat
def geodesic_dist(r1, r2):
    """
    Geodesic distance

    Input r1, r2: n x 1 vector
    """
    k = (np.linalg.norm(r1)) ** 2 + (np.linalg.norm(r2)) ** 2
    theta = 2 * np.inner(r1, r2) / k
    if theta < -1:
        theta = -1
    elif theta > 1:
        theta = 1
    return np.abs(np.arccos(np.clip(theta, -1, 1)))

def exp_map(foot_point=None, x=None, t=1):
    """
    Exponential map
    foot_point of the map: by default it's the north pole

    See https://ronnybergmann.net/mvirt/manifolds/Sn/exp.html
    x: the vector on the tangential space TxM
    t: step
    """
    if foot_point is None:
        d = len(x)
        north_pole = [0] * (d - 1)
        north_pole.append(1)
        foot_point = np.asarray(north_pole)
    eps = 1e-13
    v = np.linalg.norm(x) + eps
    return np.cos(t * v) * foot_point + np.sin(t*v) * x / v

def log_map(center=None, x=None):
    """
    Log map on sphere
    center: center of the log map, by default it's north pole
    x: target point on the sphere to be mapped
    
    """
    assert len(x.shape) == 1
    d = len(x)
    if center is None:
        north_pole = [0] * (d - 1)
        north_pole.append(1)
        center = np.asarray(north_pole)

    u = x-np.inner(center, x)*center
    vec=u/np.linalg.norm(u)*geodesic_dist(center,x)
    return vec

def log_maps(x=None, center=None):
    """
    Log map n x d matrix from a unit sphere to a linear space centered at center

    x: n x d
    center: d x 1
    """
    assert(len(x.shape) > 1)
    n, d = x.shape
    if center is None:
        north_pole = [0] * (d - 1)
        north_pole.append(1)
        center = np.asarray(north_pole)

    ret = []
    for i in range(n):
        xi = x[i, :]
        u = xi-np.inner(center, xi)*center
        vec=u/np.linalg.norm(u)*geodesic_dist(center,xi)
        ret.append(vec)
    return np.array(ret)
def exp_maps(x=None, center=None):
    """
    Exponential map a matrix-valued data

    x: n x d, tangent vectors in the tangent space TxM
    center: d x 1
    """
    n, d = x.shape
    if center is None:
        north_pole = [0] * (d - 1)
        north_pole.append(1)
        center = np.asarray(north_pole)
    ret = []
    for i in range(n):
        xi = x[i, :]

        length_v = np.linalg.norm(xi)
        exp_x_v = np.cos(length_v) * center + np.sin(length_v) * xi / length_v
        ret.append(exp_x_v)

    return np.array(ret)
def exp_north_pole(x):
    """
    EXPNP Riemannian exponential map at North pole of S^k
    returns (k+1) x n matrix where each column is a point on a
    sphere and the input v is k x n matrix where each column
    is a point on tangent  space at north pole.

    Input: d x n matrix
    """
    d, n = x.shape
    nv = np.sqrt(np.sum(x ** 2, axis=0))
    tmp = np.sin(nv) / (nv + 1e-15)
    exp_px = np.vstack((tmp * x, np.cos(nv)))
    exp_px[:, nv < 1e-16] = np.repeat(np.vstack((np.zeros((d, 1)), 1)), np.sum(nv<1e-16), axis=1)
    return exp_px

def log_north_pole(x):
    """
    LOGNP Riemannian log map at North pole of S^k
        LogNP(x) returns k x n matrix where each column is a point on tangent
        space at north pole and the input x is (k+1) x n matrix where each column
        is a point on a sphere.
    Input: d x n matrix w.r.t. the extrinsic coords system
    Output: (d-1) x n matrix w.r.t. the coords system (tangent space) origined at the NP
    """
    d, n = x.shape
    scale = np.arccos(np.clip(x[-1, :], -1, 1)) / np.sqrt(1-x[-1, :]**2)
    scale[np.isnan(scale)] = 1
    log_px = scale * x[:-1, :]
    return log_px

#print(log_north_pole(np.array([[0., 0., 1.0]]).T))