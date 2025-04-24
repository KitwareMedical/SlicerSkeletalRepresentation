"""
Given an ellipsoidal or elliptic srep, convert euclidean coordinates to srep coordinates.
Useful for looking at data in srep coords.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy import stats
from scipy import sparse
import pyacvd
import pyvista as pv

def ellipse_x_y_to_theta_t(x_arr, y_arr, a, b):
    def to_optimize(arr):
        thetas = arr[::2]
        ts = arr[1::2]
        res = np.zeros(x_arr.shape[0]*2)
        res[::2] = (a + ts*b)*np.cos(thetas) - x_arr
        res[1::2] = (b + ts*a)*np.sin(thetas) - y_arr
        return res
    
    guess = np.zeros(x_arr.shape[0]*2)
    opt = scipy.optimize.root(to_optimize, guess)
    theta_arr = opt.x[::2]
    t_arr = opt.x[1::2]
    return theta_arr, t_arr

def ellipse_x_y_to_theta_t2(x_arr, y_arr, a, b, my_guess=None):
    # print("a, b: " , a, b)
    # print("-b/a:", -b/a)
    def to_optimize(arr):
        thetas = arr[::2]
        ts = arr[1::2]
        res = np.zeros(x_arr.shape[0]*2)
        res[::2] = ((a + ts*b)*np.cos(thetas) - x_arr)**2
        res[1::2] = ((b + ts*a)*np.sin(thetas) - y_arr)**2
        return np.sum(res)
    
    guess = np.zeros(x_arr.shape[0]*2)
    if my_guess is not None:
        guess = my_guess
    bounds = [(-np.pi, np.pi)] * guess.shape[0]
    bounds[1::2] = [(-b/a, None)] * x_arr.shape[0]
    opt = scipy.optimize.minimize(to_optimize, guess, method='TNC' ,bounds=bounds)
    theta_arr = opt.x[::2]
    t_arr = opt.x[1::2]
    return theta_arr, t_arr

def ellipse_x_y_to_theta_t_lineartime(x_arr, y_arr, a, b, my_guess=None, my_method=None):
    theta_arr = np.zeros_like(x_arr)
    t_arr = np.zeros_like(x_arr)

    for i in range(x_arr.shape[0]):
        print(i)
        x = x_arr[i]
        y = y_arr[i]

        def to_optimize(arr):
            theta = arr[0]
            t = arr[1]
            res = ((a + t*b)*np.cos(theta) - x)**2
            res += ((b + t*a)*np.sin(theta) - y)**2
            return res
        
        guess = np.zeros(2)
        if my_guess is not None:
            guess = my_guess
        bounds = [(-np.pi, np.pi), (-b/a + 1e-15, None)]
        # Optimization methods that work for me:
        # TNC, COBYLA, Powell
        if my_method is None:
            opt = scipy.optimize.minimize(to_optimize, guess, method='trust-constr', bounds=bounds)
        else:
            opt = scipy.optimize.minimize(to_optimize, guess, method=my_method, bounds=bounds)

        theta_arr[i] = opt.x[0]
        t_arr[i] = opt.x[1]
    
    return theta_arr, t_arr

def ellipse_theta_t_to_x_y(theta_arr, t_arr, a, b):
    """
    t goes from -b/a to 0
    """
    x_arr = (a + t_arr * b)*np.cos(theta_arr)
    y_arr = (b + t_arr * a)*np.sin(theta_arr)
    return x_arr, y_arr

def ellipse_theta_tau_to_x_y(theta_arr, tau_arr, a, b):
    t_arr = (tau_arr - 1) * b/a
    return ellipse_theta_t_to_x_y(theta_arr=theta_arr, t_arr=t_arr, a=a, b=b)

# Now working on ellipsoid
def ellipsoid_x_y_z_to_xm_ym_tau(x_arr, y_arr, z_arr, a, b, c, my_guess=None):
    # principal radii of medial axis (ellipse)
    r1 = (a**2 - c**2)/a
    r2 = (b**2 - c**2)/b
    z_arr = np.abs(z_arr) # point on medial axis is the same for + z or -z
    def to_optimize(arr):
        xms = arr[::3]
        yms = arr[1::3]
        taus = arr[2::3]
        res = np.zeros(x_arr.shape[0]*3)
        res[::3] = (xms*(a**2 + c**2 * (taus - 1))/(a**2 - c**2) - x_arr)**2
        res[1::3] = (yms*(b**2 + c**2 * (taus - 1))/(b**2 - c**2) - y_arr)**2
        res[2::3] = (taus * c * np.sqrt(1e-15 + 1 - (xms**2)/(r1**2) - (yms**2)/(r2**2)) - z_arr)**2
        return np.sum(res)
    
    def constraint_func(arr):
        xms = arr[::3]
        yms = arr[1::3]
        return (xms**2)/(r1**2) + (yms**2)/(r2**2) - 1
    
    constraint = scipy.optimize.NonlinearConstraint(fun=constraint_func, lb=-np.inf, ub=0)

    guess = np.zeros(x_arr.shape[0]*3)
    guess[2::3] = np.ones(x_arr.shape[0])
    if my_guess is not None:
        guess = my_guess

    bounds = [(None, None)] * guess.shape[0]
    bounds[2::3] = [(0, None)] * x_arr.shape[0]
    opt = scipy.optimize.minimize(to_optimize, guess, bounds=bounds, constraints=constraint)
    xms_arr = opt.x[::3]
    yms_arr = opt.x[1::3]
    taus_arr = opt.x[2::3]
    return xms_arr, yms_arr, taus_arr

def ellipsoid_x_y_z_to_xm_ym_tau2(x_arr, y_arr, z_arr, a, b, c, my_guess=None):
    r1 = (a**2 - c**2)/a
    r2 = (b**2 - c**2)/b
    z_arr = np.abs(z_arr) # point on medial axis is the same for + z or -z
    def to_optimize(arr):
        xms = arr[::3]
        yms = arr[1::3]
        taus = arr[2::3]
        res = np.zeros(x_arr.shape[0]*3)
        res[::3] = (xms*(a**2 + c**2 * (taus - 1))/(a**2 - c**2) - x_arr)**2
        res[1::3] = (yms*(b**2 + c**2 * (taus - 1))/(b**2 - c**2) - y_arr)**2
        res[2::3] = (taus * c * np.sqrt(1e-15 + 1 - (xms**2)/(r1**2) - (yms**2)/(r2**2)) - z_arr)**2
        return np.sum(res)
    
    def constraint_func(arr):
        xms = arr[::3]
        yms = arr[1::3]
        return (xms**2)/(r1**2) + (yms**2)/(r2**2) - 1
    
    n = x_arr.shape[0]
    # 0,0,1,1,2,2,3,3,4,4,5,5, ...
    jac_row_ind = np.concatenate([
        np.arange(n)[:,None],
        np.arange(n)[:,None],
    ], axis=1).reshape(-1)
    # 0, 1, 3, 4, 6, 7, 9, 10, ...
    jac_col_ind = np.concatenate([
        (np.arange(n) * 3)[:,None],
        (np.arange(n) * 3 + 1)[:,None],
    ], axis=1).reshape(-1)

    def constraint_jac(arr):
        a1 = 2 * arr[::3] / (r1**2)
        a2 = 2 * arr[1::3] / (r2**2)
        jac_data = np.concatenate([
            a1[:,None],
            a2[:,None],
        ], axis=1).reshape(-1)
        return scipy.sparse.csr_matrix((jac_data, (jac_row_ind, jac_col_ind)), shape=(n,3*n))

    constraint = scipy.optimize.NonlinearConstraint(fun=constraint_func, lb=-np.inf, ub=0, jac=constraint_jac)

    guess = np.zeros(x_arr.shape[0]*3)
    guess[2::3] = np.ones(x_arr.shape[0])
    if my_guess is not None:
        guess = my_guess

    bounds = [(None, None)] * guess.shape[0]
    bounds[2::3] = [(0, None)] * x_arr.shape[0]
    opt = scipy.optimize.minimize(to_optimize, guess, bounds=bounds, constraints=constraint)
    xms_arr = opt.x[::3]
    yms_arr = opt.x[1::3]
    taus_arr = opt.x[2::3]
    return xms_arr, yms_arr, taus_arr

def ellipsoid_x_y_z_to_xm_ym_tau3(x_arr, y_arr, z_arr, a, b, c, my_guess=None):
    r1 = (a**2 - c**2)/a
    r2 = (b**2 - c**2)/b
    z_arr = np.abs(z_arr) # point on medial axis is the same for + z or -z
    def to_optimize(arr):
        xms = arr[::3]
        yms = arr[1::3]
        taus = arr[2::3]
        res = np.zeros(x_arr.shape[0]*3)
        res[::3] = (xms*(a**2 + c**2 * (taus - 1))/(a**2 - c**2) - x_arr)**2
        res[1::3] = (yms*(b**2 + c**2 * (taus - 1))/(b**2 - c**2) - y_arr)**2
        res[2::3] = (taus * c * np.sqrt(1e-15 + 1 - (xms**2)/(r1**2) - (yms**2)/(r2**2)) - z_arr)**2
        return np.sum(res)
    
    def constraint_func(arr):
        xms = arr[::3]
        yms = arr[1::3]
        return (xms**2)/(r1**2) + (yms**2)/(r2**2) - 1
    
    n = x_arr.shape[0]
    # 0,0,1,1,2,2,3,3,4,4,5,5, ...
    jac_row_ind = np.concatenate([
        np.arange(n)[:,None],
        np.arange(n)[:,None],
    ], axis=1).reshape(-1)
    # 0, 1, 3, 4, 6, 7, 9, 10, ...
    jac_col_ind = np.concatenate([
        (np.arange(n) * 3)[:,None],
        (np.arange(n) * 3 + 1)[:,None],
    ], axis=1).reshape(-1)
    jac_data = np.ones(2*n)

    jac_sparsity = scipy.sparse.csr_matrix((jac_data, (jac_row_ind, jac_col_ind)), shape=(n,3*n))

    constraint = scipy.optimize.NonlinearConstraint(fun=constraint_func, lb=-np.inf, ub=0, finite_diff_jac_sparsity=jac_sparsity)

    guess = np.zeros(x_arr.shape[0]*3)
    guess[2::3] = np.ones(x_arr.shape[0])
    if my_guess is not None:
        guess = my_guess

    bounds = [(None, None)] * guess.shape[0]
    bounds[2::3] = [(0, None)] * x_arr.shape[0]
    opt = scipy.optimize.minimize(to_optimize, guess, bounds=bounds, constraints=constraint)
    xms_arr = opt.x[::3]
    yms_arr = opt.x[1::3]
    taus_arr = opt.x[2::3]
    return xms_arr, yms_arr, taus_arr

def ellipsoid_x_y_z_to_xm_ym_tau_lineartime(x_arr, y_arr, z_arr, a, b, c, my_guess=None):
    # principal radii of medial axis (ellipse)
    r1 = (a**2 - c**2)/a
    r2 = (b**2 - c**2)/b
    z_arr = np.abs(z_arr) # point on medial axis is the same for + z or -z
    xms_arr = np.zeros_like(x_arr)
    yms_arr = np.zeros_like(x_arr)
    taus_arr = np.zeros_like(x_arr)

    for i, xpt in enumerate(x_arr):
        x = x_arr[i]
        y = y_arr[i]
        z = z_arr[i]

        def to_optimize(arr):
            xm = arr[0]
            ym = arr[1]
            tau = arr[2]
            res = (xm*(a**2 + c**2 * (tau - 1))/(a**2 - c**2) - x)**2
            res += (ym*(b**2 + c**2 * (tau - 1))/(b**2 - c**2) - y)**2
            in_sqrt = 1 - (xm**2)/(r1**2) - (ym**2)/(r2**2)
            if (in_sqrt < 0):
                # print("sqrt val less than 0", in_sqrt)
                res += 10*np.exp(-in_sqrt)
                in_sqrt = 0
            res += (tau * c * np.sqrt(in_sqrt) - z)**2
            return res
        
        def constraint_func(arr):
            return (arr[0]**2)/(r1**2) + (arr[1]**2)/(r2**2) - 1
        
        constraint = scipy.optimize.NonlinearConstraint(fun=constraint_func, lb=-np.inf, ub=0)

        guess = np.array([0.,0.,1.])
        if my_guess is not None:
            guess = my_guess

        bounds = [(None, None), (None, None), (0, None)]
        opt = scipy.optimize.minimize(to_optimize, guess, bounds=bounds, constraints=constraint)
        xms_arr[i] = opt.x[0]
        yms_arr[i] = opt.x[1]
        taus_arr[i] = opt.x[2]

    return xms_arr, yms_arr, taus_arr

# def ellipsoid_xm_ym_tau_to_x_y_z(xm_arr, ym_arr, tau_arr, a, b, c):
#     res = np.zeros(xm_arr.shape[0]*3)
#     res[::3] = xm_arr*(a**2 + c**2 * (tau_arr  - 1))/(a**2 - c**2)
#     res[1::3] = ym_arr*(b**2 + c**2 * (tau_arr  - 1))/(b**2 - c**2)
#     r1 = (a**2 - c**2)/a
#     r2 = (b**2 - c**2)/b
#     res[2::3] = tau_arr * c * np.sqrt(1 - (xm_arr**2)/(r1**2) - (ym_arr**2)/(r2**2))
#     return (res[::3], res[1::3], res[2::3])

def ellipsoid_xm_ym_tau_to_x_y_z(xm_arr, ym_arr, tau_arr, a, b, c):
    res = np.zeros(xm_arr.shape[0]*3)
    new_xs = xm_arr*(a**2 + c**2 * (tau_arr  - 1))/(a**2 - c**2)
    new_ys = ym_arr*(b**2 + c**2 * (tau_arr  - 1))/(b**2 - c**2)
    r1 = (a**2 - c**2)/a
    r2 = (b**2 - c**2)/b
    # might have negative value in sqrt 
    val = 1 - (xm_arr**2)/(r1**2) - (ym_arr**2)/(r2**2)
    if np.any(val < 0):
        # print("negative value in sqrt: ", min(val))
        val = val.clip(0, None)
    new_zs = tau_arr * c * np.sqrt(val)
    return (new_xs, new_ys, new_zs)
