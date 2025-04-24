import sys
from shanapy.models import Spoke, srep_io
from shanapy.utils import *
import os
import vtk
import numpy as np
from scipy import interpolate
from numpy import sin
import numpy.linalg as LA
#from viz_utils import *
from scipy.spatial.transform import Slerp
from scipy.spatial.transform import Rotation as R

# set the upper limit for recursion
sys.setrecursionlimit(10**6)
epsilon = 1e-7
def h1(s):
    return 2*(s * s * s) - 3*(s * s) + 1
def h2(s):
    return -2*(s * s * s) + 3*(s * s)
def h3(s):
    return (s * s * s) - 2*(s * s) + s
def h4(s):
    return (s * s * s) - (s * s)

### utilities for quaternions
def qcon(q):
    assert len(q) == 4
    # negate q_1 to q_3
    q = -q
    q[0] = -q[0]
    return q
def qlog(q):
    theta = np.arccos(q[0])
    result = theta * q / sin(theta)
    result[0] = 0
    return result
def qexp(q):
    assert len(q) == 4
    q0, q1, q2, q3 = q
    a = np.exp(q0)
    nv = np.sqrt(q1 ** 2, q2**2, q3**2)
    snv = sin(nv)

    ret = a * (snv * q/nv)
    return np.array([a * np.cos(nv), ret[1], ret[2], ret[3]])
def qmul(q0, q1):
    w0, x0, y0, z0 = q0
    w1, x1, y1, z1 = q1
    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                     x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)
class Interpolater(object):
    """
    Interpolate an srep read by the function readSrepFromXML (organized srep).
    Return in the function interpolate(...) a list of spokes.
    """
    def __init__(self):
        self.input_srep = None
        self.interpolated_dirs = None
        self.spacing = 0.0

    def _compute_derivative(self, input_srep, r, c, num_rows, num_cols):
        # use finite difference to compute derivatives

        num_crest_points = num_rows * 2 + (num_cols - 2) * 2
        num_cols = 1 + num_rows // 2
        id = r * num_cols + c
        pt0 = pt1 = [0] * 3 # pt1 - pt0 for dxdu
        pt0v = pt1v = [0] * 3 # for dxdv
        base_pts_array = input_srep.GetPointData().GetArray("basePoints")
        if r == 0:
            # first row
            pt0 = base_pts_array.GetTuple3(id)
            pt1 = base_pts_array.GetTuple3(id+num_cols)
            factor = 1.0
        elif r == num_crest_points - 1:
            # last row
            pt1 = base_pts_array.GetTuple3(id)
            pt0 = base_pts_array.GetTuple3(id-num_cols)
            factor = 1.0
        else:
            # otherwise
            pt1 = base_pts_array.GetTuple3(id + num_cols)
            pt0 = base_pts_array.GetTuple3(id - num_cols)
            factor = 0.5
        dxdu = [(pt1i - pt0i) * factor for pt1i, pt0i in zip(pt1, pt0)]

        if c == 0:
            pt1v = base_pts_array.GetTuple3(id+1)
            pt0v = base_pts_array.GetTuple3(id)
            factor = 1.0
        elif c == num_cols - 1:
            pt1v = base_pts_array.GetTuple3(id)
            pt0v = base_pts_array.GetTuple3(id-1)
            factor = 1.0
        else:
            pt1v = base_pts_array.GetTuple3(id + 1)
            pt0v = base_pts_array.GetTuple3(id - 1)
            factor = 0.5
        dxdv = [(pt1iv - pt0iv) * factor for pt1iv, pt0iv in zip(pt1v, pt0v)]

        return dxdu, dxdv
    def _interpolate_skeleton(self, relative_position, corner_pts, corner_deriv):
        u, v = relative_position
        dxdu11, dxdv11, dxdu21, dxdv21, dxdu12, dxdv12, dxdu22, dxdv22 = corner_deriv
        x11, x21, x22, x12 = corner_pts

        # poly_corners = wrap_points(corner_pts)
        # visualize(poly_corners, self.input_srep)

        hx = [[0 for i in range(4)] for j in range(4)]
        hy = [[1 for i in range(4)] for j in range(4)]
        hz = [[2 for i in range(4)] for j in range(4)]
        hx[0][0] = x11[0];          hx[0][1] = x12[0];
        hx[1][0] = x21[0];          hx[1][1] = x22[0];
        hx[2][0] = dxdu11[0];       hx[2][1] = dxdu12[0];
        hx[3][0] = dxdu21[0];       hx[3][1] = dxdu22[0];
        hx[0][2] = dxdv11[0];       hx[0][3] = dxdv12[0];
        hx[1][2] = dxdv21[0];       hx[1][3] = dxdv22[0];
        hx[2][2] = 0;               hx[2][3] = 0;
        hx[3][2] = 0;               hx[3][3] = 0;


        hy[0][0] = x11[1];          hy[0][1] = x12[1];
        hy[1][0] = x21[1];          hy[1][1] = x22[1];
        hy[2][0] = dxdu11[1];       hy[2][1] = dxdu12[1];
        hy[3][0] = dxdu21[1];       hy[3][1] = dxdu22[1];
        hy[0][2] = dxdv11[1];       hy[0][3] = dxdv12[1];
        hy[1][2] = dxdv21[1];       hy[1][3] = dxdv22[1];
        hy[2][2] = 0;               hy[2][3] = 0;
        hy[3][2] = 0;               hy[3][3] = 0;

        hz[0][0] = x11[2];       hz[0][1] = x12[2];
        hz[1][0] = x21[2];       hz[1][1] = x22[2];
        hz[2][0] = dxdu11[2];    hz[2][1] = dxdu12[2];
        hz[3][0] = dxdu21[2];    hz[3][1] = dxdu22[2];
        hz[0][2] = dxdv11[2];    hz[0][3] = dxdv12[2];
        hz[1][2] = dxdv21[2];    hz[1][3] = dxdv22[2];
        hz[2][2] = 0;            hz[2][3] = 0;
        hz[3][2] = 0;            hz[3][3] = 0;

        hu = [0] * 4
        huThx = [0] * 4
        huThy = [0] * 4
        huThz = [0] * 4
        hv = [0] * 4
        hu[0] = h1(u)
        hu[1] = h2(u)
        hu[2] = h3(u)
        hu[3] = h4(u)
        hv[0] = h1(v)
        hv[1] = h2(v)
        hv[2] = h3(v)
        hv[3] = h4(v)

        huThx[0] = hu[0] * hx[0][0] + hu[1] * hx[1][0] + hu[2] * hx[2][0] + hu[3] * hx[3][0]
        huThx[1] = hu[0] * hx[0][1] + hu[1] * hx[1][1] + hu[2] * hx[2][1] + hu[3] * hx[3][1]
        huThx[2] = hu[0] * hx[0][2] + hu[1] * hx[1][2] + hu[2] * hx[2][2] + hu[3] * hx[3][2]
        huThx[3] = hu[0] * hx[0][3] + hu[1] * hx[1][3] + hu[2] * hx[2][3] + hu[3] * hx[3][3]

        huThy[0] = hu[0] * hy[0][0] + hu[1] * hy[1][0] + hu[2] * hy[2][0] + hu[3] * hy[3][0]
        huThy[1] = hu[0] * hy[0][1] + hu[1] * hy[1][1] + hu[2] * hy[2][1] + hu[3] * hy[3][1]
        huThy[2] = hu[0] * hy[0][2] + hu[1] * hy[1][2] + hu[2] * hy[2][2] + hu[3] * hy[3][2]
        huThy[3] = hu[0] * hy[0][3] + hu[1] * hy[1][3] + hu[2] * hy[2][3] + hu[3] * hy[3][3]

        huThz[0] = hu[0] * hz[0][0] + hu[1] * hz[1][0] + hu[2] * hz[2][0] + hu[3] * hz[3][0]
        huThz[1] = hu[0] * hz[0][1] + hu[1] * hz[1][1] + hu[2] * hz[2][1] + hu[3] * hz[3][1]
        huThz[2] = hu[0] * hz[0][2] + hu[1] * hz[1][2] + hu[2] * hz[2][2] + hu[3] * hz[3][2]
        huThz[3] = hu[0] * hz[0][3] + hu[1] * hz[1][3] + hu[2] * hz[2][3] + hu[3] * hz[3][3]

        output = [0] * 3
        output[0] = huThx[0] * hv[0] + huThx[1] * hv[1] + huThx[2] * hv[2]
        output[1] = huThy[0] * hv[0] + huThy[1] * hv[1] + huThy[2] * hv[2]
        output[2] = huThz[0] * hv[0] + huThz[1] * hv[1] + huThz[2] * hv[2]

        return output
    def _scipy_slerp(self, start_dir, end_dir, num_interps):
        rotations_vect = np.vstack((np.array(start_dir)[np.newaxis, :], np.array(end_dir)[np.newaxis, :]))

        rot_rep = R.from_rotvec(rotations_vect)
        k_times = [0, 1]
        slerp = Slerp(k_times, rot_rep)

        dist = np.linspace(0, 1, num_interps)

        interp_rot = slerp(dist)

        rot_vec = interp_rot.as_rotvec()

        return rot_vec

    def _squad(self, start_dir, end_dir, dist):
        q0 = R.from_rotvec(np.array(start_dir)[np.newaxis, :]).as_quat()
        q1 = R.from_rotvec(np.array(end_dir)[np.newaxis, :]).as_quat()

        a0 = qmul(q0, qexp( -0.25 * ( qlog(qmul(qcon(q0),q1)) + qlog(qmul(qcon(q1),q0)) ) ) )

    def _slerp(self, start_dir, end_dir, dist):
        """
        Interpolate spoke directions by quaternions
        """

        u1Tu2 = start_dir[0] * end_dir[0] + start_dir[1] * end_dir[1] + start_dir[2] * end_dir[2]
        if u1Tu2 > 1.0:
            u1Tu2 = 1.0
        elif u1Tu2 < -1.0:
            u1Tu2 = -1.0

        phi = np.arccos(u1Tu2)

        output = [0] * 3
        output[0] = ( sin((1-dist)*phi)/sin(phi) )*start_dir[0] + ( sin(dist*phi)/sin(phi) )*end_dir[0]
        output[1] = ( sin((1-dist)*phi)/sin(phi) )*start_dir[1] + ( sin(dist*phi)/sin(phi) )*end_dir[1]
        output[2] = ( sin((1-dist)*phi)/sin(phi) )*start_dir[2] + ( sin(dist*phi)/sin(phi) )*end_dir[2]
        return np.array(output)

    def _compute_2nd_derivative(self, start_spoke, end_spoke=None, base_dir=None, d=None):
        """
        Approximate 2nd derivatives with finite difference
        """
        epsilon = 1e-5
        Upv1 = self._slerp(start_spoke.U, end_spoke.U, d+2*epsilon);
        Upv2 = self._slerp(start_spoke.U, end_spoke.U, d+epsilon);

        Upv4 = self._slerp(start_spoke.U, end_spoke.U, d-epsilon);
        Upv5 = self._slerp(start_spoke.U, end_spoke.U, d-2*epsilon);

        output = [0] * 3
        output[0] = 0.25 * (Upv5[0] + Upv1[0] - 2.0 * base_dir[0]);
        output[1] = 0.25 * (Upv5[1] + Upv1[1] - 2.0 * base_dir[1]);
        output[2] = 0.25 * (Upv5[2] + Upv1[2] - 2.0 * base_dir[2]);
        return output
    def _interpolate_middle(self, spoke_start, spoke_end, dist_from_start, absolute_u, absolute_v):
        """
        Interpolate the direction and radius of middle spoke between spoke_start and
        spoke_end, which distant from spoke_start dist_from_start
        Inpute absolute_uv are coords of the interpolated spoke. According to this coords, find the direction
        that is already interpolated outside and use it.
        """
        ## 1. 2nd derivatives at two ends
        assert(isinstance(spoke_start,Spoke) and isinstance(spoke_end, Spoke))
        Uvv_end = self._compute_2nd_derivative(spoke_start, spoke_end, spoke_end.U, dist_from_start)
        Uvv_start = self._compute_2nd_derivative(spoke_start, spoke_end, spoke_start.U, 0)

        ## 2. compute the middle spoke
        sum_rU = spoke_start.add(spoke_end)
        avg_rU = sum_rU / 2

        ## 3. compute the direction of the middle spoke
        half_dist = dist_from_start / 2
        # if start_dir == end_dir set middle_dir to the same
        start_dir = spoke_start.U
        end_dir = spoke_end.U

        if abs(start_dir[0] - end_dir[0]) < epsilon \
           and abs(start_dir[1] - end_dir[1]) < epsilon \
           and abs(start_dir[2] - end_dir[2]) < epsilon:
            middle_dir = start_dir
            Uvv_start = np.zeros_like(Uvv_start)
            Uvv_end = np.zeros_like(Uvv_end)
        else:
#            compare_dir = self._scipy_slerp(start_dir, end_dir, half_dist)
            middle_dir = self._slerp(start_dir, end_dir, half_dist)
#            middle_dir = self.interpolated_dirs[:, absolute_u, absolute_v]
#            middle_dir = self._squad(start_dir, end_dir, half_dist)
            # if np.linalg.norm(middle_dir - compare_dir) > 0.1:
            #     print('very different')

        ## 4. compute the radius of the middle spoke
        inner_prod1 = np.dot(middle_dir, avg_rU)
        inner_prod2 = np.dot(start_dir, Uvv_start)
        inner_prod3 = np.dot(end_dir, Uvv_end)
        middle_r = inner_prod1 - half_dist ** 2 * 0.25 * (inner_prod2 + inner_prod3)

        u_start, v_start = spoke_start.coords
        u_end, v_end = spoke_end.coords

        middle_coords = absolute_u, absolute_v

        return Spoke(middle_r, middle_dir, None, middle_coords)
    def _interpolate_quad(self, relative_position, corner_spokes, p_lambda):
        """
        Interpolate the center of a quad surrounded by 4 corner_spokes to approach the relative_position.
        The relative_position (type: float) specifies the target position want to interpolate,
        while p_lambda (2 ^ k, k <= 0) indicates how many subdivisions have been through, initialy 1
        """
        sp11, sp12, sp21, sp22 = corner_spokes
        u, v = relative_position
        absolute_u_11, absolute_v_11 = sp11.coords
        absolute_u_12, absolute_v_12 = sp12.coords
        absolute_u_21, absolute_v_21 = sp21.coords
        absolute_u_22, absolute_v_22 = sp22.coords
        ### 1. interpolate center positions on edges
        top_middle_spoke   = self._interpolate_middle(sp11, sp12, p_lambda, absolute_u_11, (absolute_v_11 + absolute_v_12) // 2)
        left_middle_spoke  = self._interpolate_middle(sp11, sp21, p_lambda, (absolute_u_11 + absolute_u_21) // 2, absolute_v_11)
        bot_middle_spoke   = self._interpolate_middle(sp21, sp22, p_lambda, absolute_u_21, (absolute_v_21 + absolute_v_22) // 2)
        right_middle_spoke = self._interpolate_middle(sp22, sp12, p_lambda, (absolute_u_12 + absolute_u_22) // 2, absolute_v_12)

        u_top_left, v_top_left = sp11.coords
        u_top_right, v_top_right = sp12.coords
        u_bot_left, v_bot_lef = sp21.coords
        u_quad_center = (u_top_left + u_top_right) // 2
        v_quad_center = (v_top_left + v_bot_lef) // 2

        ### 2. interpolate center of the quad
        vertical_center = self._interpolate_middle(top_middle_spoke, bot_middle_spoke, p_lambda, u_quad_center, v_quad_center)
        horizont_center = self._interpolate_middle(left_middle_spoke, right_middle_spoke, p_lambda, u_quad_center, v_quad_center)

        ### 3. average vertical_center and horizont_center to get a better estimation of center
        assert(not vertical_center.isnan())
        assert(not horizont_center.isnan())

        quad_center = vertical_center.avg(horizont_center)
        quad_center.coords = u_quad_center, v_quad_center

        ### 4. solve the spoke at relative_position if close. Subdivide otherwise.
        half_dist = 0.5 * p_lambda

        if abs(u - half_dist) <= epsilon and abs(v-half_dist) <= epsilon:
            # close the quad center
            result_spoke = quad_center
        elif abs(u) <= epsilon and abs(v) <= epsilon:
            # close to the left top corner
            result_spoke = sp11
        elif abs(u - p_lambda) < epsilon and abs(v) < epsilon:
            # close to the left bot corner
            result_spoke = sp21
        elif abs(v - p_lambda) < epsilon and abs(u) < epsilon:
            # close to the right top corner
            result_spoke = sp12
        elif abs(v - p_lambda) < epsilon and abs(u - p_lambda) < epsilon:
            # right bot 
            result_spoke = sp22
        elif abs(u - half_dist) <= epsilon and abs(v) < epsilon:
            # left_middle_spoke
            result_spoke = left_middle_spoke
        elif abs(u) < epsilon and abs(v - half_dist) < epsilon:
            # top_middle_spoke
            result_spoke = top_middle_spoke
        elif abs(u - half_dist) < epsilon and abs(v - p_lambda) < epsilon:
            # right_middle_spoke
            result_spoke = right_middle_spoke
        elif abs(u - p_lambda) < epsilon and abs(v - half_dist) < epsilon:
            # bot_middle_spoke
            result_spoke = bot_middle_spoke
        else:
            # subdivide to approach (u, v)
            interpolated_spokes = \
                top_middle_spoke, right_middle_spoke, bot_middle_spoke, left_middle_spoke, quad_center
            result_spoke = self._new_quad_interpolation(relative_position, corner_spokes, interpolated_spokes, half_dist, p_lambda)
        return result_spoke

    def _new_quad_interpolation(self, relative_position, prime_spokes, interpolated_spokes, half_dist, p_lambda):

        u, v = relative_position
        sp11, sp12, sp21, sp22 = prime_spokes
        top_middle_spoke, right_middle_spoke, bot_middle_spoke, left_middle_spoke, center = \
                                            interpolated_spokes
        if u < half_dist and v > half_dist:
            new_corner = top_middle_spoke, center, right_middle_spoke, sp12
            new_relative_position = u, v - half_dist
            return self._interpolate_quad(new_relative_position, new_corner, p_lambda / 2)
        elif u < half_dist and v < half_dist:
            new_corner = sp11, left_middle_spoke, center, top_middle_spoke
            return self._interpolate_quad(relative_position, new_corner, p_lambda / 2)
        elif u > half_dist and v < half_dist:
            new_corner = left_middle_spoke, sp21, bot_middle_spoke, center
            new_relative_position = u - half_dist, v
            return self._interpolate_quad(new_relative_position, new_corner, p_lambda / 2)
        elif u > half_dist and v > half_dist:
            new_corner = center, bot_middle_spoke, sp22, right_middle_spoke
            new_relative_position = u - half_dist, v - half_dist
            return self._interpolate_quad(new_relative_position, new_corner, p_lambda / 2)
        else:
            # interpolate on a line segment, subdivide the segment
            if abs(v - half_dist) < epsilon:
                new_corner = top_middle_spoke, bot_middle_spoke
                new_relative_position = u
                return self._interpolate_segment(new_relative_position, new_corner, 1, is_horizontal=False)
            elif abs(u - half_dist) < epsilon:
                new_corner = left_middle_spoke, right_middle_spoke
                new_relative_position = v
                return self._interpolate_segment(new_relative_position, new_corner, 1, is_horizontal=True)

        return new_corner
    def _interpolate_segment(self, dist, end_spokes, p_lambda, is_horizontal):
        """
        p_lambda represents the total distance from start to end of the segment
        dist is the target position of the desired interpolation
        """
        start_spoke, end_spoke = end_spokes
        absolute_u_start, absolute_v_start = start_spoke.coords
        absolute_u_end, absolute_v_end = end_spoke.coords

        if is_horizontal:
            absolute_u_middle = (absolute_u_start + absolute_u_end) // 2
            absolute_v_middle = absolute_v_start
        else:
            absolute_u_middle = absolute_u_start
            absolute_v_middle = (absolute_v_start + absolute_v_end) // 2
        middle_spoke = self._interpolate_middle(start_spoke, end_spoke, p_lambda, absolute_u_middle, absolute_v_middle)

        half_dist = p_lambda / 2
        if abs(dist - half_dist) < epsilon:
            result_spoke = middle_spoke
        elif dist < half_dist:
            new_ends = start_spoke, middle_spoke
            result_spoke = self._interpolate_segment(dist, new_ends, half_dist, is_horizontal)
        elif dist > half_dist:
            new_ends = middle_spoke, end_spoke
            result_spoke = self._interpolate_segment(dist-half_dist, new_ends, half_dist, is_horizontal)
        return result_spoke
    def _interpolate_dirs_within_quad(self, hor, vert, interpolate_level, num_rows, num_cols):
        num_steps = np.power(2, interpolate_level)
        step_size = 1.0 / num_steps
        num_spokes_per_row = (num_cols - 1) * num_steps + 1
        num_spokes_per_col = (num_rows - 1) * num_steps + 1

        ## matrix 3 x row x col
        mat_horizontal_interps = np.zeros((3, num_spokes_per_col, num_spokes_per_row))
        mat_vertical_interps = np.zeros((3, num_spokes_per_col, num_spokes_per_row))

        hor_ext = np.empty((0,3))
        for r in range(num_rows):
            start_idx = r * num_spokes_per_row
            end_idx = (r + 1) * num_spokes_per_row
            first_row = hor[start_idx:end_idx, :]   # n x 3
            hor_ext = np.vstack((hor_ext, first_row))

            mat_horizontal_interps[:, r * num_steps, :] = np.transpose(first_row)

            if r == num_rows - 1:
                continue
            for s in range(1, num_steps):
                base = s + r * num_steps
                indices = [base + k * num_spokes_per_col for k in range(num_cols)]
                hor_rots_subrow = vert[indices]
                hor_times = np.arange(0, num_cols, 1)
                slerp = Slerp(hor_times, R.from_rotvec(hor_rots_subrow))
                interp_times = np.arange(0, num_cols - 1 + step_size, step_size)
                interp_hor_subrow = slerp(interp_times).as_rotvec()
                hor_ext = np.vstack((hor_ext, interp_hor_subrow))
                mat_horizontal_interps[:, base, :] = np.transpose(interp_hor_subrow)

        ver_ext = np.empty((0, 3))
        for c in range(num_cols):
            start_idx = c * num_spokes_per_col
            end_idx = (c+1) * num_spokes_per_col
            first_col = vert[start_idx:end_idx, :]
            ver_ext = np.vstack((ver_ext, first_col))
            mat_vertical_interps[:, :, c * num_steps] = np.transpose(first_col)
            if c == num_cols - 1:
                continue
            for s in range(1, num_steps):
                base = s + c * num_steps
                indices = [base + k * num_spokes_per_row for k in range(num_rows)]
                vert_rots_subrow = hor[indices]
                ver_times = np.arange(0, num_rows, 1)
                slerp = Slerp(ver_times, R.from_rotvec(vert_rots_subrow))
                interp_times = np.arange(0, num_rows - 1 + step_size, step_size)
                interp_vert_subrow = slerp(interp_times).as_rotvec()
                ver_ext = np.vstack((ver_ext, interp_vert_subrow))
                mat_vertical_interps[:, :, base] = np.transpose(interp_vert_subrow)
        ### average horizontal and vertical interpolation
        avg_dirs = (mat_vertical_interps + mat_horizontal_interps) / 2

        ### return unit direction vectors
        return avg_dirs / LA.norm(avg_dirs, axis=0, keepdims=True)
    def _interpolate_bdry_directions(self, input_srep, interpolate_level, num_rows, num_cols):
        # interpolate all directions on quads' boundaries at once

        # within every interval there are #num_steps interpolations
        num_steps = np.power(2, interpolate_level)
        step_size = 1.0 / num_steps

        interpolated_hor_directions = np.empty((0,3))
        for r in range(num_rows):
            horizontal_dirs = np.empty((0,3))
            horizontal_times = []

            for c in range(num_cols):
                idx = r * num_cols + c
                dirs_da = input_srep.GetPointData().GetArray('spokeDirection')
                curr_dir = np.array([[dirs_da.GetValue(idx * 3),\
                       dirs_da.GetValue(idx * 3 + 1), \
                       dirs_da.GetValue(idx * 3 + 2)]])
                horizontal_dirs = np.vstack((horizontal_dirs, curr_dir))
                horizontal_times.append(c)

            horizontal_rots = R.from_rotvec(horizontal_dirs)
            slerp = Slerp(horizontal_times, horizontal_rots)
            interp_times = np.arange(0, num_cols + step_size - 1, step_size)

            interp_hor_rots = slerp(interp_times).as_rotvec()
            interpolated_hor_directions = np.vstack((interpolated_hor_directions, interp_hor_rots))

        interpolated_vert_directions = np.empty((0,3))
        for c in range(num_cols):
            vertical_dirs = np.empty((0,3))
            vertical_times = []

            for r in range(num_rows):
                idx = r * num_cols + c
                dirs_da = input_srep.GetPointData().GetArray('spokeDirection')
                curr_dir = np.array([[dirs_da.GetValue(idx * 3),\
                       dirs_da.GetValue(idx * 3 + 1), \
                       dirs_da.GetValue(idx * 3 + 2)]])
                vertical_dirs = np.vstack((vertical_dirs, curr_dir))
                vertical_times.append(r)

            vertical_rots = R.from_rotvec(vertical_dirs)
            slerp = Slerp(vertical_times, vertical_rots)
            interp_times = np.arange(0, num_rows + step_size - 1, step_size)
            interp_vert_rots = slerp(interp_times).as_rotvec()
            interpolated_vert_directions = np.vstack((interpolated_vert_directions, interp_vert_rots))
        return interpolated_hor_directions, interpolated_vert_directions

    def no_interpolate(self, input_srep):
        ret_spokes = []
        for i in range(0, input_srep.GetNumberOfPoints(), 2):
            base_pt = np.array(input_srep.GetPoint(i))
            bdry_pt = np.array(input_srep.GetPoint(i+1))
            ret_spokes.append(Spoke(base_pt=base_pt, bdry_pt=bdry_pt))
        return ret_spokes
    def interpolate(self, input_srep, interpolate_level, num_rows, num_cols):
        """
        main entry of interpolation
        """
        if interpolate_level == 0:
            return self.no_interpolate(input_srep)
        # steps of interpolation
        self.input_srep = input_srep
        num_steps = np.power(2, interpolate_level)
        step_size = 1.0 / num_steps
        self.spacing = step_size

        interpolated_spokes = []
        num_crest_points = num_rows * 2 + (num_cols - 2) * 2
        num_samples_outward = 1 + num_rows // 2
        for r in range(num_crest_points):
            for c in range(num_samples_outward -1):
                # if r!=1 or c != 1:
                #     continue
                # hor, vert = self._interpolate_bdry_directions(input_srep, interpolate_level, num_rows, num_cols)
                # ### 3 x rows x cols (where rows and cols are the number after interpolation)
                # self.interpolated_dirs = self._interpolate_dirs_within_quad(hor, vert, interpolate_level, num_rows, num_cols)

                next_row = r + 1 # next radial line index
                # compute derivatives for 4 corners
                if r == num_crest_points - 1:
                    dxdu11, dxdv11 = self._compute_derivative(input_srep, r, c, num_rows, num_cols)
                    dxdu21, dxdv21 = self._compute_derivative(input_srep, 0, c, num_rows, num_cols)
                    dxdu12, dxdv12 = self._compute_derivative(input_srep, r, c+1, num_rows, num_cols)
                    dxdu22, dxdv22 = self._compute_derivative(input_srep, 0, c+1, num_rows, num_cols)
                    next_row = 0

                else:
                    dxdu11, dxdv11 = self._compute_derivative(input_srep, r, c, num_rows, num_cols)
                    dxdu21, dxdv21 = self._compute_derivative(input_srep, r+1, c, num_rows, num_cols)
                    dxdu12, dxdv12 = self._compute_derivative(input_srep, r, c+1, num_rows, num_cols)
                    dxdu22, dxdv22 = self._compute_derivative(input_srep, r+1, c+1, num_rows, num_cols)

                corner_deriv = dxdu11, dxdv11, dxdu21, dxdv21, dxdu12, dxdv12, dxdu22, dxdv22
                # Form 4 corners of interpolation quad
                idx11 = r * num_samples_outward + c
                idx21 = next_row * num_samples_outward + c
                idx22 = next_row * num_samples_outward + (c+1)
                idx12 = r * num_samples_outward + (c+1)

                radii_da = input_srep.GetPointData().GetArray('spokeLength')
                dirs_da = input_srep.GetPointData().GetArray('spokeDirection')
                r0 = radii_da.GetValue(idx11)
                r1 = radii_da.GetValue(idx21)
                r2 = radii_da.GetValue(idx22)
                r3 = radii_da.GetValue(idx12)

                u11 = [dirs_da.GetValue(idx11 * 3),\
                       dirs_da.GetValue(idx11 * 3 + 1), \
                       dirs_da.GetValue(idx11 * 3 + 2)]
                u21 = [dirs_da.GetValue(idx21 * 3), \
                       dirs_da.GetValue(idx21 * 3 + 1), \
                       dirs_da.GetValue(idx21 * 3 + 2)]
                u22 = [dirs_da.GetValue(idx22 * 3), \
                       dirs_da.GetValue(idx22 * 3 + 1), \
                       dirs_da.GetValue(idx22 * 3 + 2)]
                u12 = [dirs_da.GetValue(idx12 * 3), \
                       dirs_da.GetValue(idx12 * 3 + 1), \
                       dirs_da.GetValue(idx12 * 3 + 2)]

                ### Note that input_srep stores base points and bdry points for spokes
                base_pts_array = input_srep.GetPointData().GetArray("basePoints")
                pt11 = pt12 = pt21 = pt22 = [0] * 3
                pt11 = base_pts_array.GetTuple3(idx11)
                pt21 = base_pts_array.GetTuple3(idx21)
                pt22 = base_pts_array.GetTuple3(idx22)
                pt12 = base_pts_array.GetTuple3(idx12)
                corner_pts = pt11, pt21, pt22, pt12

                absolute_uv11 = np.array([r, c]) * num_steps
                absolute_uv21 = np.array([next_row, c]) * num_steps
                absolute_uv22 = np.array([next_row, c+1]) * num_steps
                absolute_uv12 = np.array([r, c+1]) * num_steps
                # if r == num_rows - 1 or c == num_cols -1:
                #     print('last row / col')
                #     print(absolute_uv12, absolute_uv21)
                sp11 = Spoke(r0, u11, pt11, absolute_uv11)
                sp21 = Spoke(r1, u21, pt21, absolute_uv21)
                sp22 = Spoke(r2, u22, pt22, absolute_uv22)
                sp12 = Spoke(r3, u12, pt12, absolute_uv12)
                corner_spokes = sp11, sp12, sp21, sp22

                # interpolate dirs along edges, including vertices
                # top_dirs = self._scipy_slerp(sp11.U, sp12.U, num_steps + 1)
                # bot_dirs = self._scipy_slerp(sp21.U, sp22.U, num_steps + 1)
                # left_dirs= self._scipy_slerp(sp11.U, sp21.U, num_steps + 1)
                # right_dirs=self._scipy_slerp(sp12.U, sp22.U, num_steps + 1)

#                appendFilter = vtk.vtkAppendPolyData()
                # interpolate step by step in this quad, each iteration generate one spoke
                for ri in range(num_steps + 1):
                    for ci in range(num_steps + 1):
                        # if c > 0 and ci == 0:
                        #     ## otherwise, it is repetitive on bdry of cells
                        #     continue
                        # if r > 0 and ri == 0:
                        #     continue

                        relative_position = ri * step_size, ci * step_size
                        ## interpolate spoke base points
                        interp_skeletal_pt = self._interpolate_skeleton(relative_position, corner_pts, corner_deriv)

                        ## interpolate spoke directions and radii by a successive subdivision
                        interpolated_spoke = self._interpolate_quad(relative_position, corner_spokes, 1.0)
                        interpolated_spoke.p = interp_skeletal_pt
#                        interpolated_spoke.coords = absolute_row_index, absolute_col_index
#                        if ri == 0 or ci == 0 or ri == num_steps or ci == num_steps:
#                        interpolated_spoke.U = interp_dir
                        interpolated_spokes.append(interpolated_spoke)
                        # appendFilter.AddInputData(interpolated_spoke.visualize()[0])

                        # appendFilter.Update()
                        # visualize(appendFilter.GetOutput(), input_srep)

        ### correct directions of spokes at boundaries of quads
        # for r in range(num_rows - 1):
        #     for c in range(num_cols -1):
        #         for ri in range(num_steps + 1):
        #             for ci in range(num_steps + 1):
        #                 if ri == 0 or ci == 0 or ri == num_steps or ci == num_steps:

        return interpolated_spokes
    def construct_interp_implied_surface(self, srep, interp_spokes, num_rows=5, num_cols=9, interpolate_level=3):
        """
        Connect spokes' ends into triangle mesh.
        For visualization, I add points redundently on both boundary and skeleton
        Input srep: vtkPolyData for the visualization of primary spokes
        Input interp_spokes: list of spokes that are to be connected
        Input num_rows, num_cols: of the primary s-rep
        Input interpolate_level

        Return polydata for the implied surface
        Return corresponding polydata for the skeleton

        """
        num_crest_points = num_rows * 2 + (num_cols - 2) * 2
        num_samples_outward = 1 + num_rows // 2
        num_steps = np.power(2, interpolate_level)
        id_curr = -1
        implied_surf = vtk.vtkPolyData()
        implied_surf_pt = vtk.vtkPoints()
        implied_surf_connections = vtk.vtkCellArray()

        corresponding_base_poly = vtk.vtkPolyData()
        corresponding_base_pt = vtk.vtkPoints()
        corresponding_base_connections = vtk.vtkCellArray()
        vectors_skeleton_2_wavefront = []
        ## iterate over primary spokes
        for r in range(num_crest_points):
            for c in range(num_samples_outward -1):
                ## iterate within a quad
                for ri in range(num_steps + 1):
                    for ci in range(num_steps + 1):
                        id_curr += 1

                        if ci == num_steps or (ri % num_steps == 0 and ri > 0):
                            ## skip the last step of a quad and repeated radial boundary of a quad
                            continue

                        # if r == num_crest_points - 1:
                        #     ## The last radial line should connect with the first (0th) radial line
                        #     total_num = len(interp_spokes)
                        #     s = interp_spokes[id_curr%total_num]
                        #     s_next = interp_spokes[(id_curr + 1) % total_num]
                        #     s_next_col = interp_spokes[(id_curr + num_steps + 1) % total_num]
                        #     s_4th_corner = interp_spokes[(id_curr + num_steps + 2) % total_num]
                        # else:
                        s = interp_spokes[id_curr]
                        s_next = interp_spokes[id_curr + 1]
                        s_next_col = interp_spokes[id_curr + num_steps + 1]
                        s_4th_corner = interp_spokes[id_curr + num_steps + 2]

                        id_curr_pt = implied_surf_pt.InsertNextPoint(s.getB())
                        id_next_col_pt = implied_surf_pt.InsertNextPoint(s_next_col.getB())
                        id_next_pt = implied_surf_pt.InsertNextPoint(s_next.getB())
                        id_4th_corner = implied_surf_pt.InsertNextPoint(s_4th_corner.getB())

                        base_id_curr_pt = corresponding_base_pt.InsertNextPoint(s.p)
                        base_id_next_col_pt = corresponding_base_pt.InsertNextPoint(s_next_col.p)
                        base_id_next_pt = corresponding_base_pt.InsertNextPoint(s_next.p)
                        base_id_4th_corner = corresponding_base_pt.InsertNextPoint(s_4th_corner.p)
                        triangle = vtk.vtkTriangle()
                        triangle.GetPointIds().SetId(0, id_curr_pt)
                        triangle.GetPointIds().SetId(1, id_next_pt)
                        triangle.GetPointIds().SetId(2, id_next_col_pt)

                        another_triangle = vtk.vtkTriangle()
                        another_triangle.GetPointIds().SetId(0, id_next_pt)
                        another_triangle.GetPointIds().SetId(1, id_4th_corner)
                        another_triangle.GetPointIds().SetId(2, id_next_col_pt)

                        implied_surf_connections.InsertNextCell(triangle)
                        implied_surf_connections.InsertNextCell(another_triangle)

                        corresponding_base_connections.InsertNextCell(triangle)
                        corresponding_base_connections.InsertNextCell(another_triangle)

                        #### This part is very cool debug code that shows tile steps one by one
                        # implied_surf.SetPoints(implied_surf_pt)
                        # implied_surf.SetPolys(implied_surf_connections)
                        # implied_surf.Modified()

                        # appender = vtk.vtkAppendPolyData()
                        # appender.AddInputData(s.visualize()[0])

                        # appender.AddInputData(s_next_col.visualize()[0])
                        # appender.AddInputData(s_next.visualize()[0])
                        # appender.AddInputData(s_4th_corner.visualize()[0])
                        # appender.Update()
                        # overlay_polydata(implied_surf, srep, appender.GetOutput())

        implied_surf.SetPoints(implied_surf_pt)
        implied_surf.SetPolys(implied_surf_connections)
        implied_surf.Modified()

        corresponding_base_poly.SetPoints(corresponding_base_pt)
        corresponding_base_poly.SetPolys(corresponding_base_connections)
        corresponding_base_poly.Modified()

        return implied_surf, corresponding_base_poly

# if True:# __name__ == '__main__':
#     test_data_file = '/playpen/workspace/my_paper/linking/data/hipp_init_srep/103430/header.xml'
#     up_renderable, down_renderable, crest_renderable, num_rows, num_cols =\
#                                                     srep_io.readSrepFromXML(test_data_file)
#     interp = Interpolater()
#     interp_up_spokes = interp.interpolate(up_renderable.data, 2, up_renderable.rows, up_renderable.cols)
#     appendFilter = vtk.vtkAppendPolyData()
#     for s in interp_up_spokes:
#         appendFilter.AddInputData(s.visualize()[0])
#     appendFilter.Update()
#     visualize(appendFilter.GetOutput(), up_renderable.data)
#     print('done')