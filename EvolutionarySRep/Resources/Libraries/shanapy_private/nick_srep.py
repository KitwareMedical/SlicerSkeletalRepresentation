"""
Nick's version of a class for an s-rep. 
"""

import math
import vtk
import pyvista as pv
import numpy as np
from scipy.interpolate import interp1d
from ellipsoid_srep_coords import ellipse_theta_tau_to_x_y, ellipsoid_xm_ym_tau_to_x_y_z
from scipy.spatial.transform import Rotation as R
from scipy.spatial.transform import Slerp
from scipy.stats import norm

class Srep(object):
    def __init__(self):
        super()
        self.polydata = None
        self.num_spine_points = None
        self.num_skeletal_onionskins = None
        self.skeleton_point_indices = None
        self.boundary_point_indices = None
        self.crest_point_indices = None
        self.fold_point_indices = None
        self.vertex_indices = None
        self.coc_indices = None

        """
        A bit of documentation about skeletal sampling and onion skins:

        Consider the skeletal points (lying on a 2D surface within the object):
        We want to get the "onion skins" inside the skeleton.
        First 10 are the spine, next 20 are the first onion skin, next 20 are the second onion skin, etc...
        That continues until you get to the onion skin which is not yet at the fold curve.
        Then the next 10 are the spine again, the next 20 are the onion skins again, etc..
        Again that continues but you don't get all the way to the fold curve.
        Finally, the last 2*(num spine points) - 2 points are the fold curve.

        So in total there will be 2*(num_spine_points) + 2*(num_steps - 2)*(2*(num_spine_points) - 2) + 2*(num_spine_points)-2
        Those three terms accounting for:
            both sets of spine points
            both sets of all intermediate "onion skins" on the skeleton
            the single set of fold curve points
        """
    
    def from_polydata(self, pd, num_spine_points=None, num_steps=4):
        """
        Initialize all the Srep object attrs by loading a file.
        Assuming that the .vtk srep object being loaded was made by the same srep-fitter
        as here.
        """
        self.polydata = pd
        print(f"loaded srep with {self.polydata.GetNumberOfPoints()} points")
        if num_spine_points is None:
            # figure out number of spine points
            num_pts = self.polydata.GetNumberOfPoints() // 2
            num_spine_points = (num_pts + 4*num_steps - 6) // (4*num_steps - 4)
            print(f"calculated number of spine points: {num_spine_points}")

        self.skeleton_point_indices = list(range(0, self.polydata.GetNumberOfPoints(), 2))
        self.boundary_point_indices = list(range(1, self.polydata.GetNumberOfPoints(), 2))
        self.num_spine_points = num_spine_points
        self.num_skeletal_onionskins = num_steps
        self.fold_point_indices = self.skeleton_point_indices[-(2*num_spine_points - 2):]
        self.crest_point_indices = [x + 1 for x in self.fold_point_indices]
        self.coc_indices = [self.fold_point_indices[0], self.fold_point_indices[num_spine_points -1]]
        self.vertex_indices = [x + 1 for x in self.coc_indices]

    def load_polydata_file(self, path, num_spine_points=None, num_steps=4):
        """
        Initialize all the Srep object attrs by loading a file.
        Assuming that the .vtk srep object being loaded was made by the same srep-fitter
        as here.
        """
        self.polydata = pv.PolyData(path)
        print(f"loaded srep with {self.polydata.GetNumberOfPoints()} points")
        if num_spine_points is None:
            # figure out number of spine points
            num_pts = self.polydata.GetNumberOfPoints() // 2
            num_spine_points = (num_pts + 4*num_steps - 6) // (4*num_steps - 4)
            print(f"calculated number of spine points in {path}: {num_spine_points}")

        self.skeleton_point_indices = list(range(0, self.polydata.GetNumberOfPoints(), 2))
        self.boundary_point_indices = list(range(1, self.polydata.GetNumberOfPoints(), 2))
        self.num_spine_points = num_spine_points
        self.num_skeletal_onionskins = num_steps
        self.fold_point_indices = self.skeleton_point_indices[-(2*num_spine_points - 2):]
        self.crest_point_indices = [x + 1 for x in self.fold_point_indices]
        self.coc_indices = [self.fold_point_indices[0], self.fold_point_indices[num_spine_points -1]]
        self.vertex_indices = [x + 1 for x in self.coc_indices]

    def fit_ellipsoid_radii(self, radii, num_spine_points = 10, spine_distribution: callable=None, spine_distribution_is_from_kde: bool=False, num_steps = 4, invcdf: callable=None):
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

        self.num_spine_points = num_spine_points
        self.num_skeletal_onionskins = num_steps

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
                # print(max(cdf_pos))
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

        if num_steps > 2:
            tau1_arr = np.concatenate([
                np.ones_like(skeleton_thetas) * t for t in np.linspace(0,1,num_steps)[1:-1]
            ])
        else: tau1_arr = np.array([])

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

        self.polydata = spokes_poly

        self.skeleton_point_indices = list(range(0, spokes_poly.GetNumberOfPoints(), 2))

        self.boundary_point_indices = list(range(1, spokes_poly.GetNumberOfPoints(), 2))

        self.crest_point_indices = self.boundary_point_indices[-crest_x_arr.shape[0]:]
        s = len(self.crest_point_indices)//2
        self.crest_point_indices = self.crest_point_indices[:s+1] + self.crest_point_indices[::-1][0:s-1]

        self.fold_point_indices = self.skeleton_point_indices[-crest_x_arr.shape[0]:]
        s = len(self.fold_point_indices)//2
        self.fold_point_indices = self.fold_point_indices[:s+1] + self.fold_point_indices[::-1][0:s-1]

        self.vertex_indices = vertex_indices

        self.coc_indices = coc_indices

        return spokes_poly, coc_indices, vertex_indices

    def get_points_from_index_list(self, l):
        if self.polydata is None:
            return None
        points = np.array(
            [self.polydata.GetPoint(i) for i in l]
        )
        return points

    def plotter(self, mesh=None, **kwargs):
        if self.polydata is None:
            return None
        plt = pv.Plotter(kwargs)
        if mesh is not None:
            plt.add_mesh(mesh, color='white', opacity=0.2)

        plt.add_axes()
        plt.add_mesh(self.polydata, color='orange',line_width=2, opacity=0.2)

        # fold_pts = self.get_points_from_index_list(self.fold_point_indices + [self.fold_point_indices[0]])
        # fold_curve = pv.Spline(fold_pts, 1000)
        # plt.add_mesh(fold_curve, line_width=4, color='cornflowerblue')

        # crest_pts = self.get_points_from_index_list(self.crest_point_indices + [self.crest_point_indices[0]])
        # crest_curve = pv.Spline(crest_pts, 1000)
        # plt.add_mesh(crest_curve, line_width=4, color='red')

        return plt
    
    def __interpolate_spoke_direction(self, skeletal_pt_1, boundary_pt_1, skeletal_pt_2, boundary_pt_2, t):
        """
        Interpolates a spoke between two other spokes.
        Input spokes defined by skeletal (medial) and boundary points.
        Another input, t, gives you the weighting of the input spoke directions.

        Only worrying about spoke direction now since the base point
        is taken care of, and spoke length is also taken care of.
        """

        assert(0 <= t <= 1)
        
        # skeletal_pt_1 = np.array(skeletal_pt_1)
        # boundary_pt_1 = np.array(boundary_pt_1)
        # skeletal_pt_2 = np.array(skeletal_pt_2)
        # boundary_pt_2 = np.array(boundary_pt_2)

        vec1 = boundary_pt_1 - skeletal_pt_1
        vec2 = boundary_pt_2 - skeletal_pt_2

        vec1 /= np.linalg.norm(vec1)
        vec2 /= np.linalg.norm(vec2)

        v = np.cross(vec1, vec2)
        c = np.sum(vec1 * vec2)
        CPM = np.cross(v, np.eye(3) * -1) # cross product matrix

        # rotation matrix that takes vec1 to vec2
        rotmat = np.eye(3) + CPM + ((CPM @ CPM) * 1/(1+c))

        rotations = R.from_matrix([np.eye(3),rotmat])
        key_times = [0,1]
        slerp = Slerp(key_times, rotations)

        times = [t]
        interp_rots = slerp(times).as_matrix()
        interp_vec = interp_rots[0] @ vec1

        return interp_vec

    def __regularize_points_in_track(self, pts):
        """
        Take an array of points (N x 3) and 
        return a new array of points (N x 3)
        where each point has been adjusted along 
        a "track" defined by the original points
        so that the resulting points are all equidistant.

        The track is piecewise linear.

        Also returns the lengths of the adjustments 
        as a fraction of the length of the track.

        Adjustments are calculated for the spine, and applied
        to the other tracks on the onion skins and fold curve.
        """
        newpts = np.zeros_like(pts)
        indices_and_partial_lens = []

        # calculate total length
        total_length = 0
        # record current lengths along the track
        current_lengths = [0]
        for i in range(pts.shape[0] - 1):
            total_length += np.linalg.norm(pts[i+1] - pts[i])
            current_lengths.append(total_length)
        
        adjustments = []
        new_lengths = np.linspace(0,total_length,len(pts))

        adjustments = [(x - y)/total_length for x,y in zip(new_lengths, current_lengths)]
        adjustments[0] = 0
        adjustments[-1] = 0

        # go through each point and see if it needs to be adjusted
        # if it does, then calculate what point lies along the track
        # at the correct length, and move this point there.
        # Also record the adjustment that was made
        for i, l in enumerate(new_lengths):
            if i == 0:
                newpts[i] = pts[i]
                indices_and_partial_lens.append((0, 0.))
                continue

            if i == len(current_lengths) - 1:
                newpts[i] = pts[i]
                indices_and_partial_lens.append((len(current_lengths) - 2, 1.))
                continue

            if adjustments[i] == 0:
                newpts[i] = pts[i]
                indices_and_partial_lens.append((i, 0.))
                continue

            # calculate where the new point should go
            for j, _ in enumerate(current_lengths):
                if current_lengths[j] < l and current_lengths[j+1] >= l:
                    remainder = l - current_lengths[j]
                    newpt = pts[j] + (pts[j+1] - pts[j])*(remainder / (current_lengths[j+1] - current_lengths[j]))
                    newpts[i] = newpt

                    partial_dist = remainder / (current_lengths[j+1] - current_lengths[j])
                    indices_and_partial_lens.append((j, partial_dist))
                    break

        return newpts, adjustments, indices_and_partial_lens

    def __adjust_points_in_track(self, pts, adjustments):
        """
        Take a set of points (N x 3) where you can form a
        "track" on which they all lie. Then apply a set of
        adjustments (perturbations along the track) which come
        as proportions of the track size.

        Adjustments are calculated from the spine, and applied 
        to the tracks in the onion skins and fold curve.

        Return the new points after adjustments, and
        information about where those points were moved to
        along the track: indices and partial lens.

        Indices and partial lens: each adjusted point is moved along its 
        track to a place that will lie between two "old" points.
        The index of that first "old" point is kept track of, and the proportional
        length of the adjustment (between 0 and 1) between those two
        "old" points is kept track of.
        """

        newpts = np.zeros_like(pts)
        indices_and_partial_lens = []

        # calculate total length
        total_length = 0
        # record current lengths along the track
        current_lengths = [0]
        for i in range(pts.shape[0] - 1):
            total_length += np.linalg.norm(pts[i+1] - pts[i])
            current_lengths.append(total_length)

        
        # print(f'{current_lengths=}')
        
        # scale adjustments by total length
        adjustments = 0.5*np.array(adjustments) * total_length
        # print(f'{adjustments=}')

        # compute new lengths
        new_lengths = np.array(current_lengths) + adjustments
        # print(f'{new_lengths=}')



        for i, l in enumerate(new_lengths):
            if i == 0:
                newpts[i] = pts[i]
                indices_and_partial_lens.append((0, 0.))
                continue

            if i == len(current_lengths) - 1:
                newpts[i] = pts[i]
                indices_and_partial_lens.append((len(current_lengths) - 2, 1.))
                continue

            if adjustments[i] == 0:
                newpts[i] = pts[i]
                indices_and_partial_lens.append((i, 0.))
                continue

            # calculate where the new point should go
            for j, _ in enumerate(current_lengths):
                try:
                    _ = current_lengths[j+1]
                except:
                    print("------")
                    print(i,l,j)
                    print(current_lengths)
                    print(adjustments)
                    print(new_lengths)
                    quit()
                if current_lengths[j] < l and current_lengths[j+1] >= l:
                    remainder = l - current_lengths[j]
                    newpt = pts[j] + (pts[j+1] - pts[j])*(remainder / (current_lengths[j+1] - current_lengths[j]))

                    newpts[i] = newpt

                    partial_dist = remainder / (current_lengths[j+1] - current_lengths[j])
                    indices_and_partial_lens.append((j, partial_dist))
                    break

        return newpts, indices_and_partial_lens
    
    def regularize_crest(self,factor=0.5):
        # get spine inds and points
        n = self.num_spine_points
        m = self.num_skeletal_onionskins

        print(n)
        print(m)

        # top and bottom here refers to the top side (+z) and 
        # bottom size (-z) of the skeletal sheet

        # Note: no top or bottom for the fold points

        top_spine_inds = self.skeleton_point_indices[:n]
        print(top_spine_inds[1:-1])
        start = (1+2*(m-2))*n - 2*(m-2) # where all the bottom indices start
        bottom_spine_inds = self.skeleton_point_indices[start:start+n]

        top_skeletal_onion_skins_inds = []
        bottom_skeletal_onion_skins_inds = []

        fold_curve_inds = self.skeleton_point_indices[2*start:2*start + 2*n-2]
        crest_curve_inds = self.boundary_point_indices[2*start:2*start + 2*n-2]

        # Move the crest points out by half a spoke length
        for _,(fold,crest) in enumerate(zip(fold_curve_inds,crest_curve_inds)):
            fold_point = np.array(self.polydata.GetPoints().GetPoint(fold))
            crest_point = np.array(self.polydata.GetPoints().GetPoint(crest))

            v = crest_point - fold_point
            nfp = fold_point + factor*v
            self.polydata.GetPoints().SetPoint(fold,nfp)

        # Now we need to collect indices along each vein
        # inds = [256,86,38,12,62,110,280]
        # new_points, adjustments, inds_and_lens = self.__regularize_points_in_track(self.get_points_from_index_list(inds))
        # for i,(ii,pt) in enumerate(zip(inds,new_points)):
        #     self.polydata.GetPoints().SetPoint(ii,pt)

        #     if adjustments[i] == 0:
        #         continue

        #     ind,len = inds_and_lens[i]
        #     behind_pt = np.array(self.polydata.GetPoint(inds[ind]))
        #     behind_bdry = np.array(self.polydata.GetPoint(inds[ind]+1))
        #     ahead_pt = np.array(self.polydata.GetPoint(inds[ind+1]))
        #     ahead_bdry = np.array(self.polydata.GetPoint(inds[ind+1]+1))

        #     new_dir = self.__interpolate_spoke_direction(behind_pt,behind_bdry,ahead_pt,ahead_bdry,len)
        #     new_len = (np.linalg.norm(behind_bdry-behind_pt) + np.linalg.norm(ahead_bdry - ahead_pt)) / 2

        #     self.polydata.GetPoints().SetPoint(ii+1,pt+new_len*new_dir)

        # work on spine end individually
        inds = [268,98,50,24]
        new_points, adjustments, inds_and_lens = self.__regularize_points_in_track(self.get_points_from_index_list(inds))
        print(adjustments)
        for i,(ii,pt) in enumerate(zip(inds,new_points)):
            self.polydata.GetPoints().SetPoint(ii,pt)

            if adjustments[i] == 0:
                continue

            ind,len = inds_and_lens[i]
            behind_pt = np.array(self.polydata.GetPoint(inds[ind]))
            behind_bdry = np.array(self.polydata.GetPoint(inds[ind]+1))
            ahead_pt = np.array(self.polydata.GetPoint(inds[ind+1]))
            ahead_bdry = np.array(self.polydata.GetPoint(inds[ind+1]+1))

            new_dir = self.__interpolate_spoke_direction(behind_pt,behind_bdry,ahead_pt,ahead_bdry,len)
            new_len = (np.linalg.norm(behind_bdry-behind_pt) + np.linalg.norm(ahead_bdry - ahead_pt)) / 2

            self.polydata.GetPoints().SetPoint(ii+1,pt+new_len*new_dir)

        inds = [244,74,26,0]
        new_points, adjustments, inds_and_lens = self.__regularize_points_in_track(self.get_points_from_index_list(inds))
        print(adjustments)
        for i,(ii,pt) in enumerate(zip(inds,new_points)):
            self.polydata.GetPoints().SetPoint(ii,pt)

            if adjustments[i] == 0:
                continue

            ind,len = inds_and_lens[i]
            behind_pt = np.array(self.polydata.GetPoint(inds[ind]))
            behind_bdry = np.array(self.polydata.GetPoint(inds[ind]+1))
            ahead_pt = np.array(self.polydata.GetPoint(inds[ind+1]))
            ahead_bdry = np.array(self.polydata.GetPoint(inds[ind+1]+1))

            new_dir = self.__interpolate_spoke_direction(behind_pt,behind_bdry,ahead_pt,ahead_bdry,len)
            new_len = (np.linalg.norm(behind_bdry-behind_pt) + np.linalg.norm(ahead_bdry - ahead_pt)) / 2

            self.polydata.GetPoints().SetPoint(ii+1,pt+new_len*new_dir)

        for si in top_spine_inds[1:-1]:
            inds = [244+si,74+si,26+si,si,si+50,si+98,268+si]
            new_points, adjustments, inds_and_lens = self.__regularize_points_in_track(self.get_points_from_index_list(inds))

            for i,(ii,pt) in enumerate(zip(inds,new_points)):
                self.polydata.GetPoints().SetPoint(ii,pt)

                if adjustments[i] == 0:
                    continue

                ind,len = inds_and_lens[i]
                behind_pt = np.array(self.polydata.GetPoint(inds[ind]))
                behind_bdry = np.array(self.polydata.GetPoint(inds[ind]+1))
                ahead_pt = np.array(self.polydata.GetPoint(inds[ind+1]))
                ahead_bdry = np.array(self.polydata.GetPoint(inds[ind+1]+1))

                new_dir = self.__interpolate_spoke_direction(behind_pt,behind_bdry,ahead_pt,ahead_bdry,len)
                new_len = (np.linalg.norm(behind_bdry-behind_pt) + np.linalg.norm(ahead_bdry - ahead_pt)) / 2

                self.polydata.GetPoints().SetPoint(ii+1,pt+new_len*new_dir)

            

        for si in bottom_spine_inds[1:-1]:
            inds = [122+si,74+si,26+si,si,si+50,si+98,146+si]
            new_points, adjustments, inds_and_lens = self.__regularize_points_in_track(self.get_points_from_index_list(inds))
            for i,(ii,pt) in enumerate(zip(inds,new_points)):
                self.polydata.GetPoints().SetPoint(ii,pt)

                if adjustments[i] == 0:
                    continue
                
                ind,len = inds_and_lens[i]
                behind_pt = np.array(self.polydata.GetPoint(inds[ind]))
                behind_bdry = np.array(self.polydata.GetPoint(inds[ind]+1))
                ahead_pt = np.array(self.polydata.GetPoint(inds[ind+1]))
                ahead_bdry = np.array(self.polydata.GetPoint(inds[ind+1]+1))

                new_dir = self.__interpolate_spoke_direction(behind_pt,behind_bdry,ahead_pt,ahead_bdry,len)
                new_len = (np.linalg.norm(behind_bdry-behind_pt) + np.linalg.norm(ahead_bdry - ahead_pt)) / 2

                self.polydata.GetPoints().SetPoint(ii+1,pt+new_len*new_dir)

        

    def regularize_skeleton(self):
        """
        First, create "track" connecting all the spine points, then
        move interior spine points along the "track" so that
        they are equidistant from each other according to distances
        on the track.

        Then adjust the skeletal onion skin points and fold points
        similarly to how the spine was adjusted, along their own "tracks". 
        Use the adjustments of the spine and apply those to the onion skins
        along their own respective tracks.

        After adjustment, skeletal points fall into a new location between
        two of the "old" points on the track. Interpolate a new spoke
        based on the old spokes of those old points. 

        Sanity check note: spine points are ordered in ascending x coordinate.
        """
        # get spine inds and points
        n = self.num_spine_points
        m = self.num_skeletal_onionskins

        # top and bottom here refers to the top side (+z) and 
        # bottom size (-z) of the skeletal sheet

        # Note: no top or bottom for the fold points

        top_spine_inds = self.skeleton_point_indices[:n]
        start = (1+2*(m-2))*n - 2*(m-2) # where all the bottom indices start
        bottom_spine_inds = self.skeleton_point_indices[start:start+n]

        top_skeletal_onion_skins_inds = []
        bottom_skeletal_onion_skins_inds = []

        fold_curve_inds = self.skeleton_point_indices[2*start:2*start + 2*n-2]

        for i in range(m-2):
            # collect inds of skeletal onion skins, not spine or fold
            # includes "right" and "left" sides of ellipsoid
            # right side (+y) includes hubs that are inbetween spine endpoint and fold curve vertex
            # left side (-y) does not include these "vertex" hubs
            # right side is comprised of first n elements,
            # left side is the remaining n-2 elts
            top_onion_skin_inds = self.skeleton_point_indices[(1+2*(i))*n-(2*i) : (1+2*(i+1))*n-(2*(i+1))]
            bottom_onion_skin_inds = self.skeleton_point_indices[start + (1+2*(i))*n-(2*i) : start + (1+2*(i+1))*n-(2*(i+1))]
            top_skeletal_onion_skins_inds.append(top_onion_skin_inds)
            bottom_skeletal_onion_skins_inds.append(bottom_onion_skin_inds)

        # print(f'top_skeletal_onion_skins_inds: {top_skeletal_onion_skins_inds}' )
        # print(f'bottom_skeletal_onion_skins_inds: {bottom_skeletal_onion_skins_inds}' )

        # first get adjustments from spine curves.

        new_spine_points, adjustments, spine_inds_and_lens = self.__regularize_points_in_track(self.get_points_from_index_list(top_spine_inds))
        # print(f'{adjustments=}')

        # calculate interpolated spine boundary points
        new_top_boundary_pts = []
        new_bottom_boundary_pts = []
        for k, c in enumerate(zip(top_spine_inds, bottom_spine_inds)):
            # index of top and bottom spine pt
            topi, boti = c

            # no need to do any interpolation if adjustment is 0
            if adjustments[k] == 0:
                interpolated_top_boundary_pt = self.polydata.GetPoint(topi + 1)
                interpolated_bottom_boundary_pt = self.polydata.GetPoint(boti + 1)
                new_top_boundary_pts.append(interpolated_top_boundary_pt)
                new_bottom_boundary_pts.append(interpolated_bottom_boundary_pt)
                continue

            # interpolate spoke directions  
            # top
            track_idx, partial_track_len = spine_inds_and_lens[k]
            behind_pt_idx = top_spine_inds[track_idx]
            ahead_pt_idx = top_spine_inds[track_idx + 1]
            behind_skeleton_pt = np.array(self.polydata.GetPoint(behind_pt_idx))
            behind_boundary_pt = np.array(self.polydata.GetPoint(behind_pt_idx + 1))
            ahead_skeleton_pt = np.array(self.polydata.GetPoint(ahead_pt_idx))
            ahead_boundary_pt = np.array(self.polydata.GetPoint(ahead_pt_idx + 1))
            interpolated_direction = self.__interpolate_spoke_direction(
                behind_skeleton_pt,
                behind_boundary_pt,
                ahead_skeleton_pt,
                ahead_boundary_pt,
                partial_track_len,
            )
            interpolated_spoke_length = (np.linalg.norm(behind_boundary_pt - behind_skeleton_pt) + np.linalg.norm(ahead_boundary_pt - ahead_skeleton_pt)) / 2
            interpolated_top_boundary_pt = new_spine_points[k] + interpolated_direction * interpolated_spoke_length

            # bot
            track_idx, partial_track_len = spine_inds_and_lens[k]
            behind_pt_idx = bottom_spine_inds[track_idx]
            ahead_pt_idx = bottom_spine_inds[track_idx + 1]
            behind_skeleton_pt = np.array(self.polydata.GetPoint(behind_pt_idx))
            behind_boundary_pt = np.array(self.polydata.GetPoint(behind_pt_idx + 1))
            ahead_skeleton_pt = np.array(self.polydata.GetPoint(ahead_pt_idx))
            ahead_boundary_pt = np.array(self.polydata.GetPoint(ahead_pt_idx + 1))
            interpolated_direction = self.__interpolate_spoke_direction(
                behind_skeleton_pt,
                behind_boundary_pt,
                ahead_skeleton_pt,
                ahead_boundary_pt,
                partial_track_len,
            )
            interpolated_spoke_length = (np.linalg.norm(behind_boundary_pt - behind_skeleton_pt) + np.linalg.norm(ahead_boundary_pt - ahead_skeleton_pt)) / 2
            interpolated_bottom_boundary_pt = new_spine_points[k] + interpolated_direction * interpolated_spoke_length

            new_top_boundary_pts.append(interpolated_top_boundary_pt)
            new_bottom_boundary_pts.append(interpolated_bottom_boundary_pt)

        # relocate all the spine points (skeleton and boundary)

        for k, c in enumerate(zip(top_spine_inds, bottom_spine_inds)):
            # index of top and bottom spine pt
            if adjustments[k] == 0:
                continue

            topi, boti = c

            # move skeleton spine points
            self.polydata.GetPoints().SetPoint(topi, new_spine_points[k])
            self.polydata.GetPoints().SetPoint(boti, new_spine_points[k])

            # move boundary spine points
            self.polydata.GetPoints().SetPoint(topi + 1, new_top_boundary_pts[k])
            self.polydata.GetPoints().SetPoint(boti + 1, new_bottom_boundary_pts[k])

        # # get new regularized interior skeletal points and update

        for toplist, botlist in zip(top_skeletal_onion_skins_inds, bottom_skeletal_onion_skins_inds):
            top_right = toplist[:n]
            top_left = toplist[n:]
            top_left = [top_right[0]] + top_left + [top_right[-1]]

            bottom_right = botlist[:n]
            bottom_left = botlist[n:]
            bottom_left = [bottom_right[0]] + bottom_left + [bottom_right[-1]]

            top_right_pts = self.get_points_from_index_list(top_right)
            # print(f"{top_right=}")
            top_left_pts = self.get_points_from_index_list(top_left)
            # print(f"{top_left=}")

            new_right_pts, right_track_inds_and_lens = self.__adjust_points_in_track(top_right_pts, adjustments)
            # print(f'top_right_pts: {top_right_pts}')
            new_left_pts, left_track_inds_and_lens = self.__adjust_points_in_track(top_left_pts, adjustments)
            # print(f'top_right_pts: {top_left_pts}')

            # calculate interpolated skeletal boundary points
            new_top_left_boundary_pts = []
            new_top_right_boundary_pts = []
            new_bottom_left_boundary_pts = []
            new_bottom_right_boundary_pts = []
            for k, c in enumerate(zip(top_right, top_left, bottom_right, bottom_left)):
                # index of top and bottom left and right skeletal points
                tr, tl, br, bl = c

                # no need to do any interpolation if adjustment is 0
                if adjustments[k] == 0:
                    interpolated_top_right_boundary_pt = self.polydata.GetPoint(tr + 1)
                    interpolated_top_left_boundary_pt = self.polydata.GetPoint(tl + 1)
                    interpolated_bottom_right_boundary_pt = self.polydata.GetPoint(br + 1)
                    interpolated_bottom_left_boundary_pt = self.polydata.GetPoint(bl + 1)

                    new_top_right_boundary_pts.append(interpolated_top_right_boundary_pt)
                    new_top_left_boundary_pts.append(interpolated_top_left_boundary_pt)
                    new_bottom_right_boundary_pts.append(interpolated_bottom_right_boundary_pt)
                    new_bottom_left_boundary_pts.append(interpolated_bottom_left_boundary_pt)
                    continue

                # interpolate spoke directions  
                # top right
                track_idx, partial_track_len = right_track_inds_and_lens[k]
                behind_pt_idx = top_right[track_idx]
                ahead_pt_idx = top_right[track_idx + 1]
                behind_skeleton_pt = np.array(self.polydata.GetPoint(behind_pt_idx))
                behind_boundary_pt = np.array(self.polydata.GetPoint(behind_pt_idx + 1))
                ahead_skeleton_pt = np.array(self.polydata.GetPoint(ahead_pt_idx))
                ahead_boundary_pt = np.array(self.polydata.GetPoint(ahead_pt_idx + 1))
                interpolated_direction = self.__interpolate_spoke_direction(
                    behind_skeleton_pt,
                    behind_boundary_pt,
                    ahead_skeleton_pt,
                    ahead_boundary_pt,
                    partial_track_len,
                )
                interpolated_spoke_length = (np.linalg.norm(behind_boundary_pt - behind_skeleton_pt) + np.linalg.norm(ahead_boundary_pt - ahead_skeleton_pt)) / 2
                interpolated_top_right_boundary_pt = new_right_pts[k] + interpolated_direction * interpolated_spoke_length

                # top left
                track_idx, partial_track_len = left_track_inds_and_lens[k]
                behind_pt_idx = top_left[track_idx]
                ahead_pt_idx = top_left[track_idx + 1]
                behind_skeleton_pt = np.array(self.polydata.GetPoint(behind_pt_idx))
                behind_boundary_pt = np.array(self.polydata.GetPoint(behind_pt_idx + 1))
                ahead_skeleton_pt = np.array(self.polydata.GetPoint(ahead_pt_idx))
                ahead_boundary_pt = np.array(self.polydata.GetPoint(ahead_pt_idx + 1))
                interpolated_direction = self.__interpolate_spoke_direction(
                    behind_skeleton_pt,
                    behind_boundary_pt,
                    ahead_skeleton_pt,
                    ahead_boundary_pt,
                    partial_track_len,
                )
                interpolated_spoke_length = (np.linalg.norm(behind_boundary_pt - behind_skeleton_pt) + np.linalg.norm(ahead_boundary_pt - ahead_skeleton_pt)) / 2
                interpolated_top_left_boundary_pt = new_left_pts[k] + interpolated_direction * interpolated_spoke_length

                # bottom right
                track_idx, partial_track_len = right_track_inds_and_lens[k]
                behind_pt_idx = bottom_right[track_idx]
                ahead_pt_idx = bottom_right[track_idx + 1]
                behind_skeleton_pt = np.array(self.polydata.GetPoint(behind_pt_idx))
                behind_boundary_pt = np.array(self.polydata.GetPoint(behind_pt_idx + 1))
                ahead_skeleton_pt = np.array(self.polydata.GetPoint(ahead_pt_idx))
                ahead_boundary_pt = np.array(self.polydata.GetPoint(ahead_pt_idx + 1))
                interpolated_direction = self.__interpolate_spoke_direction(
                    behind_skeleton_pt,
                    behind_boundary_pt,
                    ahead_skeleton_pt,
                    ahead_boundary_pt,
                    partial_track_len,
                )
                interpolated_spoke_length = (np.linalg.norm(behind_boundary_pt - behind_skeleton_pt) + np.linalg.norm(ahead_boundary_pt - ahead_skeleton_pt)) / 2
                interpolated_bottom_right_boundary_pt = new_right_pts[k] + interpolated_direction * interpolated_spoke_length

                # bottom left
                track_idx, partial_track_len = left_track_inds_and_lens[k]
                behind_pt_idx = bottom_left[track_idx]
                ahead_pt_idx = bottom_left[track_idx + 1]
                behind_skeleton_pt = np.array(self.polydata.GetPoint(behind_pt_idx))
                behind_boundary_pt = np.array(self.polydata.GetPoint(behind_pt_idx + 1))
                ahead_skeleton_pt = np.array(self.polydata.GetPoint(ahead_pt_idx))
                ahead_boundary_pt = np.array(self.polydata.GetPoint(ahead_pt_idx + 1))
                interpolated_direction = self.__interpolate_spoke_direction(
                    behind_skeleton_pt,
                    behind_boundary_pt,
                    ahead_skeleton_pt,
                    ahead_boundary_pt,
                    partial_track_len,
                )
                interpolated_spoke_length = (np.linalg.norm(behind_boundary_pt - behind_skeleton_pt) + np.linalg.norm(ahead_boundary_pt - ahead_skeleton_pt)) / 2
                interpolated_bottom_left_boundary_pt = new_left_pts[k] + interpolated_direction * interpolated_spoke_length

                new_top_right_boundary_pts.append(interpolated_top_right_boundary_pt)
                new_top_left_boundary_pts.append(interpolated_top_left_boundary_pt)
                new_bottom_right_boundary_pts.append(interpolated_bottom_right_boundary_pt)
                new_bottom_left_boundary_pts.append(interpolated_bottom_left_boundary_pt)
                
            # now actually change the mesh to the new points
            for k, c in enumerate(zip(top_right, top_left, bottom_right, bottom_left)):
                tr, tl, br, bl = c

                # skeletal pts
                self.polydata.GetPoints().SetPoint(tr, new_right_pts[k])
                self.polydata.GetPoints().SetPoint(tl, new_left_pts[k])
                self.polydata.GetPoints().SetPoint(br, new_right_pts[k])
                self.polydata.GetPoints().SetPoint(bl, new_left_pts[k])

                # boundary pts
                self.polydata.GetPoints().SetPoint(tr + 1, new_top_right_boundary_pts[k])
                self.polydata.GetPoints().SetPoint(tl + 1, new_top_left_boundary_pts[k])
                self.polydata.GetPoints().SetPoint(br + 1, new_bottom_right_boundary_pts[k])
                self.polydata.GetPoints().SetPoint(bl + 1, new_bottom_left_boundary_pts[k])
            
        # get new regularized fold points
                

        right_fold_inds = fold_curve_inds[:n]
        print(right_fold_inds)
        left_fold_inds = fold_curve_inds[n:]
        left_fold_inds = [right_fold_inds[0]] + left_fold_inds + [right_fold_inds[-1]]

        right_fold_pts = self.get_points_from_index_list(right_fold_inds)
        left_fold_pts = self.get_points_from_index_list(left_fold_inds)

        new_right_fold_pts, right_track_inds_and_lens = self.__adjust_points_in_track(right_fold_pts, adjustments)
        new_left_fold_pts, left_track_inds_and_lens = self.__adjust_points_in_track(left_fold_pts, adjustments)

        # calculate interpolated fold boundary (crest) points
        new_left_crest_pts = []
        new_right_crest_pts = []
        for k, c in enumerate(zip(right_fold_inds, left_fold_inds)):
            # index of top and bottom left and right skeletal points
            r, l = c

            # no need to do any interpolation if adjustment is 0
            if adjustments[k] == 0:
                interpolated_right_crest_pt = self.polydata.GetPoint(r + 1)
                interpolated_left_crest_pt = self.polydata.GetPoint(l + 1)

                new_right_crest_pts.append(interpolated_right_crest_pt)
                new_left_crest_pts.append(interpolated_left_crest_pt)
                continue

            # interpolate spoke directions  
            # right
            track_idx, partial_track_len = right_track_inds_and_lens[k]
            behind_pt_idx = right_fold_inds[track_idx]
            ahead_pt_idx = right_fold_inds[track_idx + 1]
            behind_skeleton_pt = np.array(self.polydata.GetPoint(behind_pt_idx))
            behind_boundary_pt = np.array(self.polydata.GetPoint(behind_pt_idx + 1))
            ahead_skeleton_pt = np.array(self.polydata.GetPoint(ahead_pt_idx))
            ahead_boundary_pt = np.array(self.polydata.GetPoint(ahead_pt_idx + 1))
            interpolated_direction = self.__interpolate_spoke_direction(
                behind_skeleton_pt,
                behind_boundary_pt,
                ahead_skeleton_pt,
                ahead_boundary_pt,
                partial_track_len,
            )
            interpolated_spoke_length = (np.linalg.norm(behind_boundary_pt - behind_skeleton_pt) + np.linalg.norm(ahead_boundary_pt - ahead_skeleton_pt)) / 2
            interpolated_right_crest_pt = new_right_fold_pts[k] + interpolated_direction * interpolated_spoke_length

            # left
            track_idx, partial_track_len = left_track_inds_and_lens[k]
            behind_pt_idx = left_fold_inds[track_idx]
            ahead_pt_idx = left_fold_inds[track_idx + 1]
            behind_skeleton_pt = np.array(self.polydata.GetPoint(behind_pt_idx))
            behind_boundary_pt = np.array(self.polydata.GetPoint(behind_pt_idx + 1))
            ahead_skeleton_pt = np.array(self.polydata.GetPoint(ahead_pt_idx))
            ahead_boundary_pt = np.array(self.polydata.GetPoint(ahead_pt_idx + 1))
            interpolated_direction = self.__interpolate_spoke_direction(
                behind_skeleton_pt,
                behind_boundary_pt,
                ahead_skeleton_pt,
                ahead_boundary_pt,
                partial_track_len,
            )
            interpolated_spoke_length = (np.linalg.norm(behind_boundary_pt - behind_skeleton_pt) + np.linalg.norm(ahead_boundary_pt - ahead_skeleton_pt)) / 2
            interpolated_left_crest_pt = new_left_fold_pts[k] + interpolated_direction * interpolated_spoke_length

            new_right_crest_pts.append(interpolated_right_crest_pt)
            new_left_crest_pts.append(interpolated_left_crest_pt)

        for k, c in enumerate(zip(right_fold_inds, left_fold_inds)):
            r, l = c

            self.polydata.GetPoints().SetPoint(r, new_right_fold_pts[k])
            self.polydata.GetPoints().SetPoint(l, new_left_fold_pts[k])

            self.polydata.GetPoints().SetPoint(r + 1, new_right_crest_pts[k])
            self.polydata.GetPoints().SetPoint(l + 1, new_left_crest_pts[k])

        self.polydata.Modified()

    def medialize_skeleton(self):

        """
        Adjust skeletal points so that they become more medial.
        Does not include fold points.

        First collect the indices of all the skeletal points that need
        to be adjusted. Collect: index of skeletal point,
        index of its double, index of top boundary point, index of
        bottom boundary point.

        Record each perturbation as well as its location in space.

        Apply all perturbations to all skeletal points in a weighted
        sum according to the distance between the perturbation and the
        skeletal point in question.

        Weights determined by gaussian kernel with std dev = 1/4 of
        second-largest bounding box dimension.

        Note: for a skeletal point at index i, its boundary point
        is at index i+1.
        """

        n = self.num_spine_points
        m = self.num_skeletal_onionskins

        perturbations = [] # xyz of location, then xyz of perturbation vec

        top_spine_inds = self.skeleton_point_indices[:n]
        start = (1+2*(m-2))*n - 2*(m-2)
        bottom_spine_inds = self.skeleton_point_indices[start:start+n]

        top_skeletal_onion_skins_inds = []
        bottom_skeletal_onion_skins_inds = []

        for i in range(m-2):
            # collect inds of onion skins
            # includes "right" and "left" sides of ellipsoid
            # right side (+y) includes hubs that are inbetween spine endpoint and fold curve vertex
            # left side (-y) does not include these "vertex" hubs
            # right side is comprised of first n elements,
            # left side is the remaining n-2 elts
            top_onion_skin_inds = self.skeleton_point_indices[(1+2*(i))*n-(2*i):(1+2*(i+1))*n-(2*(i+1))]
            bottom_onion_skin_inds = self.skeleton_point_indices[start + (1+2*(i))*n-(2*i):start + (1+2*(i+1))*n-(2*(i+1))]
            top_skeletal_onion_skins_inds.append(top_onion_skin_inds)
            bottom_skeletal_onion_skins_inds.append(bottom_onion_skin_inds)
        
        # spine adjustments
        for top_i, bot_i in zip(top_spine_inds, bottom_spine_inds):
            m = self.polydata.GetPoint(top_i)
            assert m == self.polydata.GetPoint(bot_i)

            a = self.polydata.GetPoint(top_i + 1)
            b = self.polydata.GetPoint(bot_i + 1)

            m = np.array(m)
            a = np.array(a)
            b = np.array(b)

            # if m == ((a + b) / 2):
            #     continue

            # l is vector from b to a
            l = a - b
            # let C be the point that is the orthogonal projection of m on to B-A
            # let x be such that C = B + xL
            # solve for x here using orthogonality condition (dot product = 0)
            x = np.sum(l * (m - b)) / np.sum(l ** 2)
            adj = 0.5 - x
            new_m = m + adj * l

            # record the perturbation
            perturbations.append(np.array([ 
                m,
                adj * l,
            ]))


        # Do the same thing for non-spine non-fold skeletal points
        for top_list, bot_list in zip(top_skeletal_onion_skins_inds, bottom_skeletal_onion_skins_inds):
            for top_i, bot_i in zip(top_list, bot_list):
                m = self.polydata.GetPoint(top_i)
                # assert m == self.polydata.GetPoint(bot_i)
                print(m[0])
                print(self.polydata.GetPoint(bot_i)[0])
                assert math.isclose(m[0],self.polydata.GetPoint(bot_i)[0])
                assert math.isclose(m[1],self.polydata.GetPoint(bot_i)[1])
                assert math.isclose(m[2],self.polydata.GetPoint(bot_i)[2])

                m = np.array(m)
                a = np.array(a)
                b = np.array(b)

                # if m == ((a + b) / 2):
                #     continue

                l = a - b
                x = np.sum(l * (m - b)) / np.sum(l ** 2)
                adj = 0.5 - x
                new_m = m + adj * l

                # record the perturbation
                perturbations.append(np.array([ 
                    m,
                    adj * l,
                ]))
        

        # fold inds have zero perturbation
        fold_inds = self.skeleton_point_indices[-(2*n - 2):]
        for i in fold_inds:
            m = self.polydata.GetPoint(i)

            # record the perturbation
            perturbations.append(np.array([ 
                m,
                np.array([0,0,0]),
            ]))
        
        perturbations = np.array(perturbations)

        # compute kernel bandwith by heuristic
        skeletal_locs = np.array([self.polydata.GetPoint(i) for i in self.skeleton_point_indices])
        bbox_dims = np.max(skeletal_locs, axis=0) - np.min(skeletal_locs, axis=0)
        kernel_bandwith = sorted(bbox_dims)[1] / 4

        # now apply all perturbations to all skeletal points (spine, fold, and other skeletal points)
        kernel = norm(loc=0, scale=kernel_bandwith)
        for i in self.skeleton_point_indices:
            m = np.array(self.polydata.GetPoint(i)).reshape((1,3))
            p_locs = perturbations[:,0,:]
            ps = perturbations[:,1,:]

            dists_to_ps = np.linalg.norm(m - p_locs, axis=1)
            weights = kernel.pdf(dists_to_ps)
            if np.sum(weights) < 1e-2:
                continue
            weights /= np.sum(weights)
            weights = weights.reshape((-1,1))

            this_p = np.sum(weights * ps, axis=0)
            new_m = m[0] + this_p

            self.polydata.GetPoints().SetPoint(i, new_m)
        
        self.polydata.Modified()

    def refine_spoke_lengths(self, input_mesh):

        """
        Optimize spokes, such that
        the tips are touching the input_mesh more closely

        Don't change lengths of vertex spokes
        """

        n = self.num_spine_points
        m = self.num_skeletal_onionskins
        start = (1+2*(m-2))*n - 2*(m-2) 
        fold_curve_inds = self.skeleton_point_indices[2*start:2*start + 2*n-2]
        vertex_spoke_ids = (fold_curve_inds[0], fold_curve_inds[n-1])

        num_spokes = len(self.skeleton_point_indices)
        radii_array = np.zeros(num_spokes)
        dir_array = np.zeros((num_spokes, 3))
        base_array = np.zeros((num_spokes,3))

        ### read the parameters from s-rep
        for i in range(num_spokes):
            id_base_pt = self.skeleton_point_indices[i]
            id_bdry_pt = id_base_pt + 1
            base_pt = np.array(self.polydata.GetPoint(id_base_pt))
            bdry_pt = np.array(self.polydata.GetPoint(id_bdry_pt))

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

        import nlopt
        opt = nlopt.opt(nlopt.LN_NEWUOA, len(radii_array))
        opt.set_min_objective(obj_func)
        opt.set_maxeval(2000)
        minimizer = opt.optimize(radii_array)

        # update radii of s-rep and return the updated
        for i in range(num_spokes):
            id_base_pt = self.skeleton_point_indices[i]
            if id_base_pt in vertex_spoke_ids:
                continue # dont change lengths of vertex spokes
            id_bdry_pt = id_base_pt + 1
            base_pt = base_array[i, :]
            radius = minimizer[i]
            direction = dir_array[i, :]

            new_bdry_pt = base_pt + radius * direction
            self.polydata.GetPoints().SetPoint(id_bdry_pt, new_bdry_pt)

        self.polydata.Modified()
        return

    def orthogonalize_spokes_to_boundary(self, input_mesh, do_smoothing=False):
        """
        Change spoke directions so that they are orthogonal to the boundary
        """

        if do_smoothing:
            taubin_smooth = vtk.vtkWindowedSincPolyDataFilter()
            taubin_smooth.SetInputData(input_mesh)
            taubin_smooth.SetNumberOfIterations(20)
            taubin_smooth.BoundarySmoothingOff()
            taubin_smooth.FeatureEdgeSmoothingOff()
            taubin_smooth.SetPassBand(0.01)
            taubin_smooth.NonManifoldSmoothingOn()
            taubin_smooth.NormalizeCoordinatesOn()

            taubin_smooth.Update()
            smoothed_mesh = taubin_smooth.GetOutput()
        else:
            smoothed_mesh = input_mesh

        # print(f"num points in remesh {smoothed_mesh.GetNumberOfPoints()}")

        # get normals, and 4 measures of curvature

        normal_generator = vtk.vtkPolyDataNormals()
        normal_generator.SetInputData(smoothed_mesh)
        normal_generator.SplittingOff()
        normal_generator.ComputePointNormalsOn()
        normal_generator.ComputeCellNormalsOff()
        normal_generator.Update()
        normals = np.array(normal_generator.GetOutput().GetPointData().GetNormals())

        points = np.array([smoothed_mesh.GetPoint(i) for i in range(smoothed_mesh.GetNumberOfPoints())])


        # compute kernel bandwith by heuristic
        bbox_dims = np.max(points, axis=0) - np.min(points, axis=0)
        kernel_bandwith = sorted(bbox_dims)[0] / 3

        kernel = norm(loc=0, scale=kernel_bandwith)
        def get_smoothed_normal(x):
            x = x.reshape((1,3))

            dists_to_ps = np.linalg.norm(x - points, axis=1)
            weights = kernel.pdf(dists_to_ps)
            weights /= np.sum(weights)
            weights = weights.reshape((-1,1))

            smoothed_normal = np.sum(weights * normals, axis=0)
            smoothed_normal /= np.linalg.norm(smoothed_normal)

            return smoothed_normal
        
        # get each spoke and change the boundary pt
        n = self.num_spine_points
        m = self.num_skeletal_onionskins
        start = (1+2*(m-2))*n - 2*(m-2) 
        fold_curve_inds = self.skeleton_point_indices[2*start:2*start + 2*n-2]
        vertex_spoke_ids = (fold_curve_inds[0], fold_curve_inds[n-1])
        for id_base_pt in self.skeleton_point_indices:
            id_bdry_pt = id_base_pt + 1

            if id_base_pt in vertex_spoke_ids:
                continue # dont change vertex spokes

            base_pt = np.array(self.polydata.GetPoint(id_base_pt))
            bdry_pt = np.array(self.polydata.GetPoint(id_bdry_pt))

            radius = np.linalg.norm(bdry_pt - base_pt)
            direction = (bdry_pt - base_pt) / radius

            normal = get_smoothed_normal(bdry_pt)
            theta = np.arccos(np.sum(direction * normal))

            if math.isnan(theta):
                print(normal)
                print(direction)
                print(np.linalg.norm(normal))
                print(np.linalg.norm(direction))
                print(np.sum(normal * direction))
                raise Exception

            if abs(theta % (2*np.pi)) < np.pi/8:
                continue

            v = np.cross(normal, np.cross(normal, direction))
            v /= np.linalg.norm(v)
            new_bdry_pt = bdry_pt + (v * radius * np.sin(theta))

            self.polydata.GetPoints().SetPoint(id_bdry_pt, new_bdry_pt)

    def discretize_spokes(self, num_onionskins):
        """
        Before this method is called, we should have a srep that has 
        skeletal points and boundary points, and only those points.
        Calling this method before fit_ellipsoid_radii() will result in failure.

        Calling this method after it has already been called will not cause errors.
        Spokes will be straightened according to their skeletal and boundary points.

        Create spatial sample points along all the spokes.
        Allows spokes to "bend" during deformation, 
        and allows for the computation of onion skins and fitted frames.

        Before this method is called, the srep point indexing is such that
        if index i is a skeletal point, then i+1 is the boundary point 
        associated with that skeletal point (along the spoke).

        If s = the number of skeletal points, and i is the index of a skeletal point,
        then the indices of that skeletal points onion skin points will be
        2*s + i, 3*s + i, 4*s + i, ...
        """
        assert(self.polydata is not None)
        assert(num_onionskins >= 2) # need at least the skeleton and the boundary

        new_spokes_poly = vtk.vtkPolyData()
        new_spokes_cells = vtk.vtkCellArray()
        new_spokes_pts = vtk.vtkPoints()

        points = pv.wrap(self.polydata).points # old points, i.e. skeletal and boundary points
        num_skeletal_points = len(points) // 2

        # insert medial points, boundary points, and compute locations of onion skins
        all_onion_skin_points = []
        for i in range(num_skeletal_points):
            medial_point = points[2*i]
            boundary_point = points[2*i + 1]

            _ = new_spokes_pts.InsertNextPoint(medial_point)
            _ = new_spokes_pts.InsertNextPoint(boundary_point)

            onion_skin_points = np.linspace(medial_point, boundary_point, num=num_onionskins)
            onion_skin_points = onion_skin_points[1:-1]
            all_onion_skin_points.append(onion_skin_points)
        
        all_onion_skin_points = np.array(all_onion_skin_points)

        for i in range(num_onionskins - 2):
            for j in range(num_skeletal_points):
                _ = new_spokes_pts.InsertNextPoint(all_onion_skin_points[j,i,:])
        
        for i in range(num_skeletal_points):
            medial_pt_id = 2*i
            boundary_pt_id = 2*i + 1

            onionskin_pt_ids = [((2 + x) * num_skeletal_points + (medial_pt_id // 2)) for x in range(num_onionskins - 2)]

            spoke_point_ids = [ medial_pt_id]
            spoke_point_ids += onionskin_pt_ids
            spoke_point_ids.append(boundary_pt_id)

            points = [np.array(new_spokes_pts.GetPoint(pid)) for pid in spoke_point_ids]
            for j in range(len(spoke_point_ids) - 1):
                spoke_seg = vtk.vtkLine()
                spoke_seg.GetPointIds().SetId(0, spoke_point_ids[j])
                spoke_seg.GetPointIds().SetId(1, spoke_point_ids[j+1])
                new_spokes_cells.InsertNextCell(spoke_seg)
        
        new_spokes_poly.SetPoints(new_spokes_pts)
        new_spokes_poly.SetLines(new_spokes_cells)

        self.polydata = new_spokes_poly

        return


# # plots spine, interior skeleton, onion skins
# # TODO: incorporate this into the class

# s = len(my_srep.skeleton_point_indices)
# d = my_srep.num_spine_points

# spine_pts = my_srep.get_points_from_index_list(my_srep.skeleton_point_indices[:d])
# interior_skel_inds = []
# n = my_srep.num_spine_points
# m = my_srep.num_skeletal_onionskins
# for i in range(m-2):
#     interior_skel_inds += my_srep.skeleton_point_indices[
#         (1+2*(i))*n-(2*i):(1+2*(i+1))*n-(2*(i+1))
#     ]
# interior_skel_pts = my_srep.get_points_from_index_list(interior_skel_inds)
# start = (1+2*(m-2))*n - 2*(m-2)
# fold_inds = my_srep.skeleton_point_indices[-(2*n - 2):]
# fold_points = my_srep.get_points_from_index_list(fold_inds)

# plt = my_srep.plotter()
# plt.add_points(spine_pts, render_points_as_spheres=True, point_size=15, color="blue")
# plt.add_points(interior_skel_pts, render_points_as_spheres=True, point_size=15, color="red")
# plt.add_points(fold_points, render_points_as_spheres=True, point_size=15, color="cornflowerblue")
# plt.show()




