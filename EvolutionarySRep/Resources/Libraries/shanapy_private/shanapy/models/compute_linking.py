import numpy as np
import os
from shanapy import models as sm
from shanapy import *
from shanapy.models.srep import SrepWrapper
from shanapy.utils import *
from numpy import linalg as LA
import vtk
from sklearn.kernel_approximation import RBFSampler
interp_target_level = 0

class Linker(object):
    def __init__(self, tgt=None, nbr=None, tgt_mesh=None, nbr_mesh=None):
        ### in the current implementation as of Jun 11, 2020
        ### there is only one target and one neighbor
        self.target = tgt       ## SrepWrapper that represents an srep
        self.neibor = nbr       ## SrepWrapper
        self.link_struct = None ## list of tuples of two-spokes
        self.tgt_mesh = tgt_mesh
        self.nbr_mesh = nbr_mesh
        self.srep_rows = 5
        self.srep_cols = 9

    def _compute_delta(self, target_spoke, neibor_spoke):
        """
        delta is the distance between z and z_prime (in Lemma 2.5)
        z is extension of a target boundary point outside the target
        z_prime is the extension of a neighbor boundary point
        The goal of finding optimal links is minimizing delta and then take the global minimal extension t
        """
        y_prime = neibor_spoke.getB()
        y = target_spoke.getB()
        # eq. (2.5)
        t = np.dot((y_prime - y), target_spoke.U) / (1 - np.dot(neibor_spoke.U, target_spoke.U))

        if t < 0: return np.inf, np.inf
        delta = LA.norm(y_prime - y + t * (neibor_spoke.U - target_spoke.U))
        return delta, t
    def export_links(self, output_file_name=None):
        """
        Save features to output_file_name

        Features are {(spoke1; spoke2), ..., ()}, where spoke_n is (p, U, r)
        """

        with open(output_file_name, 'a') as f:
            ## save linked region
            for link in self.link_struct:
                target_s, neibor_s = link
                extension = 0 if np.isinf(target_s.ext) else target_s.ext
                neibor_px, neibor_py, neibor_pz, neibor_ux, neibor_uy, neibor_uz, neibor_r = 0, 0, 0, 0, 0, 0, 0
                if neibor_s is not None:
                    neibor_px, neibor_py, neibor_pz, neibor_ux, neibor_uy, neibor_uz, neibor_r = neibor_s.p[0],neibor_s.p[1],neibor_s.p[2],neibor_s.U[0],neibor_s.U[1],neibor_s.U[2],neibor_s.r
                f.write('%f,%f,%f,%f,%f,%f,%f;%f,%f,%f,%f,%f,%f,%f\n' % (target_s.p[0],target_s.p[1],target_s.p[2],target_s.U[0],target_s.U[1],target_s.U[2],target_s.r, neibor_px, neibor_py, neibor_pz, neibor_ux, neibor_uy, neibor_uz, neibor_r) )
                ### aggregate features
                # joint_link = [target_s.p[0],target_s.p[1],target_s.p[2],target_s.U[0],target_s.U[1],target_s.U[2],target_s.r, neibor_px, neibor_py, neibor_pz, neibor_ux, neibor_uy, neibor_uz, neibor_r]
                # trans_X = rbf_feature.fit_transform(np.array(joint_link)[None, :]).squeeze()
                # f.write('%f,%f,%f,%f,%f,%f,%f\n' % (trans_X[0], trans_X[1], trans_X[2], trans_X[3], trans_X[4], trans_X[5], trans_X[6]) )
#                f.write('%f,%f,%f,%f,%f,%f,%f\n' % (target_s.p[0] + neibor_px,target_s.p[1]+neibor_py,target_s.p[2]+neibor_pz,target_s.U[0] + neibor_ux,target_s.U[1] + neibor_uy,target_s.U[2] + neibor_uz,target_s.r + neibor_r))

    def link_spoke(self, spoke: sm.Spoke, srep, verbose=False):
        """
        Find the link correspondence (a spoke) in srep to the spoke
        Input spoke: a single spoke in the target object
        Input srep: a list of models.Spokes
        Return the updated spoke (extension), and the link correspondence in srep
        """

        # delta_min = 10
        # delta_til = 10
        ret_link = None
        link_spoke_id = None
        extension_min = np.inf
        for i in range(len(srep)):
            base_pt_id = i * 2
            bdry_pt_id = i * 2 + 1
            neighbor_spoke = srep[i]

            delta, extension = self._compute_delta(spoke, neighbor_spoke)
            ext_threshold = 100
            # if verbose:
            #     print("Spokes meet up?", delta < spoke.delta_min)
            #     print("Global min ext?", extension < spoke.ext)
            #     print("Cutoff extension?", extension >= ext_threshold)

            if delta < spoke.delta_min and extension < spoke.ext and extension < ext_threshold:
                ### For debug visualization
                # appender = vtk.vtkAppendPolyData()
                # appender.AddInputData(self.tgt_mesh)
                # appender.AddInputData(self.nbr_mesh)

                # appender.AddInputData(self.target.data)
                # appender.AddInputData(self.neibor.data)
                # appender.Update()
                
                # combined_mesh = appender.GetOutput()

                # print('delta_min:' + str(spoke.delta_min) + ' t_min:' + str(spoke.ext)  + ' delta: ' + str(delta) + ' t:' + str(extension))
                # overlay_polydata(form_spoke_poly(spoke.p, spoke.extend(extension)), combined_mesh, form_spoke_poly(neighbor_spoke.p, neighbor_spoke.extend(extension)))
                ####
                extension_min = extension
                spoke.ext = np.dot((neighbor_spoke.getB() - spoke.getB()), spoke.U) / (1 - np.dot(neighbor_spoke.U, spoke.U))
                neighbor_spoke.ext = spoke.ext
                neighbor_spoke.ext_pt = spoke.getB() + spoke.ext * spoke.U
                ret_link = neighbor_spoke
                link_spoke_id = i
        return spoke, ret_link, link_spoke_id
    def link_multi_side(self, neighbor1=None, neighbor2=None, neighbor3=None, target1=None):
        """
        Input self.target is one side of s-reps (up/down/crest)
        Input neibor: try multiple side linking, neighbor1 and neighbor2 are the two sides (up/down), neighbor3 is the crest
        Input target1: the other side of the target for self-linking
        """
        ### 1. Interpolate target s-rep
        from shanapy import Interpolater
        interp = Interpolater()
        target_poly = self.target.data if isinstance(self.target, SrepWrapper) else self.target
        interpolated_tgt = interp.interpolate(target_poly, interp_target_level, self.target.rows, self.target.cols)

        prev_link = []
        final_link = None
        if neighbor1 is not None:
            link1 = self.link(interpolated_tgt, neighbor1)
            interpolated_tgt = [k[0] for k in link1]
            prev_link = [k[1] for k in link1]
        if neighbor2 is not None:
            link2 = self.link(interpolated_tgt, neighbor2, orig_link=prev_link)
            interpolated_tgt = [k[0] for k in link1]
            prev_link = [k[1] for k in link1]
        if neighbor3 is not None:
            link2 = self.link(interpolated_tgt, neighbor3, dunt_interp=True, orig_link=prev_link)
            interpolated_tgt = [k[0] for k in link1]
            prev_link = [k[1] for k in link1]

        if target1 is not None:
            final_link = self.link(interpolated_tgt, target1, orig_link=prev_link)
        ## link2 is the updated link structure
        return final_link
    def link(self, interpolated_tgt=None, neighbor=None, dunt_interp=False, orig_link=None):
        """
        An accurate linking structure that bases on
        interpolated s-reps of both the target and neibor
        """
        ### 1. Interpolate target s-rep
        from shanapy import Interpolater
        interp = Interpolater()
        if interpolated_tgt is None:
            target_poly = self.target.data if isinstance(self.target, SrepWrapper) else self.target
            interpolated_tgt = interp.interpolate(target_poly, interp_target_level, self.srep_rows, self.srep_cols)

        ### 2. Interpolate the neighbor object
        interpolated_nbr, nbr_srep = None, None
        if neighbor is not None and isinstance(neighbor, SrepWrapper):
            nbr_srep = neighbor.data
        else:
            nbr_srep = self.neibor.data

        if dunt_interp:
            interp_neibor_level = 0
        else:
            interp_neibor_level = 3
        interpolated_nbr = interp.interpolate(nbr_srep, interp_neibor_level, self.srep_rows, self.srep_cols)

        ### 3. Search neighbor for linking correspondence of target
        ret_link_struct = []
        for it, targetS in enumerate(interpolated_tgt):
#            if (it + 1) % 3 != 0: continue
#            print('Target spoke id: ' + str(it))
            extend_target_spoke, nbr_spoke, neighbor_id = self.link_spoke(targetS, interpolated_nbr)

            # if nbr_spoke is None:
            #     extend_target_spoke, nbr_spoke, neighbor_id = self.link_spoke(targetS, interpolated_nbr, verbose=True)
            if nbr_spoke is None and orig_link is not None:
                # preserve the previous link_struct
                ret_link_struct.append((extend_target_spoke, orig_link[it][1], it, orig_link[it][-1]))
            else:
                ret_link_struct.append((extend_target_spoke, nbr_spoke, it, neighbor_id))
        ### 4. Output the linking structure consists of the target s-rep and between-object pts
        self.link_struct = ret_link_struct
        return ret_link_struct

    def show_links(self, target_mesh=None, nbr_mesh=None):
        if len(self.link_struct) == 0:
            print('No links were found.')
            return
        ### now only visualize between-object sheet
        between_object_poly = vtk.vtkPolyData()
        between_object_pts = vtk.vtkPoints()
        between_object_cells = vtk.vtkCellArray()

        no_link_spokes = 0
        all_links = 0
        spokes_append = vtk.vtkAppendPolyData()
        link_append = vtk.vtkAppendPolyData()
        for link in self.link_struct:
            spoke, neighbor_spoke, _, _ = link
            if abs(spoke.ext) < 1e-4 or neighbor_spoke is None:
                no_link_spokes += 1

                continue
            pt = spoke.p + spoke.U * (spoke.r + spoke.ext)
            id_pt = between_object_pts.InsertNextPoint(pt)
            between_object_cells.InsertNextCell(1)
            between_object_cells.InsertCellPoint(id_pt)
            all_links += 1

            spoke_target, link_target = spoke.visualize()
            spoke_neighbor, link_neighbor = neighbor_spoke.visualize()
            spokes_append.AddInputData(spoke_target)
            spokes_append.AddInputData(spoke_neighbor)
            link_append.AddInputData(link_target)
            link_append.AddInputData(link_neighbor)

        print('Found ' + str(all_links) + ' links.')
        print(str(no_link_spokes) + ' spokes in the target object have no links.')
        spokes_append.Update()
        link_append.Update()
        between_object_poly.SetPoints(between_object_pts)
        between_object_poly.SetVerts(between_object_cells)

        # show original s-reps
        target_poly = self.target.data if isinstance(self.target, SrepWrapper) else self.target
        neibor_poly = self.neibor.data if isinstance(self.neibor, SrepWrapper) else self.neibor
        appender = vtk.vtkAppendPolyData()
        appender.AddInputData(target_poly)
        appender.AddInputData(neibor_poly)
        appender.Update()

        combined_srep = appender.GetOutput()

        # fit a surface to points
        # splatter = vtk.vtkGaussianSplatter()
        # splatter.SetInputData(between_object_poly)
        # splatter.SetSampleDimensions(50,50,50)
        # splatter.SetRadius(0.5)
        # splatter.ScalarWarpingOff()

        # surface = vtk.vtkContourFilter()
        # surface.SetInputConnection(splatter.GetOutputPort())
        # surface.SetValue(0,0.01)
        if target_mesh and nbr_mesh:
            mesh_append = vtk.vtkAppendPolyData()
            mesh_append.AddInputData(target_mesh)
            mesh_append.AddInputData(nbr_mesh)
            mesh_append.Update()

#            visualize(mesh_append.GetOutput())
            overlay_polydata(between_object_poly, mesh_append.GetOutput(), spokes_append.GetOutput(), link_append.GetOutput())#, combined_srep)
        else:
            overlay_polydata(between_object_poly, highlight=spokes_append.GetOutput())
#        return between_object_poly # surface.GetOutput()
