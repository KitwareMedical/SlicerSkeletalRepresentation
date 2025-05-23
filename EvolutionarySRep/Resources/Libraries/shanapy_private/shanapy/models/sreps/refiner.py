"""
The refinement aims to optimize the interior geometry represented by the s-rep.
This optimization ought to consider i) the goodness of fitting to the boundary geometry and 
ii) the smoothness of interior radial distance level surfaces.
As of Dec. 27, 2021, the refinement in this python package can optimize the spokes' lengths.

Edited by Nick Tapp-Hughes Jan 24 2022.
Removed optimization of spoke length because all I need is fold curve stuff.
"""
from audioop import avg
import numpy as np
import vtk
import pyvista as pv
class Refiner:
    """This class optimize an initial s-rep to better fit to the boundary"""
    def __init__(self):
        ## TODO: Set parameters of refinement
        self.eps = np.finfo(float).eps
        pass
    def relocate(self, bdry_pt, input_mesh):
        """
        Relocate base points of a (fold) spoke with the two end points: base_pt and bdry_pt,
        such that the length of the spoke is reciprocal of boundary curvature
        """
        # find the closest mesh point to the tip of the spoke
        cell_locator = vtk.vtkCellLocator()
        cell_locator.SetDataSet(input_mesh)
        cell_locator.BuildLocator()

        cellId = vtk.reference(0)
        c = [0.0, 0.0, 0.0]
        subId = vtk.reference(0)
        d = vtk.reference(0.0)
        cell_locator.FindClosestPoint(bdry_pt, c, cellId, subId, d)
        pt_ids = vtk.vtkIdList()
        input_mesh.GetCellPoints(cellId, pt_ids)

        curvature = vtk.vtkCurvatures()
        curvature.SetInputData(input_mesh)
        curvature.SetCurvatureTypeToMean()
        curvature.Update()
        mean_curvatures = curvature.GetOutput()
        mean_curvature = []
        for i in range(pt_ids.GetNumberOfIds()):
            pt_id = pt_ids.GetId(i)
            mean_curvature.append(mean_curvatures.GetPointData().GetArray(0).GetValue(pt_id))
            
        # print(1/np.mean(mean_curvature))
        # quit()
        return 1/ np.mean(mean_curvature) / 10
    
    def refine(self, srep, input_mesh, num_crest_points=24):
        """
        The main entry of the refinement
        Input: an initial s-rep  srep
        Input: the boundary mesh input_mesh
        Return spokes poly (a set of spokes that can be visualized) and fold points
        """
        print('Refining ...')
        srep_poly = vtk.vtkPolyData()
        srep_poly.DeepCopy(srep)
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
        
        # p = pv.Plotter()
        # p.add_mesh(input_mesh, color='white', opacity=0.2)
        # p.add_mesh(srep)
        # p.add_points(base_array, render_points_as_spheres=True, point_size=30.)
        # p.show()
        # quit()

        # print(radii_array)
        # quit()

        ## update radii of s-rep and return the updated
        arr_length = vtk.vtkDoubleArray()
        arr_length.SetNumberOfComponents(1)
        arr_length.SetName("spokeLength")

        arr_dirs = vtk.vtkDoubleArray()
        arr_dirs.SetNumberOfComponents(3)
        arr_dirs.SetName("spokeDirection")

        bdry_pts = []

        for i in range(num_spokes):
            id_base_pt = i * 2
            id_bdry_pt = id_base_pt + 1
            base_pt = base_array[i, :]
            radius = radii_array[i]
            direction = dir_array[i, :]

            new_bdry_pt = base_pt + radius * direction
            bdry_pts.append(new_bdry_pt)

            arr_length.InsertNextValue(radius)
            arr_dirs.InsertNextTuple(direction)
            srep_poly.GetPoints().SetPoint(id_bdry_pt, new_bdry_pt)

            ### relocate base points for fold spokes such that their lengths are reciprocal of boundary mean curvature 
            if i >= num_spokes - num_crest_points:
                # new_radius = min(radius - 1, self.relocate(new_bdry_pt, input_mesh))
                new_radius = self.relocate(new_bdry_pt, input_mesh)
                
                # new_base_pt = new_bdry_pt + new_radius * direction
                new_base_pt = new_bdry_pt - new_radius * direction
                srep_poly.GetPoints().SetPoint(id_base_pt, new_base_pt)
        

        # bdry_pts = np.array(bdry_pts)
        # p = pv.Plotter()
        # p.add_mesh(input_mesh, color='white', opacity=0.2)
        # p.add_mesh(srep)
        # p.add_points(bdry_pts, render_points_as_spheres=True, point_size=30.)
        # p.show()
        # quit()

        srep_poly.GetPointData().AddArray(arr_length)
        srep_poly.GetPointData().AddArray(arr_dirs)
        srep_poly.Modified()
        return srep_poly