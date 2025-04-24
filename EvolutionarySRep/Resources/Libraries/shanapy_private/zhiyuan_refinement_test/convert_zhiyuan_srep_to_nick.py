"""Inverts the srep re-ordering so that we get back to Nick's original ordering."""
import vtk
import pickle
import numpy as np
import pyvista as pv

# Load the list from the file using pickle
with open('nick_to_zhiyuan_srep_index_mapping.pkl', 'rb') as f:
    mapping = pickle.load(f)

reordered_srep = pv.PolyData("reordered_srep.vtk")

# Create new srep with original ordering

srep = vtk.vtkPolyData()
srep_pt = vtk.vtkPoints()
spoke_cells = vtk.vtkCellArray()

for ind_skel, ind_bdry in zip(mapping[::2], mapping[1::2]):
    id1 = srep_pt.InsertNextPoint(reordered_srep.GetPoint(ind_skel))
    id2 = srep_pt.InsertNextPoint(reordered_srep.GetPoint(ind_bdry))

    spoke = vtk.vtkLine()
    spoke.GetPointIds().SetId(0, id1)
    spoke.GetPointIds().SetId(1, id2)
    spoke_cells.InsertNextCell(spoke)
srep.SetPoints(srep_pt)
srep.SetLines(spoke_cells)

srep_writer = vtk.vtkPolyDataWriter()
srep_writer.SetInputData(srep)
srep_writer.SetFileName('srep_ordered_back.vtk')
srep_writer.Update()