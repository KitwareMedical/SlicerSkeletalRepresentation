"""
Sometimes due to the rotational deformation, the linking correspondence will change from a pair
to another. This file helps to visualize such changes given the case id and new linked spoke id
"""

from viz import *
from compute_linking import *

data_folder = '../no_joint_data'
num_configs = 60
case_id = 4# np.random.randint(num_configs)
new_linked_spoke_id = 100

top_srep_file = os.path.join(data_folder, 's_reps/top_srep' + str(case_id)+'.vtk')
#top_srep_file = os.path.join(data_folder, 's_reps/top_ell_srep.vtk')
## target, neighbor s-reps and target spoke
reader = vtk.vtkPolyDataReader()
reader.SetFileName(top_srep_file)
reader.Update()
target_srep = reader.GetOutput()

bot_srep_file_name = os.path.join(data_folder, 's_reps/bot_srep' + str(case_id) + '.vtk')
#bot_srep_file_name = os.path.join(data_folder, 's_reps/bot_ell_srep.vtk')
bot_srep_reader = vtk.vtkPolyDataReader()
bot_srep_reader.SetFileName(bot_srep_file_name)
bot_srep_reader.Update()
neighbor_srep = bot_srep_reader.GetOutput()

top_mesh_file_name = os.path.join(data_folder, 'final_mesh/top' + str(case_id) + '_label_SPHARM.vtk')
mesh_reader = vtk.vtkPolyDataReader()
mesh_reader.SetFileName(top_mesh_file_name)
mesh_reader.Update()
target_mesh = mesh_reader.GetOutput()

bot_mesh_file_name = os.path.join(data_folder, 'final_mesh/bot' + str(case_id) + '_label_SPHARM.vtk')
bot_mesh_reader = vtk.vtkPolyDataReader()
bot_mesh_reader.SetFileName(bot_mesh_file_name)
bot_mesh_reader.Update()
bot_mesh = bot_mesh_reader.GetOutput()

spoke_id = 82
base_pt_id = spoke_id * 2
bdry_pt_id = base_pt_id + 1
base_pt = np.array(target_srep.GetPoint(base_pt_id))
bdry_pt = np.array(target_srep.GetPoint(bdry_pt_id))
s = bdry_pt - base_pt
r = np.linalg.norm(s)
u = s / r


def form_spoke_poly(np_base_pt, np_bdry_pt):
    target_spoke_poly = vtk.vtkPolyData()
    ts_pts = vtk.vtkPoints()
    vs_line = vtk.vtkCellArray()
    id0 = ts_pts.InsertNextPoint(np_base_pt)
    id1 = ts_pts.InsertNextPoint(np_bdry_pt)

    line = vtk.vtkLine()
    line.GetPointIds().SetId(0, id0)
    line.GetPointIds().SetId(1, id1)
    vs_line.InsertNextCell(line)

    target_spoke_poly.SetPoints(ts_pts)
    target_spoke_poly.SetLines(vs_line)
    target_spoke_poly.Modified()
    return target_spoke_poly

target_spoke_poly = form_spoke_poly(base_pt, bdry_pt)
append_filter = vtk.vtkAppendPolyData()
append_filter.AddInputData(target_srep)
append_filter.AddInputData(neighbor_srep)
append_filter.AddInputData(target_mesh)
append_filter.AddInputData(bot_mesh)
append_filter.Update()

def form_linked_spoke(spoke_id):
    base_pt_id = spoke_id * 2
    bdry_pt_id = base_pt_id + 1
    base_pt = np.array(neighbor_srep.GetPoint(base_pt_id))
    bdry_pt = np.array(neighbor_srep.GetPoint(bdry_pt_id))
    return form_spoke_poly(base_pt, bdry_pt)

linked_spokes_append = vtk.vtkAppendPolyData()
linked_spokes_append.AddInputData(target_spoke_poly)
linked_spokes_append.AddInputData(form_linked_spoke(82))
linked_spokes_append.Update()

linker = Linker()
target_spoke = Spoke(r, u, base_pt, None)
target_spoke, linked_spoke, linked_spoke_id = linker.link_spoke(target_spoke, neighbor_srep)

overlay_polydata(target_spoke_poly, append_filter.GetOutput(), form_linked_spoke(linked_spoke_id))
print(case_id, linked_spoke_id)

