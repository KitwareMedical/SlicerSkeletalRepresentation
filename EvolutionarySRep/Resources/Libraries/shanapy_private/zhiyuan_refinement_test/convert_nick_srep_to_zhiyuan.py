import vtk
import numpy as np
import pyvista as pv
import pickle

# tau_1 and tau_2 are on skeleton: tau_1 refers to the first level curve surrounding the spine;
#  tau_2 refers to the next circle out

srep_reader = vtk.vtkPolyDataReader()
srep_reader.SetFileName('srep_before_refinement_after_regularization.vtk')
srep_reader.Update()
srep_poly = srep_reader.GetOutput()

def reorderOnOneSide(offset=0):
    """
    
    """
    # treat the spine as a degenerated circle, allowing redundant skeletal points on spine
    spine_base = []
    spine_base_inds = []
    spine_bdry = []
    spine_bdry_inds = []
    for i in range(13):
        base_id = i + offset
        base_pt = np.array(srep_poly.GetPoint(base_id * 2))
        bdry_pt = np.array(srep_poly.GetPoint(base_id * 2 + 1))
        spine_base.append(base_pt)
        spine_base_inds.append(base_id * 2)
        spine_bdry.append(bdry_pt)
        spine_bdry_inds.append(base_id * 2 + 1)

    spine_base = np.array(spine_base + spine_base[1:-1][::-1])
    spine_base_inds = spine_base_inds + spine_base_inds[1:-1][::-1]
    spine_bdry = np.array(spine_bdry + spine_bdry[1:-1][::-1]) # 24 x 3 
    spine_bdry_inds = spine_bdry_inds + spine_bdry_inds[1:-1][::-1]

    # now tau_1 = first circle of points around the spine on the skeleton

    tau_1_base = []
    tau_1_base_inds = []
    tau_1_bdry = []
    tau_1_bdry_inds = []
    for i in range(13):
        base_id = i + 13 + offset
        base_pt = np.array(srep_poly.GetPoint(base_id * 2))
        bdry_pt = np.array(srep_poly.GetPoint(base_id * 2 + 1))
        tau_1_base.append(base_pt)
        tau_1_base_inds.append(base_id * 2)
        tau_1_bdry.append(bdry_pt)
        tau_1_bdry_inds.append(base_id * 2 + 1)

    tau_1_base = np.array(tau_1_base)
    tau_1_bdry = np.array(tau_1_bdry) # 13 x 3

    tau_1_base_2 = []
    tau_1_base_2_inds = []
    tau_1_bdry_2 = []
    tau_1_bdry_2_inds = []
    for i in range(11):
        base_id = i + 26 + offset
        base_pt = np.array(srep_poly.GetPoint(base_id * 2))
        bdry_pt = np.array(srep_poly.GetPoint(base_id * 2 + 1))
        tau_1_base_2.append(base_pt)
        tau_1_base_2_inds.append(base_id * 2)
        tau_1_bdry_2.append(bdry_pt)
        tau_1_bdry_2_inds.append(base_id * 2 + 1)
    tau_1_base_2 = np.array(tau_1_base_2)
    tau_1_bdry_2 = np.array(tau_1_bdry_2) # 11 x 3

    tau_1_base_all = np.concatenate((tau_1_base[0, :][None, :], tau_1_base_2, tau_1_base[-1:0:-1, :]))
    tau_1_base_all_inds = [tau_1_base_inds[0]] + tau_1_base_2_inds + tau_1_base_inds[-1:0:-1]
    tau_1_bdry_all = np.concatenate((tau_1_bdry[0, :][None, :], tau_1_bdry_2, tau_1_bdry[-1:0:-1, :])) # 24 x 3
    tau_1_bdry_all_inds = [tau_1_bdry_inds[0]] + tau_1_bdry_2_inds + tau_1_bdry_inds[-1:0:-1]

    # now tau_2 = next circle out on the skeleton

    tau_2_base = []
    tau_2_base_inds = []
    tau_2_bdry = []
    tau_2_bdry_inds = []
    for i in range(13):
        base_id = i + 24 + 13 + offset
        base_pt = np.array(srep_poly.GetPoint(base_id * 2))
        bdry_pt = np.array(srep_poly.GetPoint(base_id * 2 + 1))
        tau_2_base.append(base_pt)
        tau_2_base_inds.append(base_id * 2)
        tau_2_bdry.append(bdry_pt)
        tau_2_bdry_inds.append(base_id * 2 + 1)
    tau_2_base = np.array(tau_2_base)
    tau_2_bdry = np.array(tau_2_bdry) # 13 x 3

    tau_2_base_2 = []
    tau_2_base_2_inds = []
    tau_2_bdry_2 = []
    tau_2_bdry_2_inds = []
    for i in range(11):
        base_id = i + 50 + offset
        base_pt = np.array(srep_poly.GetPoint(base_id * 2))
        bdry_pt = np.array(srep_poly.GetPoint(base_id * 2 + 1))
        tau_2_base_2.append(base_pt)
        tau_2_base_2_inds.append(base_id * 2)
        tau_2_bdry_2.append(bdry_pt)
        tau_2_bdry_2_inds.append(base_id * 2 + 1)
    tau_2_base_2 = np.array(tau_2_base_2)
    tau_2_bdry_2 = np.array(tau_2_bdry_2) # 11 x 3

    tau_2_base_all = np.concatenate((tau_2_base[0, :][None, :], tau_2_base_2, tau_2_base[-1:0:-1, :]))
    tau_2_base_all_inds = [tau_2_base_inds[0]] + tau_2_base_2_inds + tau_2_base_inds[-1:0:-1]
    tau_2_bdry_all = np.concatenate((tau_2_bdry[0, :][None, :], tau_2_bdry_2, tau_2_bdry[-1:0:-1, :])) # 24 x 3
    tau_2_bdry_all_inds = [tau_2_bdry_inds[0]] + tau_2_bdry_2_inds + tau_2_bdry_inds[-1:0:-1]

    # now putting spine, tau_1, tau_2 all together

    skeletal_pts_all = []
    skeletal_pts_all_inds = []
    bdry_pts_all = []
    bdry_pts_all_inds = []
    for i in range(24):
        skeletal_pts_all.append(spine_base[i, :])
        skeletal_pts_all_inds.append(spine_base_inds[i])
        skeletal_pts_all.append(tau_1_base_all[i, :])
        skeletal_pts_all_inds.append(tau_1_base_all_inds[i])
        skeletal_pts_all.append(tau_2_base_all[i, :])
        skeletal_pts_all_inds.append(tau_2_base_all_inds[i])

        bdry_pts_all.append(spine_bdry[i, :])
        bdry_pts_all_inds.append(spine_bdry_inds[i])
        bdry_pts_all.append(tau_1_bdry_all[i, :])
        bdry_pts_all_inds.append(tau_1_bdry_all_inds[i])
        bdry_pts_all.append(tau_2_bdry_all[i, :])
        bdry_pts_all_inds.append(tau_2_bdry_all_inds[i])

    skeletal_pts_all = np.array(skeletal_pts_all)
    bdry_pts_all = np.array(bdry_pts_all)

    return skeletal_pts_all, bdry_pts_all, skeletal_pts_all_inds, bdry_pts_all_inds

top_skeletal_pts, top_bdry_pts, top_skeletal_pts_inds, top_bdry_pts_inds = reorderOnOneSide()
bot_skeletal_pts, bot_bdry_pts, bot_skeletal_pts_inds, bot_bdry_pts_inds = reorderOnOneSide(13+24+24)

## fold spokes
fold_base_pts = []
fold_base_pts_inds = []
fold_bdry_pts = []
fold_bdry_pts_inds = []
for i in range(24):
    base_id = i + 2*(13+24+24)
    base_pt = np.array(srep_poly.GetPoint(base_id * 2))
    bdry_pt = np.array(srep_poly.GetPoint(base_id * 2+1))
    fold_base_pts.append(base_pt)
    fold_base_pts_inds.append(base_id * 2)
    fold_bdry_pts.append(bdry_pt)
    fold_bdry_pts_inds.append(base_id * 2 + 1)

fold_base_pts = np.array([fold_base_pts[0]] + fold_base_pts[13:] + fold_base_pts[1:13][::-1])
fold_base_pts_inds = [fold_base_pts_inds[0]] + fold_base_pts_inds[13:] + fold_base_pts_inds[1:13][::-1]
fold_bdry_pts = np.array([fold_bdry_pts[0]] + fold_bdry_pts[13:] + fold_bdry_pts[1:13][::-1])
fold_bdry_pts_inds = [fold_bdry_pts_inds[0]] + fold_bdry_pts_inds[13:] + fold_bdry_pts_inds[1:13][::-1]

all_skeletal_pts = np.concatenate((top_skeletal_pts, bot_skeletal_pts, fold_base_pts))
all_skeletal_pts_inds = top_skeletal_pts_inds + bot_skeletal_pts_inds + fold_base_pts_inds
all_bdry_pts = np.concatenate((top_bdry_pts, bot_bdry_pts, fold_bdry_pts))
all_bdry_pts_inds = top_bdry_pts_inds + bot_bdry_pts_inds + fold_bdry_pts_inds

# Now create new re-ordered srep with new information

srep = vtk.vtkPolyData()
srep_pt = vtk.vtkPoints()
spoke_cells = vtk.vtkCellArray()
for i in range(all_skeletal_pts.shape[0]):
    id1 = srep_pt.InsertNextPoint(all_skeletal_pts[i, :])
    id2 = srep_pt.InsertNextPoint(all_bdry_pts[i, :])

    spoke = vtk.vtkLine()
    spoke.GetPointIds().SetId(0, id1)
    spoke.GetPointIds().SetId(1, id2)
    spoke_cells.InsertNextCell(spoke)
srep.SetPoints(srep_pt)
srep.SetLines(spoke_cells)

srep_writer = vtk.vtkPolyDataWriter()
srep_writer.SetInputData(srep)
srep_writer.SetFileName('reordered_srep.vtk')
srep_writer.Update()

# Save index mapping lists
"""
Currently, the all_skeletal_pts_inds list has the property that

n_srep.GetPoint(all_skeletal_pts_inds[z_srep_ind // 2])
==
z_srep.GetPoint(z_srep_ind)

i.e. the elements of the list are n_srep indices, and the indices of the list are the z_srep indices.

But we want to convert it so that we have a new list 'mapping' with

z_srep.GetPoint(mapping[n_srep_ind // 2])
==
n_srep.GetPoint(n_srep_ind)

So I do that here:
"""
assert len(all_skeletal_pts_inds) == len(all_bdry_pts_inds)
assert len(all_skeletal_pts_inds) == (srep.GetNumberOfPoints() // 2)

mapping = -np.ones((len(all_skeletal_pts_inds) * 2)).astype(int)
for i, oinds in enumerate(zip(all_skeletal_pts_inds, all_bdry_pts_inds)):
    rind_skel = 2*i
    rind_bdry = 2*i + 1
    oind_skel = oinds[0]
    oind_bdry = oinds[1]

    try:
        d = np.linalg.norm(np.array(srep_poly.GetPoint(oind_skel)) - np.array(srep.GetPoint(rind_skel)))
        assert d < 1e-3
    except:
        print(rind_skel)
        print(oind_skel)
        print(d)
        quit()

    try:
        d = np.linalg.norm(np.array(srep_poly.GetPoint(oind_bdry)) - np.array(srep.GetPoint(rind_bdry)))
        assert d < 1e-3
    except:
        print(rind_bdry)
        print(oind_bdry)
        print(d)
        quit()

    mapping[oind_skel] = rind_skel
    mapping[oind_bdry] = rind_bdry

mapping = list(mapping)

# get rid of any extra -1s since the num points in srep is greater than that in srep_poly
mapping = [x for x in mapping if x != -1]

with open('nick_to_zhiyuan_srep_index_mapping.pkl', 'wb') as f:
    pickle.dump(mapping, f)

print('done')
