import vtk
import math
import numpy as np
import sys
import os
# from curvedSrep import curvedSrep as cs
sys.path.append(os.path.abspath("."))
import matplotlib.pyplot as plt
from numpy import random
from scipy.spatial import distance
from shanapy.models import mean_curvature_flow as mcf, curvedSrep as cs, srep_fitter
from shanapy.utils import viz
import pyvista as pv
# method to find closest points on 2D skeleton
def form_frames(spokeT, finU, finT, finN, fitFrames_ends, fitFrames_lines):

    id0 = fitFrames_ends.InsertNextPoint(tuple(spokeT))
    id1 = fitFrames_ends.InsertNextPoint(tuple(finU))
    spoke_line = vtk.vtkLine()
    spoke_line.GetPointIds().SetId(0, id0)
    spoke_line.GetPointIds().SetId(1, id1)
    fitFrames_lines.InsertNextCell(spoke_line)

    id0 = fitFrames_ends.InsertNextPoint(tuple(spokeT))
    id1 = fitFrames_ends.InsertNextPoint(tuple(finT))
    spoke_line = vtk.vtkLine()
    spoke_line.GetPointIds().SetId(0, id0)
    spoke_line.GetPointIds().SetId(1, id1)
    fitFrames_lines.InsertNextCell(spoke_line)

    id0 = fitFrames_ends.InsertNextPoint(tuple(spokeT))
    id1 = fitFrames_ends.InsertNextPoint(tuple(finN))
    spoke_line = vtk.vtkLine()
    spoke_line.GetPointIds().SetId(0, id0)
    spoke_line.GetPointIds().SetId(1, id1)
    fitFrames_lines.InsertNextCell(spoke_line)
    return fitFrames_ends, fitFrames_lines

def findClosetPt(pt, pts, sampling):
    # test = pt[2]
    dist = (pt[0]-pts[0][0])**2+(pt[1]-pts[0][1])**2+(pt[2]-pts[0][2])**2
    ind0 = 0
    # print(pts)
    for i in range(0, len(pts)):
        val = (pt[0]-pts[i][0])**2+(pt[1]-pts[i][1])**2+(pt[2]-pts[i][2])**2
        if val < dist:
            dist = val
            ind0 = i
    
    nextPt = (ind0+ (sampling*2)) % len(pts)
    prevPt = ind0-(sampling*2)
    if ind0 -  sampling < 0:
        prevPt = ind0+sampling
    elif ind0 - (sampling*2) < 0:
        prevPt = ind0-sampling
    if ind0 + sampling > len(pts):
        nextPt = ind0 - (sampling*2)
        prevPt = ind0 - sampling
    elif ind0 + (sampling*2) > len(pts):
        nextPt = ind0 + sampling
    distComp1 = (pt[0]-pts[nextPt][0])**2+(pt[1]-pts[nextPt][1])**2+(pt[2]-pts[nextPt][2])**2
    distComp2 = (pt[0]-pts[prevPt][0])**2+(pt[1]-pts[prevPt][1])**2+(pt[2]-pts[prevPt][2])**2

    if distComp1 < distComp2:
        ind1 = nextPt
    else:
        ind1 = prevPt

    if pts[ind0][0] == pts[ind1][0] and pts[ind0][1] == pts[ind1][1] and pts[ind0][2] == pts[ind1][2]:
        if ind1 == nextPt:
            ind1 = prevPt
        else:
            ind1 = nextPt

    # print([ind0, ind1])
    return dist, [ind0, ind1]


def fit_local_frames(mesh_file_path = '/playpen/workspace/Simulate_Shapes/data/107524.vtk', output='../../data/frames.vtk'):
    
    #### Here is an example usage of mcf for srep fitting -- Nick
    # mesh_file_path = '/playpen/workspace/Simulate_Shapes/data/ell_107524.vtk'
    #### fit ellipsoid srep 
    # bot_mesh, tps_list = mcf.flow(mesh_file_path)

    # eigen_vectors, init_bot_srep = mcf.fit_srep_to_quasi_ellipsoid(bot_mesh, auto_flip=False)
    # bot_srep = srep_fitter.refine_srep(init_bot_srep, bot_mesh)

    # ell_mesh_writer = vtk.vtkPolyDataWriter()
    # ell_mesh_writer.SetFileName('/playpen/workspace/Simulate_Shapes/data/ell_107524.vtk')
    # ell_mesh_writer.SetInputData(bot_mesh)
    # ell_mesh_writer.Update()

    # new_srep_ell = vtk.vtkPolyData()
    # spokes_poly = vtk.vtkCellArray()
    # new_srep_pts = vtk.vtkPoints()
    # for i in range(0, bot_srep.GetNumberOfPoints(), 2):
    #     new_srep_pts.InsertNextPoint(bot_srep.GetPoint(i))
    #     new_srep_pts.InsertNextPoint(bot_srep.GetPoint(i+1))
    #     spoke_line = vtk.vtkLine()
    #     spoke_line.GetPointIds().SetId(0, i)
    #     spoke_line.GetPointIds().SetId(1, i + 1)
    #     spokes_poly.InsertNextCell(spoke_line)
    # new_srep_ell.SetLines(spokes_poly)
    # new_srep_ell.SetPoints(new_srep_pts)
    # new_srep_ell.Modified()
    # ell_srep_writer = vtk.vtkPolyDataWriter()
    # ell_srep_writer.SetFileName('/playpen/workspace/Simulate_Shapes/data/ell_srep_107524.vtk')
    # ell_srep_writer.SetInputData(new_srep_ell)
    # ell_srep_writer.Update()

    # bot_srep = new_srep_ell
    bot_mesh_reader = vtk.vtkPolyDataReader()
    bot_mesh_reader.SetFileName('/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/107524/ellipsoid.vtk')
    bot_mesh_reader.Update()
    bot_mesh = bot_mesh_reader.GetOutput()

    bot_srep_reader = vtk.vtkPolyDataReader()
    bot_srep_reader.SetFileName('/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/107524/ell_srep.vtk')
    bot_srep_reader.Update()
    bot_srep = bot_srep_reader.GetOutput()
    # p = pv.Plotter()
    # p.add_mesh(bot_mesh, opacity=0.2)
    # p.add_mesh(bot_srep)

    # p.show()

    # get all the skeletal and boundary points of the ellipsoid s-rep
    # ax = plt.subplot(111, projection='3d')
    source_pts = vtk.vtkPoints()
    for i in range(bot_srep.GetNumberOfCells()):
        base_pt_id = i * 2
        bdry_pt_id = i * 2 + 1
        s_pt = bot_srep.GetPoint(base_pt_id)
        b_pt = bot_srep.GetPoint(bdry_pt_id)
        # # ax.scatter(s_pt[0], s_pt[1], s_pt[2], color='k')
        source_pts.InsertNextPoint([s_pt[0], s_pt[1], s_pt[2]])
    source_pts.Modified()

    # getting 2D s-rep of the ellipsoid skeleton
    rXs, rYs, rZs, rSamPts = cs.curvedSrep(source_pts)
    numSamp = len(rSamPts[0])
    sampNum = numSamp
    print(numSamp)

    iXs = []
    iYs = []
    iZs = []
    iSamPts = []
    srepPts = vtk.vtkPoints()
    pPts = vtk.vtkPoints()

    # viewing the 2D skeleton
    ptsOnSkel = []
    for i in range(0, len(rXs)):
        intPt = (rXs[i], rYs[i], rZs[i])
        #ax.scatter(rXs[i], rYs[i], rZs[i], color='b')
        iXs.append(intPt[0])
        iYs.append(intPt[1])
        iZs.append(intPt[2])

        if i == 0:
            iSpoke = []
            srepPts.InsertNextPoint([iXs[i], iYs[i], iZs[i]])
            ptsOnSkel.append([iXs[i], iYs[i], iZs[i]])
            for j in range(0, numSamp):
                intPt = rSamPts[i][j]
                iSpoke.append([intPt[0], intPt[1], intPt[2]])
                srepPts.InsertNextPoint([intPt[0], intPt[1], intPt[2]])
                ptsOnSkel.append([intPt[0], intPt[1], intPt[2]])
                # ax.scatter(intPt[0], intPt[1], intPt[2], color='b')
                # if j == 0:
                #     # ax.plot([iXs[i], iSpoke[0][0]], [iYs[i], iSpoke[0][1]], [iZs[i], iSpoke[0][2]], 'r')
                # else:
                #     # ax.plot([iSpoke[j-1][0], iSpoke[j][0]], [iSpoke[j-1][1], iSpoke[j][1]], [iSpoke[j-1][2], iSpoke[j][2]], 'r')
            iSamPts.append(iSpoke)
        elif i == len(rXs)-1:
            iSpoke = []
            srepPts.InsertNextPoint([iXs[i], iYs[i], iZs[i]])
            ptsOnSkel.append([iXs[i], iYs[i], iZs[i]])
            for j in range(0, numSamp):
                intPt = rSamPts[-1][j]
                iSpoke.append([intPt[0], intPt[1], intPt[2]])
                srepPts.InsertNextPoint([intPt[0], intPt[1], intPt[2]])
                ptsOnSkel.append([intPt[0], intPt[1], intPt[2]])
                # ax.scatter(intPt[0], intPt[1], intPt[2], color='b')
                # if j == 0:
                #     # ax.plot([iXs[i], iSpoke[0][0]], [iYs[i], iSpoke[0][1]], [iZs[i], iSpoke[0][2]], 'r')
                # else:
                #     # ax.plot([iSpoke[j-1][0], iSpoke[j][0]], [iSpoke[j-1][1], iSpoke[j][1]], [iSpoke[j-1][2], iSpoke[j][2]], 'r')
            iSamPts.append(iSpoke)
        else:
            iSpoke = []
            srepPts.InsertNextPoint([iXs[i], iYs[i], iZs[i]])
            ptsOnSkel.append([iXs[i], iYs[i], iZs[i]])
            for j in range(0, numSamp):
                intPt = rSamPts[2*i-1][j]
                iSpoke.append([intPt[0], intPt[1], intPt[2]])
                srepPts.InsertNextPoint([intPt[0], intPt[1], intPt[2]])
                ptsOnSkel.append([intPt[0], intPt[1], intPt[2]])
                # ax.scatter(intPt[0], intPt[1], intPt[2], color='b')
                # if j == 0:
                #     # ax.plot([iXs[i], iSpoke[0][0]], [iYs[i], iSpoke[0][1]], [iZs[i], iSpoke[0][2]], 'r')
                # else:
                #     # ax.plot([iSpoke[j-1][0], iSpoke[j][0]], [iSpoke[j-1][1], iSpoke[j][1]], [iSpoke[j-1][2], iSpoke[j][2]], 'r')
            iSamPts.append(iSpoke)

            iSpoke = []
            srepPts.InsertNextPoint([iXs[i], iYs[i], iZs[i]])
            ptsOnSkel.append([iXs[i], iYs[i], iZs[i]])
            for j in range(0, numSamp):
                intPt = rSamPts[2*i][j]
                iSpoke.append([intPt[0], intPt[1], intPt[2]])
                srepPts.InsertNextPoint([intPt[0], intPt[1], intPt[2]])
                ptsOnSkel.append([intPt[0], intPt[1], intPt[2]])
                # ax.scatter(intPt[0], intPt[1], intPt[2], color='b')
                # if j == 0:
                #     # ax.plot([iXs[i], iSpoke[0][0]], [iYs[i], iSpoke[0][1]], [iZs[i], iSpoke[0][2]], 'r')
                # else:
                #     # ax.plot([iSpoke[j-1][0], iSpoke[j][0]], [iSpoke[j-1][1], iSpoke[j][1]], [iSpoke[j-1][2], iSpoke[j][2]], 'r')
            iSamPts.append(iSpoke)
    srepSet = vtk.vtkPolyData()
    srepSet.SetPoints(srepPts)
    srepSet.Modified()
    # ax.plot(iXs, iYs, iZs, 'r')

    print(len(rXs))
    print(len(ptsOnSkel))

    # using thin plate splines to be able to interpolate 3D spokes on the ellipsoid s-rep
    source_pts2 = vtk.vtkPoints()
    target_pts2 = vtk.vtkPoints()
    for i in range(math.floor((bot_srep.GetNumberOfCells() -24) / 2)):
        base_pt_id = i * 2
        bdry_pt_id = i * 2 + 1
        s_pt = bot_srep.GetPoint(base_pt_id)
        b_pt = bot_srep.GetPoint(bdry_pt_id)
        source_pts2.InsertNextPoint(s_pt)
        target_pts2.InsertNextPoint(b_pt)
    tps = vtk.vtkThinPlateSplineTransform()
    tps.SetSourceLandmarks(source_pts2)
    tps.SetTargetLandmarks(target_pts2)
    tps.SetBasisToR()
    tps.Modified()

    # ellipsoid mesh
    # reader3 = vtk.vtkPolyDataReader()
    # reader3.SetFileName('sreps/data/recenter.vtk')
    # reader3.Update()
    # bot_mesh = reader3.GetOutput()
    # mesh of target object
    reader4 = vtk.vtkPolyDataReader()
    # reader4.SetFileName('sreps/control/FinTopMesh1.vtk')
    reader4.SetFileName(mesh_file_path)
    reader4.Update()
    top_mesh = reader4.GetOutput()

    source_pts3 = vtk.vtkPoints()
    target_pts3 = vtk.vtkPoints()

    ellPts = bot_mesh.GetNumberOfPoints()

    # getting points on ellipsoid and target object mesh
    for i in range(0, ellPts):
        pt = [0] * 3
        bot_mesh.GetPoint(i, pt)
        # print(pt)
        source_pts3.InsertNextPoint(pt)

        top_mesh.GetPoint(i, pt)
        # print(pt)
        target_pts3.InsertNextPoint(pt)
        # if i == 5:
        #     break
    source_pts3.Modified()
    target_pts3.Modified()

    ### Interpolate deformation with thin-plate-spline
    tps2 = vtk.vtkThinPlateSplineTransform()
    tps2.SetSourceLandmarks(source_pts3)
    tps2.SetTargetLandmarks(target_pts3)
    tps2.SetBasisToR()
    tps2.Modified()

    fitFrames = vtk.vtkPolyData()
    fitFrames_ends = vtk.vtkPoints()
    fitFrames_lines = vtk.vtkCellArray()
    special_fitFrames = vtk.vtkPolyData()
    special_fitFrames_ends = vtk.vtkPoints()
    special_fitFrames_lines = vtk.vtkCellArray()

    highlight_fitFrames = vtk.vtkPolyData()
    highlight_fitFrames_ends = vtk.vtkPoints()
    highlight_fitFrames_lines = vtk.vtkCellArray()
    highlight_special_fitFrames = vtk.vtkPolyData()
    highlight_special_fitFrames_ends = vtk.vtkPoints()
    highlight_special_fitFrames_lines = vtk.vtkCellArray()

#    special_ids = [0, 6*3, 6*3 + 2, 12 * 3, 18*3 + 2]
    special_ids = [0, 6*3, 12 * 3] # 2 ends of the spine and a center point
    special_fold_ids = [24//4, 24 // 4 * 3] ## 2 sides of the center point
    source_b = vtk.vtkPoints()
    target_b = vtk.vtkPoints()
    for i in range(math.floor((bot_srep.GetNumberOfCells() -24) / 2), bot_srep.GetNumberOfCells()):
        base_pt_id = i * 2
        bdry_pt_id = i * 2 + 1
        s_pt = bot_srep.GetPoint(base_pt_id)
        b_pt = bot_srep.GetPoint(bdry_pt_id)
        source_b.InsertNextPoint(s_pt)
        target_b.InsertNextPoint(b_pt)
    tpsB = vtk.vtkThinPlateSplineTransform()
    tpsB.SetSourceLandmarks(source_b)
    tpsB.SetTargetLandmarks(target_b)
    tpsB.SetBasisToR()
    tpsB.Modified()

    # finding closest 2D s-rep points on 2D s-rep
    numSamp = numSamp + 1
    eps = 0.5
    for i in range(0, bot_srep.GetNumberOfCells()):# math.floor((bot_srep.GetNumberOfCells() -24) / 2)):
        base_pt_id = i * 2
        bdry_pt_id = i * 2 + 1
        s_pt = bot_srep.GetPoint(base_pt_id)
        b_pt = bot_srep.GetPoint(bdry_pt_id)
        base_pt = np.array(s_pt)
        bdry_pt = np.array(b_pt)
        radius = np.linalg.norm(bdry_pt - base_pt)
        direction = (bdry_pt - base_pt) / radius

        dist, ind = findClosetPt([s_pt[0], s_pt[1], s_pt[2]], ptsOnSkel, numSamp)
        # # ax.scatter(s_pt[0], s_pt[1], s_pt[2], color='r')

        spk0 = math.floor(ind[0] / numSamp)
        spk1 = math.floor(ind[1] / numSamp)

        # if ind[0] % numSamp == 0:
        #     continue

        # if spk0 == len(rXs) -1:
        #     # ax.scatter(ptsOnSkel[ind[0]][0], ptsOnSkel[ind[0]][1], ptsOnSkel[ind[0]][2], color='r')
        #     # ax.scatter(ptsOnSkel[ind[1]][0], ptsOnSkel[ind[1]][1], ptsOnSkel[ind[1]][2], color='r')


        if spk0 % 2 == 1:
            if spk1 < spk0:
                uDir = [ptsOnSkel[ind[1]][0] - ptsOnSkel[ind[0]][0], ptsOnSkel[ind[1]][1] - ptsOnSkel[ind[0]][1], ptsOnSkel[ind[1]][2] - ptsOnSkel[ind[0]][2]]
            else:
                uDir = [ptsOnSkel[ind[0]][0] - ptsOnSkel[ind[1]][0], ptsOnSkel[ind[0]][1] - ptsOnSkel[ind[1]][1], ptsOnSkel[ind[0]][2] - ptsOnSkel[ind[1]][2]]
        else:
            if spk1 > spk0:
                uDir = [ptsOnSkel[ind[1]][0] - ptsOnSkel[ind[0]][0], ptsOnSkel[ind[1]][1] - ptsOnSkel[ind[0]][1], ptsOnSkel[ind[1]][2] - ptsOnSkel[ind[0]][2]]
            else:
                uDir = [ptsOnSkel[ind[0]][0] - ptsOnSkel[ind[1]][0], ptsOnSkel[ind[0]][1] - ptsOnSkel[ind[1]][1], ptsOnSkel[ind[0]][2] - ptsOnSkel[ind[1]][2]]
        if spk0 == 0:
            if spk1 % 2 == 1:
                uDir = [ptsOnSkel[ind[0]][0] - ptsOnSkel[ind[1]][0], ptsOnSkel[ind[0]][1] - ptsOnSkel[ind[1]][1], ptsOnSkel[ind[0]][2] - ptsOnSkel[ind[1]][2]]
            else:
                uDir = [ptsOnSkel[ind[1]][0] - ptsOnSkel[ind[0]][0], ptsOnSkel[ind[1]][1] - ptsOnSkel[ind[0]][1], ptsOnSkel[ind[1]][2] - ptsOnSkel[ind[0]][2]]           
        if ind[0] >= len(ptsOnSkel) - numSamp:
            if spk1 % 2 == 1:
                uDir = [ptsOnSkel[ind[1]][0] - ptsOnSkel[ind[0]][0], ptsOnSkel[ind[1]][1] - ptsOnSkel[ind[0]][1], ptsOnSkel[ind[1]][2] - ptsOnSkel[ind[0]][2]]
            else:
                uDir = [ptsOnSkel[ind[0]][0] - ptsOnSkel[ind[1]][0], ptsOnSkel[ind[0]][1] - ptsOnSkel[ind[1]][1], ptsOnSkel[ind[0]][2] - ptsOnSkel[ind[1]][2]]


        # uDir = [ptsOnSkel[ind[1]][0] - ptsOnSkel[ind[0]][0], ptsOnSkel[ind[1]][1] - ptsOnSkel[ind[0]][1], ptsOnSkel[ind[1]][2] - ptsOnSkel[ind[0]][2]]
        # print(ptsOnSkel[125])
        # print(ptsOnSkel[130])
        length = math.sqrt(uDir[0]**2 + uDir[1]**2 + uDir[2]**2)
        uDir[0] = uDir[0]*eps / length
        uDir[1] = uDir[1]*eps / length
        uDir[2] = uDir[2]*eps / length


        if ind[0] % numSamp != numSamp-1 and ind[0] % numSamp != 0:
            tao = [ptsOnSkel[ind[0]+1][0] - ptsOnSkel[ind[0]-1][0], ptsOnSkel[ind[0]+1][1] - ptsOnSkel[ind[0]-1][1], ptsOnSkel[ind[0]+1][2] - ptsOnSkel[ind[0]-1][2]]
        if ind[0] % numSamp == 0:
            tao = [ptsOnSkel[ind[0]+1][0] - ptsOnSkel[ind[0]][0], ptsOnSkel[ind[0]+1][1] - ptsOnSkel[ind[0]][1], ptsOnSkel[ind[0]+1][2] - ptsOnSkel[ind[0]][2]]
        if ind[0] % numSamp == numSamp-1:
            tao = [ptsOnSkel[ind[0]][0] - ptsOnSkel[ind[0]-1][0], ptsOnSkel[ind[0]][1] - ptsOnSkel[ind[0]-1][1], ptsOnSkel[ind[0]][2] - ptsOnSkel[ind[0]-1][2]]

        if ind[1] % numSamp != numSamp-1 and ind[1] % numSamp != 0:
            tao1 = [ptsOnSkel[ind[1]+1][0] - ptsOnSkel[ind[1]-1][0], ptsOnSkel[ind[1]+1][1] - ptsOnSkel[ind[1]-1][1], ptsOnSkel[ind[1]+1][2] - ptsOnSkel[ind[1]-1][2]]
        if ind[1] % numSamp == 0:
            tao1 = [ptsOnSkel[ind[1]+1][0] - ptsOnSkel[ind[1]][0], ptsOnSkel[ind[1]+1][1] - ptsOnSkel[ind[1]][1], ptsOnSkel[ind[1]+1][2] - ptsOnSkel[ind[1]][2]]
        if ind[1] % numSamp == numSamp-1:
            tao1 = [ptsOnSkel[ind[1]][0] - ptsOnSkel[ind[1]-1][0], ptsOnSkel[ind[1]][1] - ptsOnSkel[ind[1]-1][1], ptsOnSkel[ind[1]][2] - ptsOnSkel[ind[1]-1][2]]

        dist1 = (s_pt[0]-ptsOnSkel[ind[0]][0])**2+(s_pt[1]-ptsOnSkel[ind[0]][1])**2+(s_pt[2]-ptsOnSkel[ind[0]][2])**2
        dist2 = (s_pt[0]-ptsOnSkel[ind[1]][0])**2+(s_pt[1]-ptsOnSkel[ind[1]][1])**2+(s_pt[2]-ptsOnSkel[ind[1]][2])**2
        ratio1 = dist1 / (dist1+dist2)
        ratio2 = dist2 / (dist1+dist2)
        finTao = [(ratio1*tao[0]+ratio2*tao1[0])/2, (ratio1*tao[1]+ratio2*tao1[1])/2, (ratio1*tao[2]+ratio2*tao1[2])/2]
        lenT = math.sqrt(finTao[0]**2 + finTao[1]**2 + finTao[2]**2)
        finTao[0] = finTao[0]*eps / lenT
        finTao[1] = finTao[1]*eps / lenT
        finTao[2] = finTao[2]*eps / lenT

        tU = [s_pt[0]+finTao[0], s_pt[1]+finTao[1], s_pt[2]+finTao[2]]
        tD = [s_pt[0]-finTao[0], s_pt[1]-finTao[1], s_pt[2]-finTao[2]]
        uR =  [s_pt[0]+uDir[0], s_pt[1]+uDir[1], s_pt[2]+uDir[2]]
        uL =  [s_pt[0]-uDir[0], s_pt[1]-uDir[1], s_pt[2]-uDir[2]]

        spokeTU = tps.TransformPoint(tU)
        if i >= math.floor((bot_srep.GetNumberOfCells() -24) / 2):
            spokeTU = tpsB.TransformPoint(tU)
            # ax.scatter(spokeTU[0], spokeTU[1], spokeTU[2], color='k')
        vecTU = [spokeTU[0] - tU[0], spokeTU[1] - tU[1], spokeTU[2] - tU[2]]
        lengthTU = math.sqrt(vecTU[0]**2 + vecTU[1]**2 + vecTU[2]**2)
        vecTU = [vecTU[0]/lengthTU, vecTU[1]/lengthTU, vecTU[2]/lengthTU]

        spokeTD = tps.TransformPoint(tD)
        if i >= math.floor((bot_srep.GetNumberOfCells() -24) / 2):
            spokeTD = tpsB.TransformPoint(tD)
        vecTD = [spokeTD[0] - tD[0], spokeTD[1] - tD[1], spokeTD[2] - tD[2]]
        lengthTD = math.sqrt(vecTD[0]**2 + vecTD[1]**2 + vecTD[2]**2)
        vecTD = [vecTD[0]/lengthTD, vecTD[1]/lengthTD, vecTD[2]/lengthTD]

        spoke = tps.TransformPoint(s_pt)
        if i >= math.floor((bot_srep.GetNumberOfCells() -24) / 2):
            spoke = tpsB.TransformPoint(s_pt)
        vec = [spoke[0] - s_pt[0], spoke[1] - s_pt[1], spoke[2] - s_pt[2]]
        length = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
        vec = [vec[0]/length, vec[1]/length, vec[2]/length]

        spokeUR = tps.TransformPoint(uR)
        if i >= math.floor((bot_srep.GetNumberOfCells() -24) / 2):
            spokeUR = tpsB.TransformPoint(uR)
        vecUR = [spokeUR[0] - uR[0], spokeUR[1] - uR[1], spokeUR[2] - uR[2]]
        lengthUR = math.sqrt(vecUR[0]**2 + vecUR[1]**2 + vecUR[2]**2)
        vecUR = [vecUR[0]/lengthUR, vecUR[1]/lengthUR, vecUR[2]/lengthUR]

        spokeUL = tps.TransformPoint(uL)
        if i >= math.floor((bot_srep.GetNumberOfCells() -24) / 2):
            spokeUL = tpsB.TransformPoint(uL)
        vecUL = [spokeUL[0] - uL[0], spokeUL[1] - uL[1], spokeUL[2] - uL[2]]
        lengthUL = math.sqrt(vecUL[0]**2 + vecUL[1]**2 + vecUL[2]**2)
        vecUL = [vecUL[0]/lengthUL, vecUL[1]/lengthUL, vecUL[2]/lengthUL]

        spokeSam = 2

        for j in range(0, spokeSam):
            # if j == 0:
            #     continue
            unit = length/(spokeSam-1)
            spokePt = [s_pt[0] + vec[0]*j*unit, s_pt[1] + vec[1]*j*unit, s_pt[2] + vec[2]*j*unit]

            unit = lengthTU/(spokeSam-1)
            tUpt = [tU[0] + vecTU[0]*j*unit, tU[1] + vecTU[1]*j*unit, tU[2] + vecTU[2]*j*unit]

            unit = lengthTD/(spokeSam-1)
            tDpt = [tD[0] + vecTD[0]*j*unit, tD[1] + vecTD[1]*j*unit, tD[2] + vecTD[2]*j*unit]

            unit = lengthUR/(spokeSam-1)
            uRpt = [uR[0] + vecUR[0]*j*unit, uR[1] + vecUR[1]*j*unit, uR[2] + vecUR[2]*j*unit]

            unit = lengthUL/(spokeSam-1)
            uLpt = [uL[0] + vecUL[0]*j*unit, uL[1] + vecUL[1]*j*unit, uL[2] + vecUL[2]*j*unit]

            # if spk0 == len(rXs) -1:
            #     # ax.scatter(spokePt[0], spokePt[1], spokePt[2], color='k')

            spokeT = tps2.TransformPoint(spokePt)
            tUtrans = tps2.TransformPoint(tUpt)
            tDtrans = tps2.TransformPoint(tDpt)
            uRtrans = tps2.TransformPoint(uRpt)
            uLtrans = tps2.TransformPoint(uLpt)

            # # ax.scatter(tDpt[0], tDpt[1], tDpt[2], color='k')
            # # ax.scatter(tUpt[0], tUpt[1], tUpt[2], color='k')
            # # ax.scatter(uRpt[0], uRpt[1], uRpt[2], color='k')
            # # ax.scatter(uLpt[0], uLpt[1], uLpt[2], color='k')
            # # ax.scatter(tUtrans[0], tUtrans[1], tUtrans[2], color='k')
            # # ax.scatter(tDtrans[0], tDtrans[1], tDtrans[2], color='r')
            # if j == 4:
            #     # ax.scatter(uRtrans[0], uRtrans[1], uRtrans[2], color='k')
            #     # ax.scatter(uLtrans[0], uLtrans[1], uLtrans[2], color='r')

            # # ax.scatter(tU[0]+ tU[0]*4*unit, tU[1]+ vecTU[1]*4*unit, tU[2]+ vecTU[2]*4*unit, color='k')
            # # ax.scatter(ptsOnSkel[ind[0]][0], ptsOnSkel[ind[0]][1], ptsOnSkel[ind[0]][2], color='r')
            # # ax.scatter(ptsOnSkel[ind[0]][0], ptsOnSkel[ind[0]][1], ptsOnSkel[ind[0]][2], color='r')
            # # ax.scatter(ptsOnSkel[ind[0]][0], ptsOnSkel[ind[0]][1], ptsOnSkel[ind[0]][2], color='r')

            uVec = [uRtrans[0] - uLtrans[0], uRtrans[1] - uLtrans[1], uRtrans[2] - uLtrans[2]]
            uDist = np.linalg.norm(uVec)
            uVec = [uVec[0]/uDist, uVec[1]/uDist, uVec[2]/uDist]

            tVec = [tUtrans[0] - tDtrans[0], tUtrans[1] - tDtrans[1], tUtrans[2] - tDtrans[2]]
            tDist = np.linalg.norm(tVec)
            tVec = [tVec[0]/tDist, tVec[1]/tDist, tVec[2]/tDist]

            norm1 = np.cross(uVec, tVec)
            norm2 = np.cross(tVec, uVec)


            finU = [spokeT[0] + uVec[0], spokeT[1] + uVec[1], spokeT[2] + uVec[2]]

            finN1 = [spokeT[0] + norm1[0], spokeT[1] + norm1[1], spokeT[2] + norm1[2]]
            finN2 = [spokeT[0] + norm2[0], spokeT[1] + norm2[1], spokeT[2] + norm2[2]]
            # # ax.scatter(finN1[0], finN1[1], finN1[2], color='b')
            # # ax.scatter(finN2[0], finN2[1], finN2[2], color='r')


            if j == spokeSam-1:
                nextPt = [s_pt[0] + vec[0]*(j-1)*unit, s_pt[1] + vec[1]*(j-1)*unit, s_pt[2] + vec[2]*(j-1)*unit]
                nextPt = tps2.TransformPoint(nextPt)
            else:
                nextPt = [s_pt[0] + vec[0]*(j+1)*unit, s_pt[1] + vec[1]*(j+1)*unit, s_pt[2] + vec[2]*(j+1)*unit]
                nextPt = tps2.TransformPoint(nextPt)
                # # ax.scatter(nextPt[0], nextPt[1], nextPt[2], color='k')
            normToPt1 = [nextPt[0] - finN1[0], nextPt[1] - finN1[1], nextPt[2] - finN1[2]]
            normDist1 = np.linalg.norm(normToPt1)
            normToPt2 = [nextPt[0] - finN2[0], nextPt[1] - finN2[1], nextPt[2] - finN2[2]]
            normDist2 = np.linalg.norm(normToPt2)

            boo1 = False
            if normDist1 < normDist2:
                finN = finN1
                vec3 = np.cross(uVec, norm1)
                boo1 = True
            else:
                finN = finN2
                vec3 = np.cross(uVec, norm2)

            if j == spokeSam-1:
                if boo1:
                    finN = finN2
                    vec3 = np.cross(uVec, norm2)
                else:
                    finN = finN1
                    vec3 = np.cross(uVec, norm1)


            finT1 = [spokeT[0] + vec3[0], spokeT[1] + vec3[1], spokeT[2] + vec3[2]]
            normToPt1 = [tUtrans[0] - finT1[0], tUtrans[1] - finT1[1], tUtrans[2] - finT1[2]]
            normDist1 = np.linalg.norm(normToPt1)

            finT2 = [spokeT[0] - vec3[0], spokeT[1] - vec3[1], spokeT[2] - vec3[2]]
            normToPt2 = [tUtrans[0] - finT2[0], tUtrans[1] - finT2[1], tUtrans[2] - finT2[2]]
            normDist2 = np.linalg.norm(normToPt2)

            boo1 = False
            if normDist1 < normDist2:
                finT = finT1
            else:
                finT = finT2
                vec3 = np.cross(uVec, norm2)

            if j == 1:
                continue
            if i >= math.floor((bot_srep.GetNumberOfCells() -24) / 2):
                highlight_fitFrames_ends, highlight_fitFrames_lines = form_frames(spokeT, finU, finT, finN, highlight_fitFrames_ends, highlight_fitFrames_lines)
                if int(i-math.floor((bot_srep.GetNumberOfCells() -24) / 2)) in special_ids:
                    highlight_special_fitFrames_ends, highlight_special_fitFrames_lines = form_frames(spokeT, finU, finT, finN, highlight_special_fitFrames_ends, highlight_special_fitFrames_lines)
                elif int(i-math.floor((bot_srep.GetNumberOfCells() -24))) in special_fold_ids:
                    special_fitFrames_ends, special_fitFrames_lines = form_frames(spokeT, finU, finT, finN, special_fitFrames_ends, special_fitFrames_lines)
            else:
                fitFrames_ends, fitFrames_lines = form_frames(spokeT, finU, finT, finN, fitFrames_ends, fitFrames_lines)
                if i in special_ids:
                    special_fitFrames_ends, special_fitFrames_lines = form_frames(spokeT, finU, finT, finN, special_fitFrames_ends, special_fitFrames_lines)

    fitFrames.SetPoints(fitFrames_ends)
    fitFrames.SetLines(fitFrames_lines)
    fitFrames.Modified()
    special_fitFrames.SetPoints(special_fitFrames_ends)
    special_fitFrames.SetLines(special_fitFrames_lines)
    special_fitFrames.Modified()

    highlight_fitFrames.SetPoints(highlight_fitFrames_ends)
    highlight_fitFrames.SetLines(highlight_fitFrames_lines)
    highlight_fitFrames.Modified()

    highlight_special_fitFrames.SetPoints(highlight_special_fitFrames_ends)
    highlight_special_fitFrames.SetLines(highlight_special_fitFrames_lines)
    highlight_special_fitFrames.Modified()

    frames_appender = vtk.vtkAppendPolyData()
    frames_appender.AddInputData(special_fitFrames)
    frames_appender.AddInputData(highlight_special_fitFrames)
    frames_appender.Update()

    frame_collection = frames_appender.GetOutput()
    writer2 = vtk.vtkPolyDataWriter()
    writer2.SetInputData(fitFrames)
    writer2.SetFileName(output)
    writer2.Write()

    other_frames_appender = vtk.vtkAppendPolyData()
    other_frames_appender.AddInputData(fitFrames)
    other_frames_appender.AddInputData(highlight_fitFrames)
    other_frames_appender.Update()
#    plt.show()

#    viz.overlay_polydata(frame_collection, top_mesh)#, special_fitFrames)#, appendum=other_frames_appender.GetOutput())
#    viz.overlay_polydata(highlight_special_fitFrames, top_mesh, special_fitFrames, appendum=other_frames_appender.GetOutput())

    return frame_collection
if __name__ == "__main__":
    fit_local_frames(output='/playpen/workspace/Simulate_Shapes/data/frames2.vtk')
    print('Done')