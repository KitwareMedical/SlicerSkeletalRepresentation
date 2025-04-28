"""
Compute surface mesh vertices and their centers of curvature (coc)
"""
import vtk
import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib
import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt
from glob import glob
import pyvista as pv
import vtk
import pyacvd
from scipy.stats import norm
from glob import glob

def compute_coc_and_vertices_from_mesh(inputmesh, previous_vertices=None, do_clustering=False, do_smoothing=True):
    """
    Get the desired skeletal spine endpoints from a mesh.
    First find the "vertices" of the mesh (points of maximal gaussian curvature)
    Then get the normal vector and the smaller principal
    curvature at those points to find the focal points,
    which are the spine endpoints.

    clustering: remeshing
    smoothing: taubin smoothing

    return recommended spine endpoints, lens to vertices, and then vertices
    """
    if do_clustering:

        clus = pyacvd.Clustering(pv.PolyData(inputmesh))
        clus.subdivide(3)
        clus.cluster(inputmesh.GetNumberOfPoints() * 3)
        inputmesh = clus.create_mesh()

    if do_smoothing:
        taubin_smooth = vtk.vtkWindowedSincPolyDataFilter()
        taubin_smooth.SetInputData(inputmesh)
        taubin_smooth.SetNumberOfIterations(20)
        taubin_smooth.BoundarySmoothingOff()
        taubin_smooth.FeatureEdgeSmoothingOff()
        taubin_smooth.SetPassBand(0.01)
        taubin_smooth.NonManifoldSmoothingOn()
        taubin_smooth.NormalizeCoordinatesOn()

        taubin_smooth.Update()
        smoothed_mesh = taubin_smooth.GetOutput()
    else:
        smoothed_mesh = inputmesh

    # print(f"num points in remesh {smoothed_mesh.GetNumberOfPoints()}")

    # get normals, and 4 measures of curvature

    normal_generator = vtk.vtkPolyDataNormals()
    normal_generator.SetInputData(smoothed_mesh)
    normal_generator.SplittingOff()
    normal_generator.ComputePointNormalsOn()
    normal_generator.ComputeCellNormalsOff()
    normal_generator.Update()
    normals = normal_generator.GetOutput().GetPointData().GetNormals()

    curvatures = vtk.vtkCurvatures()
    curvatures.SetInputData(smoothed_mesh)

    curvatures.SetCurvatureTypeToGaussian()
    curvatures.Update()
    gauss_curvatures = curvatures.GetOutput().GetPointData().GetScalars()

    curvatures.SetCurvatureTypeToMean()
    curvatures.Update()
    mean_curvatures = curvatures.GetOutput().GetPointData().GetScalars()

    curvatures.SetCurvatureTypeToMaximum()
    curvatures.Update()
    max_curvatures = curvatures.GetOutput().GetPointData().GetScalars()

    curvatures.SetCurvatureTypeToMinimum()
    curvatures.Update()
    min_curvatures = curvatures.GetOutput().GetPointData().GetScalars()

    print(gauss_curvatures.GetValue(0), mean_curvatures.GetValue(0))
    # quit()

    smoothed_mesh.BuildLinks()

    # get edge binary matrix
    edge_matrix = np.zeros((smoothed_mesh.GetNumberOfPoints(), smoothed_mesh.GetNumberOfPoints()))
    for i in range(len(edge_matrix)):
        for j in range(len(edge_matrix)):
            edge_matrix[i,j] = int(smoothed_mesh.IsEdge(i,j))

    # smooth mean curvature
    mean_curvatures_arr = np.array(mean_curvatures)

    for iterations in range(2):
        mean_curvatures_arr = (mean_curvatures_arr / 2) + (((edge_matrix @ mean_curvatures_arr.reshape(-1,1)).reshape(-1) / np.sum(edge_matrix, axis=1)) / 2)


    """
    The vertex-finding algorithm is summarized here. Recall that a vertex is 
    a maximum of convex curvature. Also, for the mandible / hotdog shape, one
    vertex will have negative x coord and one will have positive x coord.
    We assume one vertex will have positive x coord and the other will have negative 
    x coord from here.

    1. Find all local maxima of curvature
    2. Separate into two buckets: postive x coordinate, and negative x coordinate
    3. Get absolute maxima of curvature in those buckets - those are vertices
    3.5 If we have previous vertices, choose points that are closer to those previous vertices
    4. Use curvature at vertices to approximate centers of curvature

    """

    inds_curvatures = []
    # Find all the relative maxima of gaussian curvature and store them in inds_curvatures
    for i in range(smoothed_mesh.GetNumberOfPoints()):
        is_relmax = True
        G = gauss_curvatures.GetValue(i)
        H = mean_curvatures_arr[i]
        if H == 0.0: continue # perfect saddle or flat point
        # check if point is more curved than all its neighbors
        for j in range(smoothed_mesh.GetNumberOfPoints()):
            if i == j: continue
            if edge_matrix[i, j]:
                # is_relmax *= G > gauss_curvatures.GetValue(j)
                is_relmax *= H > abs(mean_curvatures_arr[j])
        # if is_relmax:
        if True:
            k1 = max_curvatures.GetValue(i)
            k2 = min_curvatures.GetValue(i)
            # k1 = H + np.sqrt(H**2 - this_curvature)
            # k2 = H - np.sqrt(H**2 - this_curvature)
            inds_curvatures.append([i, G, H, k1, k2])


    if len(inds_curvatures) == 0:
        print("failed to find relative maxima of Gaussian curvature")

    # Now we take into account the distance to the previous vertex    
    print("before: ",len(inds_curvatures))

    if previous_vertices is not None:
        # print("using previous vertices")
        # print(previous_vertices)
        old_inds_curvatures = [x for x in inds_curvatures]
        inds_curvatures = []
        tol = 0.1
        for row in old_inds_curvatures:
            pt = np.array(smoothed_mesh.GetPoint(row[0]))[None,:]
            print(np.min(np.linalg.norm(previous_vertices-pt, axis=1)))
            if np.min(np.linalg.norm(previous_vertices-pt, axis=1)) < tol:
                inds_curvatures.append(row)

    print("after: ",len(inds_curvatures))

    # inds_curvatures = sorted(inds_curvatures, key=lambda x: -x[1]) #sort in decending gaussian curvature
    inds_curvatures = sorted(inds_curvatures, key=lambda x: -x[2]) #sort in decending mean curvature

    # locs = [smoothed_mesh.GetPoint(x[0]) for x in inds_curvatures]
    locs = np.array([smoothed_mesh.GetPoint(i) for i in range(smoothed_mesh.GetNumberOfPoints())])

    # separate points with positive x and points with neg x and pos x
    # for hippo, this is actually y, not x. hotdog uses x
    separate_dim = 0
    inds_curvatures_negx = []
    inds_curvatures_posx = []
    for ic in inds_curvatures:
        j = ic[0]
        if locs[j,separate_dim] < 0:
            inds_curvatures_negx.append(ic)
        else:
            inds_curvatures_posx.append(ic)

    # candidate_inds_???x are still sorted at this point
    print("-------------------------")
    print(len(inds_curvatures_negx))
    print(len(inds_curvatures_posx))
    candidate_inds_locs = [
        [inds_curvatures_negx[0][0], smoothed_mesh.GetPoint(inds_curvatures_negx[0][0])],
        [inds_curvatures_posx[0][0], smoothed_mesh.GetPoint(inds_curvatures_posx[0][0])]
        ]
    points = [smoothed_mesh.GetPoint(i[0]) for i in candidate_inds_locs]
    points = np.asarray(points)

    # Now find centers of curvature from vertices

    all_locs = [smoothed_mesh.GetPoint(i) for i in range(smoothed_mesh.GetNumberOfPoints())]
    all_locs = np.array(all_locs)
    nearby_inds = []
    tol = 0.1
    for cl in candidate_inds_locs:
        arr = np.linalg.norm(all_locs - cl[1], axis=1) < tol
        nearby_inds.append(
            np.argwhere(arr)
        )

    nearby_normals = []
    for ni in nearby_inds:
        nearby_normals.append(
            np.array([normals.GetTuple3(i[0]) for i in ni])
        )

    avg_normals = [np.mean(nn, axis=0) for nn in nearby_normals]

    candidate_avg_normals = avg_normals[:2]
    candidate_avg_normals = np.array(candidate_avg_normals)

    # "newpts" are COC, "points" are vertices.

    newpts = points - 0.07 * candidate_avg_normals
    len_from_vertex = 0.07 * np.linalg.norm(candidate_avg_normals, axis=1)
    len_from_vertex = len_from_vertex.reshape((2,1))

    if newpts[0,separate_dim] > newpts[-1,separate_dim]:
        newpts = newpts[::-1]
    if points[0,separate_dim] > points[-1,separate_dim]:
        points = points[::-1]

    # coc, lens between coc and vertices, vertices, smoothed mesh
    assert(points[0,separate_dim] < points[1,separate_dim])
    assert(newpts[0,separate_dim] < newpts[1,separate_dim])
    return newpts, len_from_vertex, points, smoothed_mesh

def compute_coc_and_vertices_from_mesh_by_density(inputmesh, previous_vertices=None, do_clustering=False, do_smoothing=True):
    """
    Get the desired skeletal spine endpoints from a mesh.
    First find the mesh node elements which have high density
    of other mesh node elements.
    Then get the normal vector and the smaller principal
    curvature at those points to find the focal points,
    which are the spine endpoints.

    smoothing: taubin smoothing

    return recommended spine endpoints, lens to vertices, and then vertices
    """
    if do_clustering:

        clus = pyacvd.Clustering(pv.PolyData(inputmesh))
        clus.subdivide(3)
        clus.cluster(inputmesh.GetNumberOfPoints() * 3)
        inputmesh = clus.create_mesh()

    if do_smoothing:
        taubin_smooth = vtk.vtkWindowedSincPolyDataFilter()
        taubin_smooth.SetInputData(inputmesh)
        taubin_smooth.SetNumberOfIterations(20)
        taubin_smooth.BoundarySmoothingOff()
        taubin_smooth.FeatureEdgeSmoothingOff()
        taubin_smooth.SetPassBand(0.01)
        taubin_smooth.NonManifoldSmoothingOn()
        taubin_smooth.NormalizeCoordinatesOn()

        taubin_smooth.Update()
        smoothed_mesh = taubin_smooth.GetOutput()
    else:
        smoothed_mesh = inputmesh

    # print(f"num points in remesh {smoothed_mesh.GetNumberOfPoints()}")

    # get normals, and 4 measures of curvature

    normal_generator = vtk.vtkPolyDataNormals()
    normal_generator.SetInputData(smoothed_mesh)
    normal_generator.SplittingOff()
    normal_generator.ComputePointNormalsOn()
    normal_generator.ComputeCellNormalsOff()
    normal_generator.Update()
    normals = normal_generator.GetOutput().GetPointData().GetNormals()

    curvatures = vtk.vtkCurvatures()
    curvatures.SetInputData(smoothed_mesh)

    curvatures.SetCurvatureTypeToGaussian()
    curvatures.Update()
    gauss_curvatures = curvatures.GetOutput().GetPointData().GetScalars()

    curvatures.SetCurvatureTypeToMean()
    curvatures.Update()
    mean_curvatures = curvatures.GetOutput().GetPointData().GetScalars()

    curvatures.SetCurvatureTypeToMaximum()
    curvatures.Update()
    max_curvatures = curvatures.GetOutput().GetPointData().GetScalars()

    curvatures.SetCurvatureTypeToMinimum()
    curvatures.Update()
    min_curvatures = curvatures.GetOutput().GetPointData().GetScalars()

    # print(gauss_curvatures.GetValue(0), mean_curvatures.GetValue(0))
    # quit()

    smoothed_mesh.BuildLinks()

    # get edge binary matrix
    edge_matrix = np.zeros((smoothed_mesh.GetNumberOfPoints(), smoothed_mesh.GetNumberOfPoints()))
    for i in range(len(edge_matrix)):
        for j in range(len(edge_matrix)):
            edge_matrix[i,j] = int(smoothed_mesh.IsEdge(i,j))

    # smooth mean curvature
    mean_curvatures_arr = np.array(mean_curvatures)

    for iterations in range(2):
        mean_curvatures_arr = (mean_curvatures_arr / 2) + (((edge_matrix @ mean_curvatures_arr.reshape(-1,1)).reshape(-1) / np.sum(edge_matrix, axis=1)) / 2)


    """
    The vertex-finding algorithm is summarized here. Recall that a vertex is 
    a maximum of convex curvature. Also, for the mandible / hotdog shape, one
    vertex will have negative x coord and one will have positive x coord.
    We assume one vertex will have positive x coord and the other will have negative 
    x coord from here.

    1. Find all local maxima of curvature
    2. Separate into two buckets: postive x coordinate, and negative x coordinate
    3. Get absolute maxima of curvature in those buckets - those are vertices
    3.5 If we have previous vertices, choose points that are closer to those previous vertices
    4. Use curvature at vertices to approximate centers of curvature

    """

    locs = np.array([smoothed_mesh.GetPoint(i) for i in range(smoothed_mesh.GetNumberOfPoints())])

    # compute kernel bandwith by heuristic
    bbox_dims = np.max(locs, axis=0) - np.min(locs, axis=0)
    kernel_bandwith = sorted(bbox_dims)[0] / 20
    kernel = norm(loc=0, scale=kernel_bandwith)

    inds_density_scores = []
    for i in range(locs.shape[0]):
        m = locs[i].reshape((1,3))
        dists_to_locs = np.linalg.norm(m - locs, axis=1)
        score = np.sum(kernel.pdf(dists_to_locs))
        inds_density_scores.append([i, score])

    inds_density_scores = sorted(inds_density_scores, key=lambda x: -x[1]) #sort in decending density score

    # separate points with positive x and points with neg x and pos x
    inds_scores_negx = []
    inds_scores_posx = []
    for j, score in inds_density_scores:
        if locs[j,0] < 0:
            inds_scores_negx.append([j, score])
        else:
            inds_scores_posx.append([j, score])

    # candidate_inds_???x are still sorted at this point
    candidate_inds_locs = [
        [inds_scores_negx[0][0], smoothed_mesh.GetPoint(inds_scores_negx[0][0])],
        [inds_scores_posx[0][0], smoothed_mesh.GetPoint(inds_scores_posx[0][0])]
        ]

    points = [smoothed_mesh.GetPoint(i[0]) for i in candidate_inds_locs]
    points = np.asarray(points)

    # Now find centers of curvature from vertices
    newpts, len_from_vertex = compute_coc_from_mesh_and_vertices(smoothed_mesh, candidate_inds_locs, normals)

    # coc, lens between coc and vertices, vertices
    return newpts, len_from_vertex, points, smoothed_mesh, candidate_inds_locs

def compute_coc_from_mesh_and_vertices(smoothed_mesh, candidate_inds_locs, normals):
    """
    compute coc using vertices and mesh 
    """
    all_locs = [smoothed_mesh.GetPoint(i) for i in range(smoothed_mesh.GetNumberOfPoints())]
    all_locs = np.array(all_locs)

    points = [smoothed_mesh.GetPoint(i[0]) for i in candidate_inds_locs]
    points = np.asarray(points)

    nearby_inds = []
    tol = 0.1
    for cl in candidate_inds_locs:
        arr = np.linalg.norm(all_locs - cl[1], axis=1) < tol
        nearby_inds.append(
            np.argwhere(arr)
        )

    nearby_normals = []
    for ni in nearby_inds:
        nearby_normals.append(
            np.array([normals.GetTuple3(i) for i in ni])
        )

    avg_normals = [np.mean(nn, axis=0) for nn in nearby_normals]

    candidate_avg_normals = avg_normals[:2]
    candidate_avg_normals = np.array(candidate_avg_normals)

    # "newpts" are COC, "points" are vertices.

    newpts = points - 0.1 * candidate_avg_normals
    len_from_vertex = 0.1 * np.linalg.norm(candidate_avg_normals, axis=1)
    len_from_vertex = len_from_vertex.reshape((2,1))

    if newpts[0,0] > newpts[-1,0]:
        newpts = newpts[::-1]
    if points[0,0] > points[-1,0]:
        points = points[::-1]

    # coc, lens between coc and vertices, vertices
    return newpts, len_from_vertex
