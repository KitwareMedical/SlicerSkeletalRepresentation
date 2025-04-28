"""
Used for generating the perfect ellipsoid and its srep, vertices, and crest from the quasi ellipsoid mesh.
"""
try:
    from nick_srep import Srep
except ImportError:
    from shanapy_private.nick_srep import Srep
import numpy as np
import pyvista as pv
import vtk

def get_crest_points_of_perfect_ellipsoid(radii, pe_mesh):
    """
    Here we find the indices of the crest curve by a closest point analysis
    """
    r1 = radii[0]
    r2 = radii[1]
    t = np.linspace(0, 2*np.pi, 1000)
    p_crest_curve = np.array([r1 * np.cos(t), r2 * np.sin(t), np.zeros_like(t)])
    p_crest_curve = np.transpose(p_crest_curve)
    crest_indices = []
    pts = []
    for c in p_crest_curve:
        d = np.inf
        candidate_ind = None

        for i, p in enumerate(pe_mesh.points):
            new = np.linalg.norm(c - p)
            if new < d:
                d = new
                candidate_ind = i

        if crest_indices == []:
            crest_indices = [candidate_ind]
        elif candidate_ind not in crest_indices:
            # check if the x,y coords are the same as another point (directly on top of )
            is_unique = True
            for idx in crest_indices:
                if np.linalg.norm(pe_mesh.points[idx][:2]-pe_mesh.points[candidate_ind][:2])<1e-3:
                    is_unique = False
            if is_unique: 
                crest_indices.append(candidate_ind)
    
    return crest_indices, p_crest_curve

def get_perfect_ellipsoid_and_srep_from_quasi(mesh, num_spine_points=20, num_skeletal_onionskins=4):
    """
    Takes a quasi-ellipsoid mesh and performs a PCA analysis
    to find principal radii of best-fitting "perfect" ellipsoid, then uses a 
    bounding box method to find the lengths of those radii.


    """

    points = np.copy(mesh.points)
    center = np.mean(points, axis=0)[np.newaxis, :]
    centered_points = points - center
    w, v = np.linalg.eig(np.matmul(points.T, points))
    idx = w.argsort()[::-1]
    w = w[idx]
    v = v[:, idx]

    # By this point, w contains the principal values in deacreasing order
    # and v contains the associated eigenvectors as columns

    rotated_points = centered_points @ v

    # bounding box method

    r1, r2, r3 = (np.max(rotated_points, axis=0) - np.min(rotated_points, axis=0)) / 2

    # scale radii so that ellipsoid volume matches mesh volume
    ellipsoid_volume = 4 / 3.0 * np.pi * r1 * r2 * r3
    pv_vol = mesh.volume
    volume_factor = pow(pv_vol / ellipsoid_volume, 1.0 / 3.0)
    r1 *= volume_factor
    r2 *= volume_factor
    r3 *= volume_factor

    radii = (r1, r2, r3)
    print("best fitting ellipsoid radii:", radii)
    dist = None
    # dist = lambda x: x**2 + 1/6
    # dist = lambda x: 1

    # use pyvista to get a mesh for the ellipsoid
    perfect_ellipsoid = pv.ParametricEllipsoid(*radii)
    crest_indices, p_crest_curve = get_crest_points_of_perfect_ellipsoid(radii, perfect_ellipsoid)
    # crest_indices, p_crest_curve = None, None

    # use Srep class to fit an srep to the ellipsoid and get vertices and centers of curvature
    srep = Srep()
    srep_pd, coc_indices, vertex_indices = srep.fit_ellipsoid_radii(radii, num_spine_points, num_steps=num_skeletal_onionskins)

    # p=pv.Plotter()
    # p.add_mesh(srep_pd, color='blue')
    # p.add_mesh(ellipsoid, opacity=0.2, color="white")
    # p.add_points(coc_points, render_points_as_spheres=True, point_size=10.0, color='white')
    # p.add_points(vertex_points, render_points_as_spheres=True, point_size=10.0, color='red')
    # p.show()

    # rotate srep back so that it fits original mesh
    srep_rotated_back = pv.wrap(vtk.vtkPolyData())
    srep_rotated_back.DeepCopy(srep_pd)
    srep_rotated_back.points = (srep_rotated_back.points @ np.linalg.inv(v)) + center
    srep_pd = srep_rotated_back

    perfect_ellipsoid_rotated_back = pv.wrap(vtk.vtkPolyData())
    perfect_ellipsoid_rotated_back.DeepCopy(perfect_ellipsoid)
    perfect_ellipsoid_rotated_back.points = (perfect_ellipsoid_rotated_back.points @ np.linalg.inv(v)) + center
    perfect_ellipsoid = perfect_ellipsoid_rotated_back

    # compute coc and vertices from indices
    coc_points = [srep_pd.GetPoint(coc_indices[0]), srep_pd.GetPoint(coc_indices[1])]
    coc_points = sorted(coc_points, key=lambda x: x[0]) # sort from low to high x coordinate
    coc_points = np.array(coc_points)
    vertex_points = [srep_pd.GetPoint(vertex_indices[0]), srep_pd.GetPoint(vertex_indices[1])]
    vertex_points = sorted(vertex_points, key=lambda x: x[0]) # sort from low to high x coordinate
    vertex_points = np.array(vertex_points)

    # p_crest_curve = (p_crest_curve @ np.linalg.inv(v)) + center

    return perfect_ellipsoid, srep_pd, coc_points, vertex_points, crest_indices, p_crest_curve
