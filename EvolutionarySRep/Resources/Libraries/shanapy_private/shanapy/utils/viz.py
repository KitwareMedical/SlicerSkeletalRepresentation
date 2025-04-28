import vtk
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from sklearn.decomposition import PCA
import plotly.express as px
from dwd import DWD, KernGDWD
from sklearn import metrics
def vtk_show(renderer, w=1280, h=900):
    """
    Takes vtkRenderer instance and returns an IPython Image with the rendering.
    """
    renderWindow = vtk.vtkRenderWindow()
#     renderWindow.SetOffScreenRendering(1)
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(w, h)
    renderWindow.Render()
     
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.Update()
     
    writer = vtk.vtkPNGWriter()
    writer.SetWriteToMemory(1)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()
    data = bytes(memoryview(writer.GetResult()))
    
    from IPython.display import Image
    return Image(data)

def viz_joint_variation_on_sphere(X, Y, x_joint, y_joint=None, x_joint_gt=None, title=""):
    """
    Draw directional data x, y and their joint components x_joint, y_joint respectively

    X: n x 3
    x_joint_gt: ground truth of the joint component of x
    """
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = 1 * np.outer(np.cos(u), np.sin(v))
    y = 1 * np.outer(np.sin(u), np.sin(v))
    z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
    #for i in range(2):
    #    ax.plot_surface(x+random.randint(-5,5), y+random.randint(-5,5), z+random.randint(-5,5),  rstride=4, cstride=4, color='b', linewidth=0, alpha=0.5)
    elev = 10.0
    rot = 80.0 / 180 * np.pi
    ax.plot_surface(x, y, z,  rstride=4, cstride=4, cmap=cm.Greys, linewidth=0, alpha=0.2)
    #calculate vectors for "vertical" circle
    a = np.array([-np.sin(elev / 180 * np.pi), 0, np.cos(elev / 180 * np.pi)])
    b = np.array([0, 1, 0])
    b = b * np.cos(rot) + np.cross(a, b) * np.sin(rot) + a * np.dot(a, b) * (1 - np.cos(rot))
    ax.plot(np.sin(u),np.cos(u),0,color='k', linestyle = 'dashed')
    horiz_front = np.linspace(0, np.pi, 100)
    ax.plot(np.sin(horiz_front),np.cos(horiz_front),0,color='k')
    # vert_front = np.linspace(np.pi / 2, 3 * np.pi / 2, 100)
    # ax.plot(a[0] * np.sin(u) + b[0] * np.cos(u), b[1] * np.cos(u), a[2] * np.sin(u) + b[2] * np.cos(u),color='k', linestyle = 'dashed')
    # ax.plot(a[0] * np.sin(vert_front) + b[0] * np.cos(vert_front), b[1] * np.cos(vert_front), a[2] * np.sin(vert_front) + b[2] * np.cos(vert_front),color='k')

    ax.view_init(elev = elev, azim = 0)

    ax.scatter(X[:, 0], X[:, 1], X[:, 2], marker='o', label='X')
    ax.scatter(Y[:, 0], Y[:, 1], Y[:, 2], marker='v', label='Y')
    ax.scatter(x_joint[:, 0], x_joint[:, 1], x_joint[:, 2], marker='*', label='X_joint_AJIVE')
    ax.scatter(y_joint[:, 0], y_joint[:, 1], y_joint[:, 2], marker='*', label='Y_joint_AJIVE')
    if x_joint_gt is not None:
        ax.scatter(x_joint_gt[:, 0], x_joint_gt[:, 1], x_joint_gt[:, 2], marker='*', label='X_joint_gt')
    # for i in range(X.shape[0]):
    #     px, py, pz = X[i, :]

    #     px2, py2, pz2 = Y[i, :]
    #     pt_x_joint = x_joint[i, :]
    #     pt_y_joint = y_joint[i, :]
    #     ax.scatter(px, py, pz, c='r', label='X')
    #     ax.scatter(px2, py2, pz2, c='b', marker='v', label='Y')
    #     # ax.plot(pt_x_joint[0], pt_x_joint[1], pt_x_joint[2], c='#17becf', marker='*-*', label='X_joint')
    #     # ax.plot(pt_y_joint[0], pt_y_joint[1], pt_y_joint[2], c='tab:brown', marker='*-*')


    plt.axis('off')
    plt.legend()
    plt.title(title)
    plt.show()
def form_spokes_poly(list_spokes):
    """
    Convert a list of Spoke to polydata
    """
    ret_poly = vtk.vtkPolyData()
    ret_pts = vtk.vtkPoints()
    ret_line = vtk.vtkCellArray()
    for i in range(len(list_spokes)):
        id0 = ret_pts.InsertNextPoint(list_spokes[i].p)
        id1 = ret_pts.InsertNextPoint(list_spokes[i].getB())
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, id0)
        line.GetPointIds().SetId(1, id1)
        ret_line.InsertNextCell(line)
    ret_poly.SetPoints(ret_pts)
    ret_poly.SetLines(ret_line)
    ret_poly.Modified()
    return ret_poly
def form_spoke_poly(np_base_pt, np_bdry_pt):
    """
    Convert a spoke to polydata for visualization
    """
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
def form_strata_from_spokes(deformed_srep):
    bdry_poly = vtk.vtkPolyData()
    bdry_pts = vtk.vtkPoints()
    skeletal_poly = vtk.vtkPolyData()
    skeletal_pts = vtk.vtkPoints()
    for j in range(deformed_srep.GetNumberOfCells()):
        bdry_pt_id = j * 2
        bdry_pts.InsertNextPoint(deformed_srep.GetPoint(bdry_pt_id))
        skeletal_pts.InsertNextPoint(deformed_srep.GetPoint(bdry_pt_id + 1))
    bdry_poly.SetPoints(bdry_pts)
    skeletal_poly.SetPoints(skeletal_pts)
    bdry_poly.Modified()
    skeletal_poly.Modified()
    return bdry_poly, skeletal_poly

def form_spokes_from_strata(skeleton, wave_front_surf):
    """
    Connect boundary surface and skeletal sheet to form polydata for spokes
    """
    connection_poly = vtk.vtkPolyData()
    tail_head_pt = vtk.vtkPoints()
    connection_line = vtk.vtkCellArray()
    for i in range(skeleton.GetNumberOfPoints()):
        id_tail = tail_head_pt.InsertNextPoint(skeleton.GetPoint(i))
        id_head = tail_head_pt.InsertNextPoint(wave_front_surf.GetPoint(i))

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, id_tail)
        line.GetPointIds().SetId(1, id_head)
        connection_line.InsertNextCell(line)

    connection_poly.SetPoints(tail_head_pt)
    connection_poly.SetLines(connection_line)
    connection_poly.Modified()
    return connection_poly
    
def viz_directions_distribution(pts, pts2, mu1=None, mu2=None, pts_joint=None, title='', show_now=True):
    """
    Draw directional data distributed on a unit sphere
    Input pts: n x 3
    Input pts2: nx 3
    Input pts_joint: the bending directions of pts n x 3
    Input mu1/mu2: 1 x 3
    """
    assert pts.shape[0] == pts2.shape[0], "Two blocks have to be of same number of samples"
    ### draw unit sphere
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection='3d')
    #ax.set_aspect('equal')

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = 1 * np.outer(np.cos(u), np.sin(v))
    y = 1 * np.outer(np.sin(u), np.sin(v))
    z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
    #for i in range(2):
    #    ax.plot_surface(x+random.randint(-5,5), y+random.randint(-5,5), z+random.randint(-5,5),  rstride=4, cstride=4, color='b', linewidth=0, alpha=0.5)
    elev = 10.0
    rot = 80.0 / 180 * np.pi
    ax.plot_surface(x, y, z,  rstride=4, cstride=4, cmap=cm.Greys, linewidth=0, alpha=0.3)
    #calculate vectors for "vertical" circle
    a = np.array([-np.sin(elev / 180 * np.pi), 0, np.cos(elev / 180 * np.pi)])
    b = np.array([0, 1, 0])
    b = b * np.cos(rot) + np.cross(a, b) * np.sin(rot) + a * np.dot(a, b) * (1 - np.cos(rot))
    ax.plot(np.sin(u),np.cos(u),0,color='k', linestyle = 'dashed')
    horiz_front = np.linspace(0, np.pi, 100)
    ax.plot(np.sin(horiz_front),np.cos(horiz_front),0,color='k')
    # vert_front = np.linspace(np.pi / 2, 3 * np.pi / 2, 100)
    # ax.plot(a[0] * np.sin(u) + b[0] * np.cos(u), b[1] * np.cos(u), a[2] * np.sin(u) + b[2] * np.cos(u),color='k', linestyle = 'dashed')
    # ax.plot(a[0] * np.sin(vert_front) + b[0] * np.cos(vert_front), b[1] * np.cos(vert_front), a[2] * np.sin(vert_front) + b[2] * np.cos(vert_front),color='k')

    ax.view_init(elev = elev, azim = 0)

    for i in range(pts.shape[0]):
        px, py, pz = pts[i, :]

        px2, py2, pz2 = pts2[i, :]


        ax.scatter(px, py, pz, c='r')
        ax.scatter(px2, py2, pz2, c='b', marker='v')
        if pts_joint is not None:
            px_j, py_j, pz_j = pts_joint[i, :]
            ax.scatter(px_j, py_j, pz_j, c='y', marker='*')

    if mu1 is not None:
        ax.scatter(mu1[0], mu1[1], mu1[2], c='r', marker='x')
    if mu2 is not None:
        ax.scatter(mu2[0], mu2[1], mu2[2], c='b', marker='x')
    plt.axis('off')
    plt.title(title)
    if show_now:
        plt.show()

def compare_first_spokes(poly_a, poly_b, background=None):
    pt_a_0 = np.array(poly_a.GetPoint(0))
    pt_a_1 = np.array(poly_a.GetPoint(1))
    spoke_a = form_spoke_poly(pt_a_0, pt_a_1)

    pt_b_0 = np.array(poly_b.GetPoint(0))
    pt_b_1 = np.array(poly_b.GetPoint(1))
    spoke_b = form_spoke_poly(pt_b_0, pt_b_1)

    overlay_polydata(spoke_a, background, highlight=spoke_b)
def compare_last_spokes(poly_a, poly_b, background=None):
    pt_a_0 = np.array(poly_a.GetPoint(poly_a.GetNumberOfPoints() - 4))
    pt_a_1 = np.array(poly_a.GetPoint(poly_a.GetNumberOfPoints() - 3))
    spoke_a = form_spoke_poly(pt_a_0, pt_a_1)

    pt_b_0 = np.array(poly_b.GetPoint(poly_b.GetNumberOfPoints() - 4))
    pt_b_1 = np.array(poly_b.GetPoint(poly_b.GetNumberOfPoints() - 3))
    spoke_b = form_spoke_poly(pt_b_0, pt_b_1)

    overlay_polydata(spoke_a, background, highlight=spoke_b)
def close_window(iren):
    render_window = iren.GetRenderWindow()
    render_window.Finalize()
    iren.TerminateApp()
def overlay_polydata(foreground = None, background = None, highlight = None, appendum=None, appendum2=None, inline=False, flash=False, highlight_color=(255, 0, 0)):
    """
    Show two vtk polydata together
    """
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(foreground)
    colors = vtk.vtkNamedColors()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(5)
    actor.GetProperty().SetColor(colors.GetColor3d("Blue"))
    actor.GetProperty().SetLineWidth(2)

    bk_mapper = vtk.vtkPolyDataMapper()
    bk_mapper.SetInputData(background)
    bk_actor = vtk.vtkActor()
    bk_actor.SetMapper(bk_mapper)
    bk_actor.GetProperty().SetColor(colors.GetColor3d("dim_grey"))
    bk_actor.GetProperty().SetLineWidth(1)
    bk_actor.GetProperty().SetOpacity(0.2)
#    bk_actor.GetProperty().EdgeVisibilityOn()

    ren1 = vtk.vtkRenderer()
    ren1.AddActor(actor)
    ren1.AddActor(bk_actor)

    if highlight is not None:
        hl_mapper = vtk.vtkPolyDataMapper()
        hl_mapper.SetInputData(highlight)
        hl_actor = vtk.vtkActor()
        hl_actor.SetMapper(hl_mapper)
        hl_actor.GetProperty().SetColor(highlight_color)#colors.GetColor3d("Red"))
        hl_actor.GetProperty().SetLineWidth(2)

        ren1.AddActor(hl_actor)
    if appendum is not None:
        appendum_mapper = vtk.vtkPolyDataMapper()
        appendum_mapper.SetInputData(appendum)
        appendum_actor = vtk.vtkActor()
        appendum_actor.SetMapper(appendum_mapper)
        appendum_actor.GetProperty().SetColor(colors.GetColor3d("Blue"))
        appendum_actor.GetProperty().SetLineWidth(1)

        ren1.AddActor(appendum_actor)
    if appendum2 is not None:
        appendum2_mapper = vtk.vtkPolyDataMapper()
        appendum2_mapper.SetInputData(appendum2)
        appendum2_actor = vtk.vtkActor()
        appendum2_actor.SetMapper(appendum2_mapper)
        appendum2_actor.GetProperty().SetColor(colors.GetColor3d("banana"))
        appendum2_actor.GetProperty().SetLineWidth(1)

        ren1.AddActor(appendum2_actor)

    ren1.SetBackground(colors.GetColor3d("ivory"))
    ren1.GetActiveCamera().ParallelProjectionOn()

    ren1.ResetCamera()
    if inline:
        return ren1
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren1)
    renWin.SetSize(1280, 900)
    style = vtk.vtkInteractorStyleTrackballCamera()

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetInteractorStyle(style)
    iren.SetRenderWindow(renWin)

    renWin.GetInteractor().SetInteractorStyle(style)
    renWin.Render()

    iren.Initialize()
    iren.Start()

    if flash:
        print('Closing window')
        close_window(iren)
        del renWin, iren

def visualize(polydata, background=None):
    colors = vtk.vtkNamedColors()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    colors = vtk.vtkNamedColors()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(5)
    actor.GetProperty().SetColor(colors.GetColor3d("dim_grey"))
    actor.GetProperty().SetLineWidth(5)
    ren1 = vtk.vtkRenderer()
    ren1.AddActor(actor)

    if background is not None:
        bk_mapper = vtk.vtkPolyDataMapper()
        bk_mapper.SetInputData(background)
        bk_colors = vtk.vtkNamedColors()
        bk_actor = vtk.vtkActor()
        bk_actor.SetMapper(bk_mapper)
        bk_actor.GetProperty().SetColor(bk_colors.GetColor3d("Tomato"))
        bk_actor.GetProperty().SetLineWidth(1)
        bk_actor.GetProperty().SetOpacity(0.2)
        ren1.AddActor(bk_actor)
    ren1.SetBackground(colors.GetColor3d("ivory"))

    ren1.GetActiveCamera().ParallelProjectionOn()

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren1)
    renWin.SetSize(1280, 900)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    renWin.GetInteractor().SetInteractorStyle(style)

    ren1.ResetCamera()
    renWin.Render()

    iren.Initialize()
    iren.Start()
def pca_explained_ratio_plots(x):
    pcamodel = PCA(n_components=5)
    pca = pcamodel.fit_transform(x)
    plt.bar(range(1,len(pcamodel.explained_variance_ratio_ )+1),pcamodel.explained_variance_ratio_ )
    plt.ylabel('Explained variance')
    plt.xlabel('Components')
    plt.plot(range(1,len(pcamodel.explained_variance_ratio_ )+1),
             np.cumsum(pcamodel.explained_variance_ratio_),
             c='red',
             label="Cumulative Explained Variance")
    plt.legend(loc='upper right')
    plt.show()
def draw_compactness(spoke_feats, link_feats):
    pca_model = PCA(n_components=100)
    pca_spoke = pca_model.fit(spoke_feats)

    pca_link = pca_model.fit(link_feats)
    k_feat_spoke = 1.0/spoke_feats.shape[1]
    k_link_spoke = 1.0/link_feats.shape[1]
    plt.plot(range(1,len(pca_spoke.explained_variance_ )+1),
             np.cumsum(pca_spoke.explained_variance_),
             c='red',
             label="Concate s-reps")
    plt.plot(range(1,len(pca_link.explained_variance_ )+1),
             np.cumsum(pca_link.explained_variance_),
             c='blue',
                 label="Link")
    plt.legend()
    plt.title('Compactness')
    plt.show()
def draw_bar_graph(spoke_ent, link_ent, link_ext_ent):
    labels = ['Real data', 'Simulated data']

    x = np.arange(len(labels))  # the label locations
    width = 0.3  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width, spoke_ent, width, label='Concate s-reps')
    rects2 = ax.bar(x, link_ent, width, label='Link')
    rects3 = ax.bar(x + width, link_ext_ent, width, label='w/ link extension*' )

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Entropy')
    ax.set_title('Entropy of 3 models')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    plt.show()
def pca_scatter_plot(features, class_labels):
    """
    features: n x d feature matrix
    """
    print(features.shape)
    
    from sklearn.decomposition import PCA
    pca = PCA()
    components = pca.fit_transform(features)
    labels = {
        str(i): f"PC {i+1} ({var:.1f}%)"
        for i, var in enumerate(pca.explained_variance_ratio_ * 100)
    }

    fig = px.scatter_matrix(
        components,
        labels=labels,
        dimensions=range(4),
        color=class_labels
    )
    fig.update_traces(diagonal_visible=True)
    fig.show()
def scatter_matrix_plot(components, class_labels, dimensions=range(4)):
    fig = px.scatter_matrix(
        components,
        dimensions=dimensions,
        color=class_labels
    )
    fig.update_traces(diagonal_visible=True)
    fig.show()
def scatter_pts(pts, pts2):
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2])
    ax.scatter(pts2[:, 0], pts2[:, 1], pts2[:, 2])
    plt.show()
def project_to_decision_boundary(X, class_labels):
    """
    Input X matrix of dimension n x d
    """
    dwd = DWD().fit(X, class_labels)
    sep_dir = dwd.coef_
    x = np.dot(X, sep_dir.T)

    remains = X - x
    pca = PCA(n_components=2)
    pca.fit(remains)
    ortho_dir = pca.components_[0, :][None, :]
    y = np.dot(X, ortho_dir.T)

    plt.scatter(x[class_labels==1, :], y[class_labels==1, :], label='Positive', c='red')
    plt.scatter(x[class_labels==0, :], y[class_labels==0, :], label='Negative', c='blue')
    plt.xlabel('DWD direction')
    plt.tick_params(labelleft=False)
    plt.legend()

    plt.show()
    coef =dwd.coef_
    plt.plot(coef.squeeze())
    plt.show()
    y_pred = dwd.decision_function(X)
    fpr, tpr, thresholds = metrics.roc_curve(class_labels, y_pred, pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    print("AUC:", roc_auc)
    return dwd
