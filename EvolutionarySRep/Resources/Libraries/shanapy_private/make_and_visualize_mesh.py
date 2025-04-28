import pyvista as pv
import vtk
import numpy as np

points = np.random.rand(1000,3)

vp = vtk.vtkPoints()
vs = vtk.vtkCellArray()
pd = vtk.vtkPolyData()

point_ids = []
for i in range(len(points)):
    pid = vp.InsertNextPoint(points[i])
    point_ids.append(i)

for j in range(len(points) - 1):
    spoke = vtk.vtkLine()
    spoke.GetPointIds().SetId(0, point_ids[j])
    spoke.GetPointIds().SetId(1, point_ids[j+1])
    vs.InsertNextCell(spoke)

pd.SetPoints(vp)
pd.SetLines(vs)

plt = pv.Plotter()
plt.add_mesh(pd)
# plt.show(screenshot="make_and_visualize_mesh.png")
plt.show()

    

    