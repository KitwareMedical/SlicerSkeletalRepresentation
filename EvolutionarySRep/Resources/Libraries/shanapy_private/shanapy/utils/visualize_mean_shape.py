import vtk

top_reader = vtk.vtkPolyDataReader()
top_reader.SetFileName('../data_diff_mean/final_mesh/top_mean_label_SPHARM.vtk')
top_reader.Update()
top_vtk = top_reader.GetOutput()

bot_reader = vtk.vtkPolyDataReader()
bot_reader.SetFileName('../data_diff_mean/final_mesh/bot_mean_label_SPHARM.vtk')
bot_reader.Update()
bot_vtk = bot_reader.GetOutput()

transformer = vtk.vtkTransform()
transformer.Translate((0, 0, -50)) ## Note: consider the distance in the image space
transformer.Update()
transform_filter = vtk.vtkTransformPolyDataFilter()
transform_filter.SetInputData(top_vtk)
transform_filter.SetTransform(transformer)
transform_filter.Update()
trans_std_ell = transform_filter.GetOutput()

dist_filter = vtk.vtkDistancePolyDataFilter()
dist_filter.SetInputData(0, bot_vtk)
dist_filter.SetInputData(1, trans_std_ell)
dist_filter.Update()

top_dist_filter = vtk.vtkDistancePolyDataFilter()
top_dist_filter.SetInputData(0, top_vtk)
top_dist_filter.SetInputData(1, top_vtk)
top_dist_filter.Update()

append_filter = vtk.vtkAppendPolyData()
append_filter.AddInputData(dist_filter.GetOutput())
append_filter.AddInputData(top_dist_filter.GetOutput())
append_filter.Update()
writer = vtk.vtkPolyDataWriter()
writer.SetFileName('../data_diff_mean/final_mesh/heat_maps_top_mean.vtk')
writer.SetInputData(top_dist_filter.GetOutput())
writer.Update()

writer = vtk.vtkPolyDataWriter()
writer.SetFileName('../data_diff_mean/final_mesh/heat_maps_bot_mean.vtk')
writer.SetInputData(dist_filter.GetOutput())
writer.Update()
print('done')
# mapper = vtk.vtkPolyDataMapper()
# mapper.SetInputData(append_filter.GetOutput())
# actor = vtk.vtkActor()
# actor.SetMapper(mapper)
# ren1 = vtk.vtkRenderer()
# ren1.AddActor(actor)

# # Assign our actor to the renderer.
# colors = vtk.vtkNamedColors()
# ren1.SetBackground(colors.GetColor3d("Beige"))
# ren1.GetActiveCamera().ParallelProjectionOn()
# renWin = vtk.vtkRenderWindow()
# renWin.AddRenderer(ren1)

# renWin.SetSize(640, 480)
# iren = vtk.vtkRenderWindowInteractor()
# iren.SetRenderWindow(renWin)

# style = vtk.vtkInteractorStyleTrackballCamera()
# iren.SetInteractorStyle(style)
# renWin.GetInteractor().SetInteractorStyle(style)

# ren1.ResetCamera()
# renWin.Render()

# iren.Initialize()
# iren.Start()
