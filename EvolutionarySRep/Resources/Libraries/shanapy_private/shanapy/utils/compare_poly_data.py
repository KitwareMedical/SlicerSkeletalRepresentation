import vtk
import numpy as np
def compare_polydata(poly1, poly2):
    """
    Compare if two surfaces/s-reps are the same.

    Return True if same
    """
    num_pts = poly1.GetNumberOfPoints()
    if num_pts != poly2.GetNumberOfPoints():
        return False

    for i in range(num_pts):
        pt1 = np.array(poly1.GetPoint(i))
        pt2 = np.array(poly2.GetPoint(i))

        norm = np.linalg.norm(pt1 - pt2)
        if norm > 0.1:
            return False

    return True

if __name__ == '__main__':
    case_id = np.random.randint(60)
    reader1 = vtk.vtkPolyDataReader()
    reader1.SetFileName('../data/final_mesh/bot' + str(case_id) + '_label_SPHARM.vtk')
    reader1.Update()

    reader2 = vtk.vtkPolyDataReader()
    reader2.SetFileName('../data/final_mesh/top' + str(case_id) + '_label_SPHARM.vtk')
    reader2.Update()

    print(case_id)
    print(compare_polydata(reader1.GetOutput(), reader2.GetOutput()))